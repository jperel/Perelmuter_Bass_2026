#!/usr/bin/env python3
"""
convert_oir_to_tiff.py
========================
Stage 1 of 6 (suite2p preprocessing) -- step 1 of 3 in this folder.

Converts the 36 single-plane Olympus .oir functional recordings to TIFF stacks
that suite2p can read natively (suite2p.io has no OIR reader). Each recording
is a separate ~2-minute acquisition at a distinct native z-depth, not a
simultaneous multi-plane volume.

Reads:
    Raw .oir files from RAW_OIR_DIR (see below). These are the original
    acquisition files and are not part of the public data release (raw 2P
    movies are large); point RAW_OIR_DIR at your own copy of the recordings,
    or override it with the RAW_OIR_DIR environment variable.

Writes:
    One (T, Y, X) TIFF stack per input file into TIFF_DIR, same base name as
    the source file (just .tif instead of .oir), so the acquisition z-value
    stays attached to the filename for the next script in this folder.
    TIFF_DIR is also local/intermediate (not part of the public data release)
    and can be overridden with the FUNCTIONAL_TIFF_DIR environment variable.

Run order in this folder:
    1. convert_oir_to_tiff.py       (this script)
    2. run_suite2p_batch.py
    3. export_mean_images.py

Uses Bio-Formats (the standard, actively-maintained library with official
Olympus OIR support), driven via jpype against a self-contained
bioformats_package.jar (all Java dependencies bundled in one jar, which
avoids the multi-repository Maven dependency resolution that dynamic-download
approaches such as `bioio-bioformats` require and that can be unreliable in
restricted or offline environments).

Requires:
    - A JVM (JDK/JRE 8+). By default this script asks jpype to locate one
      (jpype.getDefaultJVMPath()), which requires a system Java installation
      or JAVA_HOME to be set. Override with the BIOFORMATS_JVM_PATH
      environment variable to point at a specific jvm.dll/libjvm.
    - bioformats_package.jar (the self-contained Bio-Formats distribution,
      available from the official Bio-Formats releases). Place it at
      bioformats_jar/bioformats_package.jar next to this script, or point at
      it with the BIOFORMATS_JAR_PATH environment variable.
    - Python packages: jpype1, numpy, tifffile.
"""

import glob
import os

import jpype
import numpy as np
import tifffile

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))
# DATA_ROOT is resolved here for consistency with the other pipeline scripts,
# but this script's own inputs/outputs are local (raw/intermediate data, not
# part of the public release) -- see RAW_OIR_DIR / TIFF_DIR below.

HERE = os.path.dirname(os.path.abspath(__file__))

# Raw .oir recordings and the intermediate per-plane TIFF stacks are large
# raw/intermediate acquisition data and are not included in the public data
# release (only the exported mean images in data_for_upload/ are). Both
# directories default to local folders next to this script and can be
# redirected via environment variables.
RAW_OIR_DIR = os.environ.get("RAW_OIR_DIR", os.path.join(HERE, "functional_z-planes_oir"))
TIFF_DIR = os.environ.get("FUNCTIONAL_TIFF_DIR", os.path.join(HERE, "functional_z-planes_tiff"))

JAR_PATH = os.environ.get(
    "BIOFORMATS_JAR_PATH", os.path.join(HERE, "bioformats_jar", "bioformats_package.jar"))
JVM_PATH = os.environ.get("BIOFORMATS_JVM_PATH") or jpype.getDefaultJVMPath()

# Bio-Formats FormatTools pixel type constants (loci.formats.FormatTools)
PIXEL_TYPE_TO_DTYPE = {
    0: (np.int8, 1), 1: (np.uint8, 1),
    2: (np.int16, 2), 3: (np.uint16, 2),
    4: (np.int32, 4), 5: (np.uint32, 4),
    6: (np.float32, 4), 7: (np.float64, 8),
}


def read_oir_as_stack(reader, path):
    """Returns (stack[T,Y,X], metadata dict)."""
    reader.setId(path)
    size_x, size_y = reader.getSizeX(), reader.getSizeY()
    size_z, size_t, size_c = reader.getSizeZ(), reader.getSizeT(), reader.getSizeC()
    pixel_type = reader.getPixelType()
    little_endian = bool(reader.isLittleEndian())
    dtype, nbytes = PIXEL_TYPE_TO_DTYPE[pixel_type]

    meta = {
        "sizeX": size_x, "sizeY": size_y, "sizeZ": size_z, "sizeT": size_t, "sizeC": size_c,
        "bitsPerPixel": reader.getBitsPerPixel(), "pixelType": pixel_type,
        "dimensionOrder": str(reader.getDimensionOrder()),
    }

    n_planes = reader.getImageCount()
    stack = np.empty((n_planes, size_y, size_x), dtype=dtype)
    byteorder = "<" if little_endian else ">"
    np_dtype = np.dtype(dtype).newbyteorder(byteorder)
    for i in range(n_planes):
        raw = bytes(reader.openBytes(i))
        plane = np.frombuffer(raw, dtype=np_dtype).reshape(size_y, size_x)
        stack[i] = plane.astype(dtype)  # normalize to native byte order in memory

    return stack, meta


def main():
    os.makedirs(TIFF_DIR, exist_ok=True)
    if not jpype.isJVMStarted():
        jpype.startJVM(JVM_PATH, classpath=[JAR_PATH], convertStrings=False)

    ImageReader = jpype.JClass("loci.formats.ImageReader")
    DebugTools = jpype.JClass("loci.common.DebugTools")
    DebugTools.enableLogging("ERROR")  # suppress verbose Bio-Formats debug-level logging

    files = sorted(glob.glob(os.path.join(RAW_OIR_DIR, "*.oir")))
    print(f"Found {len(files)} .oir files")

    reader = ImageReader()
    summary = []
    for path in files:
        name = os.path.splitext(os.path.basename(path))[0]
        out_path = os.path.join(TIFF_DIR, name + ".tif")
        stack, meta = read_oir_as_stack(reader, path)
        tifffile.imwrite(out_path, stack, imagej=True, metadata={"axes": "TYX"})
        print(f"{name}: shape={stack.shape} dtype={stack.dtype} meta={meta} -> {out_path}")
        summary.append((name, stack.shape, meta))
        reader.close()

    print("\n=== Summary ===")
    shapes = {s[1] for s in summary}
    print(f"Distinct (T,Y,X) shapes across all files: {shapes}")


if __name__ == "__main__":
    main()
