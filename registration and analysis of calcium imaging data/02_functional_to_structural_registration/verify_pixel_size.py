#!/usr/bin/env python3
"""
verify_pixel_size.py
======================
Stage 2 of 6 (functional-to-structural registration) -- Step 0.

Standalone verification utility. Re-derives the functional 2P camera's true
physical pixel size directly from the OME-XML metadata embedded in the raw
Olympus .oir acquisition files, rather than trusting a hardcoded comment.

Why this exists: the "1.40625 um/px" figure used throughout this stage (see
build_functional_ministack.py) was originally just an unsourced comment --
stage 01's convert_oir_to_tiff.py never actually captured PhysicalSizeX/Y from
the OME metadata, only sizeX/Y/Z/T/C, bitsPerPixel, pixelType, and
dimensionOrder. This project independently hit costly mm-vs-micron unit bugs
in other parts of the pipeline (a mask-building step and a template-geometry
step each assumed the wrong linear unit at one point), so as a general rule
here, any pixel-size number used for registration geometry gets re-derived
from source metadata rather than trusted from a comment. This script re-opens
raw .oir files and reads PhysicalSizeX/Y/Z directly from the OME-XML metadata
store, independent of any value quoted elsewhere in the codebase. The value
confirmed here (1.40625 um/px in X and Y, 10.0 um between planes in Z) is what
build_functional_ministack.py uses for the NIfTI voxel spacing.

Reuses the same Bio-Formats/jpype bridge as stage 01's convert_oir_to_tiff.py
(same bundled bioformats_package.jar, same JVM) since that is the
already-proven-working Bio-Formats setup in this project.

Inputs (not part of data_for_upload -- these are large raw acquisition files
and a local Java runtime, not redistributed data):
  - the raw .oir functional recordings (point FUNCTIONAL_OIR_DIR at a local copy)
  - a local bioformats_package.jar (BIOFORMATS_JAR_PATH)
  - a local JVM shared library, e.g. jvm.dll/libjvm.so (BIOFORMATS_JVM_PATH)

Output: none written to disk. Prints the PhysicalSizeX/Y/Z found in each file
so it can be compared by eye against the constants used downstream.

Environment: same as stage 01 (a "suite2p" env with jpype installed).

Run:
    python verify_pixel_size.py
with FUNCTIONAL_OIR_DIR, BIOFORMATS_JAR_PATH, BIOFORMATS_JVM_PATH set in the
environment.
"""

import glob
import os

import jpype

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

# Raw .oir acquisition files and the Bio-Formats JVM/jar are local, machine-specific
# resources that live outside data_for_upload -- configure via environment variables.
OIR_DIR = os.environ.get("FUNCTIONAL_OIR_DIR")
JAR_PATH = os.environ.get("BIOFORMATS_JAR_PATH")
JVM_PATH = os.environ.get("BIOFORMATS_JVM_PATH")


def main():
    if not OIR_DIR or not JAR_PATH or not JVM_PATH:
        raise SystemExit(
            "Set FUNCTIONAL_OIR_DIR (folder of raw .oir files), BIOFORMATS_JAR_PATH "
            "(path to bioformats_package.jar), and BIOFORMATS_JVM_PATH (path to a "
            "JVM shared library) in the environment before running this script."
        )

    if not jpype.isJVMStarted():
        jpype.startJVM(JVM_PATH, classpath=[JAR_PATH], convertStrings=False)

    ImageReader = jpype.JClass("loci.formats.ImageReader")
    MetadataTools = jpype.JClass("loci.formats.MetadataTools")
    DebugTools = jpype.JClass("loci.common.DebugTools")
    UNITS = jpype.JClass("ome.units.UNITS")
    DebugTools.enableLogging("ERROR")

    files = sorted(glob.glob(os.path.join(OIR_DIR, "*.oir")))
    print(f"Found {len(files)} .oir files; checking physical pixel size on each\n")

    reader = ImageReader()
    results = []
    for path in files:
        name = os.path.basename(path)
        ome_meta = MetadataTools.createOMEXMLMetadata()
        reader.setMetadataStore(ome_meta)
        reader.setId(path)

        px = ome_meta.getPixelsPhysicalSizeX(0)
        py = ome_meta.getPixelsPhysicalSizeY(0)
        pz = ome_meta.getPixelsPhysicalSizeZ(0)

        x_um = float(px.value(UNITS.MICROMETER).doubleValue()) if px is not None else None
        y_um = float(py.value(UNITS.MICROMETER).doubleValue()) if py is not None else None
        z_um = float(pz.value(UNITS.MICROMETER).doubleValue()) if pz is not None else None

        print(f"{name}: PhysicalSizeX={x_um} um  PhysicalSizeY={y_um} um  PhysicalSizeZ={z_um} um")
        results.append((name, x_um, y_um, z_um))
        reader.close()

    xs = {r[1] for r in results}
    ys = {r[2] for r in results}
    print(f"\nDistinct PhysicalSizeX values across all files: {xs}")
    print(f"Distinct PhysicalSizeY values across all files: {ys}")
    print("\nCompare against the XY_SPACING_UM / Z_SPACING_UM constants used in "
          "build_functional_ministack.py (1.40625 um and 10.0 um).")


if __name__ == "__main__":
    main()
