#!/usr/bin/env python3
"""
This script reads a geo file and a PMT pos_table file (e.g. "DOMINFO_5str"),
clusters every 31 PMTs into a single cluster by averaging their positions,
and then outputs a GDML file with each cluster represented as a sphere with a
hard-coded outer radius of 400 mm. The volume hierarchy is defined from inner
(child) to outer (mother) volumes.

Usage:
    python3 geo_to_gdml.py mygeo.txt mypmtinfo.txt
"""

import sys
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom

# --- Helper Functions ---

def remove_comments(text):
    """Remove C/C++ style comments from the input text."""
    # Remove // comments
    return re.sub(r'//.*', '', text)

def parse_custom_block(text):
    """
    Very simple parser to convert a block with syntax like:
      { key: value, key: value, ... }
    into a Python dict.
    It assumes that strings (for keys and string values) are optionally quoted.
    """
    text = remove_comments(text)
    # Remove braces and newlines
    text = text.strip().lstrip("{").rstrip("}")
    # Split by commas that are not inside brackets
    tokens = re.split(r',(?![^\[]*\])', text)
    result = {}
    for token in tokens:
        if token.strip() == "":
            continue
        # Split on first colon
        parts = token.split(":", 1)
        if len(parts) != 2:
            continue
        key = parts[0].strip().strip('"')
        value = parts[1].strip()
        # If value starts with a bracket, try to parse as list of numbers
        if value.startswith("[") and value.endswith("]"):
            # Remove brackets and split by commas
            inner = value[1:-1]
            # Convert to float if possible
            values = []
            for num in inner.split(","):
                try:
                    values.append(float(num))
                except ValueError:
                    values.append(num.strip().strip('"'))
            result[key] = values
        else:
            # Try to remove quotes from value
            value = value.strip().strip('"')
            # Try to convert to float/int if possible
            try:
                if "." in value:
                    result[key] = float(value)
                else:
                    result[key] = int(value)
            except ValueError:
                result[key] = value
    return result

def parse_file(filename):
    """
    Parses the input file that is expected to contain one or more blocks
    (delimited by curly braces). Returns a list of dictionaries.
    """
    with open(filename, "r") as f:
        content = f.read()
    # Split the content on blocks – each block starts with '{' and ends with '}'
    blocks = re.findall(r"\{[^}]*\}", content, re.DOTALL)
    parsed = [parse_custom_block(block) for block in blocks]
    return parsed

def average(lst):
    """Compute the average of a list."""
    return sum(lst)/len(lst) if lst else 0.0

def prettify(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

# --- GDML Generation Functions ---

def create_gdml(geo_blocks, pmt_tables):
    """
    Create a GDML structure as an XML tree.
    Assumes that the pmt_tables is a list of dicts from which we pick the pos_table
    named (for example) "DOMINFO_5str".
    """
    # Find the pos_table we need
    pos_table = None
    for block in pmt_tables:
        if block.get("name", "") in ["DOMINFO_5str", "PMTINFO_bottom"]:
            pos_table = block
            break
    if pos_table is None:
        print("Error: Could not find a pos_table named 'DOMINFO_5str' or 'PMTINFO_bottom'")
        sys.exit(1)

    # Extract PMT positions from pos_table: expect keys "x", "y", "z"
    xs = pos_table.get("x", [])
    ys = pos_table.get("y", [])
    zs = pos_table.get("z", [])
    num_pmts = len(xs)
    if not (len(ys)==num_pmts and len(zs)==num_pmts):
        print("Error: Mismatched PMT table dimensions")
        sys.exit(1)

    # Cluster every 31 PMTs into one cluster.
    clusters = []
    for i in range(0, num_pmts, 31):
        chunk_x = xs[i:i+31]
        chunk_y = ys[i:i+31]
        chunk_z = zs[i:i+31]
        # Only add a cluster if there is at least one PMT
        if chunk_x:
            clusters.append( (average(chunk_x), average(chunk_y), average(chunk_z)) )

    # --- Build GDML XML ---
    gdml = ET.Element("gdml")

    # Define materials (example: water)
    materials = ET.SubElement(gdml, "materials")
    water = ET.SubElement(materials, "material", name="water")
    ET.SubElement(water, "D").text = "1.0"  # dummy density

    # Define solids
    solids = ET.SubElement(gdml, "solids")

    # Create a sphere solid for PMT clusters: always hard-coded to a sphere of 400 mm radius.
    # GDML sphere size: [rmin, rmax, startPhi, deltaPhi, startTheta, deltaTheta]
    cluster_solid = ET.SubElement(solids, "sphere", name="clusterSphere", lunit="mm",
                                  rmin="0.0", rmax="400.0",
                                  startphi="0.0", deltaphi="360.0",
                                  starttheta="0.0", deltatheta="180.0")

    # (If your geo file defines additional solids, you could parse them here)

    # Define volumes. We define daughter volumes first.
    # 1. Define each PMT cluster volume (child volumes)
    volumes = ET.SubElement(gdml, "structure")
    cluster_volumes = []
    for idx, (x, y, z) in enumerate(clusters):
        vol_name = f"cluster_{idx}"
        volume = ET.SubElement(volumes, "volume", name=vol_name)
        ET.SubElement(volume, "materialref", ref="water")
        ET.SubElement(volume, "solidref", ref="clusterSphere")
        # Define the position of this cluster as a physical volume placement in its mother later.
        # Here we store the average positions as attributes in a custom tag.
        pos = ET.SubElement(volume, "auxiliary", auxtype="ClusterPosition", autype="XYZ")
        pos.text = f"{x} {y} {z}"
        cluster_volumes.append(vol_name)

    # 2. Define detector volume that will contain the clusters.
    # For detector, we use a box. You can change size as needed.
    detector = ET.SubElement(volumes, "volume", name="vol_detector")
    ET.SubElement(detector, "materialref", ref="water")
    # Define a box solid for the detector volume:
    det_box = ET.SubElement(solids, "box", name="detectorBox", lunit="mm",
                             x="20000", y="20000", z="200000")
    ET.SubElement(detector, "solidref", ref="detectorBox")

    # Place each cluster volume inside detector.
    for idx, vol_name in enumerate(cluster_volumes):
        # Get cluster position from our earlier stored averages.
        x, y, z = clusters[idx]
        physvol = ET.SubElement(detector, "physvol")
        ET.SubElement(physvol, "volumeref", ref=vol_name)
        # Create a translation element for the cluster placement.
        # Here we assume a translation; rotation is identity.
        pos_elem = ET.SubElement(physvol, "position", name=f"pos_{vol_name}",
                                 unit="mm", x=str(x), y=str(y), z=str(z))
        # For rotation, if needed:
        ET.SubElement(physvol, "rotation", name=f"rot_{vol_name}",
                      unit="deg", x="0", y="0", z="0")

    # 3. Define world volume (outermost) that contains the detector.
    world = ET.SubElement(volumes, "volume", name="vol_world")
    ET.SubElement(world, "materialref", ref="water")
    # Use a box solid for world:
    world_box = ET.SubElement(solids, "box", name="worldBox", lunit="mm",
                              x="50000", y="50000", z="250000")
    ET.SubElement(world, "solidref", ref="worldBox")
    # Place detector in world:
    world_physvol = ET.SubElement(world, "physvol")
    ET.SubElement(world_physvol, "volumeref", ref="vol_detector")
    ET.SubElement(world_physvol, "position", name="pos_detector",
                  unit="mm", x="0", y="0", z="0")
    ET.SubElement(world_physvol, "rotation", name="rot_detector",
                  unit="deg", x="0", y="0", z="0")

    # Add setup element referencing the world volume.
    setup = ET.SubElement(gdml, "setup", name="Default", version="1.0")
    ET.SubElement(setup, "world", ref="vol_world")

    return gdml

# --- Main Function ---

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 geo_to_gdml.py <geo_file> <pmt_info_file>")
        sys.exit(1)
    geo_filename = sys.argv[1]
    pmt_info_filename = sys.argv[2]

    # Parse input files
    try:
        geo_blocks = parse_file(geo_filename)
    except Exception as e:
        print("Error parsing geo file:", e)
        sys.exit(1)

    try:
        pmt_tables = parse_file(pmt_info_filename)
    except Exception as e:
        print("Error parsing pmt info file:", e)
        sys.exit(1)

    gdml_tree = create_gdml(geo_blocks, pmt_tables)

    gdml_str = prettify(gdml_tree)
    out_filename = "gdml_out.gdml"
    with open(out_filename, "w") as f:
        f.write(gdml_str)
    print("GDML output written to", out_filename)

if __name__ == "__main__":
    main()

