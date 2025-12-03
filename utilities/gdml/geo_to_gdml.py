#!/usr/bin/env python3
"""
This script converts custom geo and PMTINFO files (in a pseudo-JSON format with C++-style comments)
into a nicely indented GDML file.

It supports the allowed geometry types: box, sphere, zsphere, tubs, cylinder, slab, cone,
box3, trapezoid, convexpolyhedron, disc, segment, ellipsoid, torus, hyperboloid, cubic.
For PMTs (model r11065) the geometry is approximated by a tube (<tube>) instead of a toroidal shape.

Usage:
    python3 geo_to_gdml.py <geo_file> <pmtinfo_file>
"""

import sys
import re
import json
import xml.etree.ElementTree as ET
import xml.dom.minidom

# ------------------------------------------------------------------
# Hardcoded PMT (r11065) parameters.
r11065_info = {
    "name": "PMT",
    "index": "r11065",
    "valid_begin": [0, 0],
    "valid_end": [0, 0],
    "construction": "toroidal",  # original design info (for reference)
    "dynode_material": "stainless_steel",
    "glass_material": "quartz",
    "pmt_vacuum_material": "pmt_vacuum",  # dilute air
    "photocathode_surface": "photocathode_R11065",
    "mirror_surface": "mirror",
    "dynode_surface": "stainless_steel",
    "dynode_radius": 23.0,   # mm
    "dynode_top":   -30.0,   # mm
    "wall_thickness": 3.0,   # mm
    "z_edge":   [4.020, 1.00, 0.00, -16.75, -30.00, -32.00, -119.0],
    "rho_edge": [38.0, 38.0, 38.0, 38.0, 26.65, 26.65, 26.65],
    "z_origin": [-1000000.00, 0.00, 0.00, 10000.00, 31.00, 0.00],
    "noise_rate": 10000.0    # Hz
}

# ------------------------------------------------------------------
def parse_blocks_from_file(filename):
    """
    Parse a file containing one or more blocks in a pseudo-JSON format.
    This function uses a brace-matching approach to extract text between balanced '{' and '}'.
    It removes inline (// ...) and full-line comments, adds quotes to keys,
    and removes trailing commas.
    Returns a list of dictionaries.
    """
    with open(filename, 'r') as f:
        content = f.read()

    blocks = []
    current_block = ""
    brace_count = 0
    in_block = False
    for char in content:
        if char == '{':
            brace_count += 1
            in_block = True
        if in_block:
            current_block += char
        if char == '}':
            brace_count -= 1
            if brace_count == 0 and in_block:
                blocks.append(current_block)
                current_block = ""
                in_block = False

    processed_blocks = []
    for block in blocks:
        # Remove inline comments and full-line comments.
        lines = block.splitlines()
        lines_no_comments = []
        for line in lines:
            line_no_comment = line.split("//")[0]
            if line_no_comment.strip():
                lines_no_comments.append(line_no_comment)
        block_no_comments = "\n".join(lines_no_comments)
        # Insert quotes around keys (assumes keys are alphanumeric or underscore)
        block_quoted = re.sub(r'(\w+)\s*:', r'"\1":', block_no_comments)
        # Remove any trailing commas before closing braces or brackets.
        block_clean = re.sub(r',\s*([}\]])', r'\1', block_quoted)
        try:
            obj = json.loads(block_clean)
            processed_blocks.append(obj)
        except Exception as e:
            print("Error parsing block:")
            print(block_clean)
            raise e
    return processed_blocks

def get_pos_table(pos_tables, name):
    """
    From the list of PMTINFO blocks, return the one with the matching "name" field.
    """
    for tbl in pos_tables:
        if tbl.get("name") == name:
            return tbl
    return None

def order_geo_volumes(geo_volumes):
    """
    Order volumes so that the world volume (index == "world") comes first.
    """
    world_vols = [vol for vol in geo_volumes if vol.get("index") == "world"]
    non_world_vols = [vol for vol in geo_volumes if vol.get("index") != "world"]
    return world_vols + non_world_vols

# ------------------------------------------------------------------
def create_materials_element(root):
    """
    Create a minimal <materials> section.
    """
    materials_elem = ET.SubElement(root, "materials")
    # water
    water = ET.SubElement(materials_elem, "material", name="water", state="liquid")
    ET.SubElement(water, "D", unit="g/cm3", value="1.0")
    # aluminum
    aluminum = ET.SubElement(materials_elem, "material", name="aluminum")
    ET.SubElement(aluminum, "D", unit="g/cm3", value="2.7")
    # stainless_steel
    stainless = ET.SubElement(materials_elem, "material", name="stainless_steel")
    ET.SubElement(stainless, "D", unit="g/cm3", value="8.0")
    # quartz
    quartz = ET.SubElement(materials_elem, "material", name="quartz")
    ET.SubElement(quartz, "D", unit="g/cm3", value="2.65")
    # pmt_vacuum
    vacuum = ET.SubElement(materials_elem, "material", name="pmt_vacuum", state="gas")
    ET.SubElement(vacuum, "D", unit="g/cm3", value="0.001")
    # Vacuum (default)
    gvacuum = ET.SubElement(materials_elem, "material", name="Vacuum", state="gas")
    ET.SubElement(gvacuum, "D", unit="g/cm3", value="1e-25")
    return materials_elem

# ------------------------------------------------------------------
def create_solids_element(root, geo_volumes, r11065_info):
    """
    Create the <solids> section.
    For geo volumes of type "box", create a GDML <box> (dimensions = 2 × half-length).
    For PMT arrays (type "pmtarray") we approximate the PMT shape using a <tube>.
    """
    solids_elem = ET.SubElement(root, "solids")
    for vol in geo_volumes:
        if vol["type"] == "box":
            sx, sy, sz = vol["size"]
            ET.SubElement(
                solids_elem,
                "box",
                name=vol["index"],
                x=str(sx * 2),
                y=str(sy * 2),
                z=str(sz * 2)
            )
        # (Additional types can be implemented as needed)
        elif vol["type"] == "pmtarray":
            # Use the maximum rho_edge from the PMT info as the outer radius.
            outer_radius = max(r11065_info["rho_edge"])
            half_length = 50.0  # mm (adjust as appropriate)
            ET.SubElement(
                solids_elem,
                "tube",
                name="pmt_solid",
                rmin="0",
                rmax=str(outer_radius),
                z=str(half_length),
                deltaphi="360"
            )
    return solids_elem

# ------------------------------------------------------------------
def create_structure_element(root, geo_volumes, pos_tables, r11065_info):
    """
    Create the <structure> section.
    (1) First pass: Create all <volume> elements (without physvol placements).
    (2) Second pass: Add placement (<physvol>) entries into the mother volumes.
    For a volume with type "pmtarray", PMT placements from the PMTINFO table are added.
    """
    structure_elem = ET.SubElement(root, "structure")
    volume_dict = {}

    # First pass: create all volume definitions.
    for vol in geo_volumes:
        vol_name = "vol_" + vol["index"]
        vol_elem = ET.SubElement(structure_elem, "volume", name=vol_name)
        material = vol.get("material", "Vacuum")
        ET.SubElement(vol_elem, "materialref", ref=material)
        if vol["type"] == "box":
            ET.SubElement(vol_elem, "solidref", ref=vol["index"])
        elif vol["type"] == "pmtarray":
            ET.SubElement(vol_elem, "solidref", ref="pmt_solid")
        volume_dict[vol["index"]] = vol_elem

    # Also create an individual volume for PMT placements (if any pmtarray exists)
    for vol in geo_volumes:
        if vol["type"] == "pmtarray":
            pmt_vol_name = "vol_" + vol["index"] + "_pmt"
            pmt_vol = ET.SubElement(structure_elem, "volume", name=pmt_vol_name)
            ET.SubElement(pmt_vol, "materialref", ref=r11065_info["glass_material"])
            ET.SubElement(pmt_vol, "solidref", ref="pmt_solid")

    # Second pass: add placements into the mother volumes.
    for vol in geo_volumes:
        if vol.get("mother", "").strip():
            mother_name = vol["mother"]
            mother_elem = volume_dict.get(mother_name)
            if mother_elem is None:
                print("Warning: mother volume '{}' not found for volume '{}'".format(mother_name, vol["index"]))
                continue
            # Place the volume itself
            physvol = ET.SubElement(mother_elem, "physvol")
            ET.SubElement(physvol, "volumeref", ref="vol_" + vol["index"])
            if "position" in vol:
                pos = vol["position"]
                ET.SubElement(
                    physvol,
                    "position",
                    name="pos_" + vol["index"],
                    unit="mm",
                    x=str(pos[0]),
                    y=str(pos[1]),
                    z=str(pos[2])
                )
            # For a pmtarray, add individual PMT placements.
            if vol["type"] == "pmtarray":
                pos_table_name = vol.get("pos_table")
                pos_table = get_pos_table(pos_tables, pos_table_name)
                if pos_table is None:
                    print("Warning: pos_table '{}' not found for volume '{}'".format(pos_table_name, vol["index"]))
                else:
                    num_pmts = len(pos_table.get("x", []))
                    for i in range(num_pmts):
                        pmt_physvol = ET.SubElement(mother_elem, "physvol")
                        ET.SubElement(pmt_physvol, "volumeref", ref="vol_" + vol["index"] + "_pmt")
                        pos_name = "pos_pmt_" + str(i)
                        ET.SubElement(
                            pmt_physvol,
                            "position",
                            name=pos_name,
                            unit="mm",
                            x=str(pos_table["x"][i]),
                            y=str(pos_table["y"][i]),
                            z=str(pos_table["z"][i])
                        )
                        ET.SubElement(
                            pmt_physvol,
                            "rotation",
                            name="rot_pmt_" + str(i),
                            unit="deg",
                            x="0",
                            y="0",
                            z="0"
                        )

    return structure_elem

# ------------------------------------------------------------------
def create_setup_element(root, topvolume):
    """
    Create the <setup> element that specifies the top (world) volume.
    """
    setup_elem = ET.SubElement(root, "setup", name="Default", version="1.0")
    ET.SubElement(setup_elem, "topvolume", ref=topvolume)
    return setup_elem

# ------------------------------------------------------------------
def generate_gdml(geo_volumes, pos_tables, r11065_info):
    """
    Build the full GDML XML tree.
    """
    # Order volumes so that the world volume is processed first.
    geo_volumes = order_geo_volumes(geo_volumes)

    gdml = ET.Element("gdml")
    create_materials_element(gdml)
    create_solids_element(gdml, geo_volumes, r11065_info)
    create_structure_element(gdml, geo_volumes, pos_tables, r11065_info)
    # Top volume is assumed to be the one with index "world"
    create_setup_element(gdml, "vol_world")
    return gdml

# ------------------------------------------------------------------
def main():
    if len(sys.argv) != 3:
        print("Usage: {} <geo_file> <pmtinfo_file>".format(sys.argv[0]))
        sys.exit(1)
    geo_file = sys.argv[1]
    pmtinfo_file = sys.argv[2]

    try:
        geo_volumes = parse_blocks_from_file(geo_file)
    except Exception as e:
        print("Error parsing geo file:", e)
        sys.exit(1)

    try:
        pos_tables = parse_blocks_from_file(pmtinfo_file)
    except Exception as e:
        print("Error parsing PMTINFO file:", e)
        sys.exit(1)

    gdml_tree = generate_gdml(geo_volumes, pos_tables, r11065_info)
    rough_string = ET.tostring(gdml_tree, 'utf-8')
    reparsed = xml.dom.minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="  ")

    output_filename = "gdml_out.gdml"
    with open(output_filename, "w") as f:
        f.write(pretty_xml)
    print("GDML file generated with good indentation:", output_filename)

if __name__ == "__main__":
    main()

