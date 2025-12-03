#!/usr/bin/env python3
"""
This script converts custom geo and PMTINFO files into a compliant GDML file.
"""
import sys
import re
import json
import xml.etree.ElementTree as ET
import xml.dom.minidom
from math import acos, atan2, degrees

# --- Data Structures ---

ATOMIC_MASSES = {}

R11065_INFO = {
    "name": "PMT_r11065",
    "solid_name": "solid_PMT_r11065",
    "volume_name": "vol_PMT_r11065",
    "material": "quartz",
    "solid": {
        "type": "sphere",
        "radius": 25.4,  # 1-inch radius
    }
}

# --- 1. File Parsing ---
def parse_ratdb_file(filename):
    """Parses a RATDB file into a list of Python dictionaries."""
    try:
        with open(filename, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Error: File not found '{filename}'", file=sys.stderr)
        sys.exit(1)

    # Remove comments and handle trailing commas before parsing
    content = re.sub(r'//.*?\n', '\n', content)
    content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)
    content = re.sub(r',\s*([}\]])', r'\1', content)

    # Use a regex to find all top-level JSON-like objects
    blocks = []
    for match in re.finditer(r'{\s*name:\s*"([^"]+)"\s*,\s*[^}]+}', content):
        block_text = match.group(0)
        # Make it valid JSON
        json_text = re.sub(r'(\w+)(?=\s*:)', r'"\1"', block_text)
        try:
            data = json.loads(json_text)
            blocks.append(data)
        except json.JSONDecodeError as e:
            print(f"Warning: Could not parse block in '{filename}'. Skipping. Error: {e}", file=sys.stderr)
            print("Problematic block text:", block_text)

    return blocks

# --- 2. GDML Component Definition Functions ---

def define_materials(gdml_root):
    """Creates the <materials> section with necessary materials."""
    materials_elem = ET.SubElement(gdml_root, "materials")
    
    # Define some basic elements
    define_elem = ET.SubElement(gdml_root, "define")
    ET.SubElement(define_elem, "element", name="H", formula="H", Z="1", atom="1.008")
    ET.SubElement(define_elem, "element", name="O", formula="O", Z="8", atom="16.00")

    # Water
    water_mat = ET.SubElement(materials_elem, "material", name="water", state="liquid")
    ET.SubElement(water_mat, "D", value="1.0", unit="g/cm3")
    ET.SubElement(water_mat, "composite", n="2", ref="H")
    ET.SubElement(water_mat, "composite", n="1", ref="O")

    # Quartz (for PMTs)
    ET.SubElement(materials_elem, "material", name="quartz", state="solid").set("Z", "14")
    
    # Air (for world volume)
    air_mat = ET.SubElement(materials_elem, "material", name="Air", state="gas")
    ET.SubElement(air_mat, "D", value="0.00129", unit="g/cm3")
    ET.SubElement(air_mat, "fraction", n="0.7", ref="N")
    ET.SubElement(define_elem, "element", name="N", formula="N", Z="7", atom="14.01")
    ET.SubElement(air_mat, "fraction", n="0.3", ref="O")


def define_solids(gdml_root, geo_volumes):
    """Creates the <solids> section."""
    solids_elem = ET.SubElement(gdml_root, "solids")
    for vol in geo_volumes:
        if vol.get("type") == "box":
            ET.SubElement(solids_elem, "box", name=vol["index"],
                          x=str(vol["size"][0]*2), y=str(vol["size"][1]*2), z=str(vol["size"][2]*2),
                          lunit="mm")

    # PMT solid
    pmt_solid_info = R11065_INFO["solid"]
    ET.SubElement(solids_elem, "sphere", name=R11065_INFO["solid_name"],
                  rmax=str(pmt_solid_info["radius"]), lunit="mm")


def get_rotation(dir_x, dir_y, dir_z):
    """Calculate rotation angles from a direction vector."""
    theta = degrees(acos(dir_z))
    phi = degrees(atan2(dir_y, dir_x))
    return {"z": str(phi), "y": str(theta), "x": "0"}


def define_structure(gdml_root, geo_volumes, pmt_pos_table):
    """Creates the <structure> section."""
    structure_elem = ET.SubElement(gdml_root, "structure")

    # Create logical volumes first
    for vol in geo_volumes:
        vol_elem = ET.SubElement(structure_elem, "volume", name=f"vol_{vol['index']}")
        ET.SubElement(vol_elem, "materialref", ref=vol.get("material", "Air"))
        ET.SubElement(vol_elem, "solidref", ref=vol['index'])

    # PMT logical volume
    pmt_vol = ET.SubElement(structure_elem, "volume", name=R11065_INFO["volume_name"])
    ET.SubElement(pmt_vol, "materialref", ref=R11065_INFO["material"])
    ET.SubElement(pmt_vol, "solidref", ref=R11065_INFO["solid_name"])

    # Now, place physical volumes
    for vol in sorted(geo_volumes, key=lambda v: v.get('mother', '')):
        if "mother" in vol:
            mother_vol = structure_elem.find(f".//volume[@name='vol_{vol['mother']}']")
            if mother_vol is not None:
                physvol = ET.SubElement(mother_vol, "physvol", name=f"phys_{vol['index']}")
                ET.SubElement(physvol, "volumeref", ref=f"vol_{vol['index']}")
                if "position" in vol:
                    ET.SubElement(physvol, "position", name=f"pos_{vol['index']}", unit="mm",
                                  x=str(vol["position"][0]), y=str(vol["position"][1]), z=str(vol["position"][2]))

    # Place PMTs
    mother_of_pmts_name = next((v["mother"] for v in geo_volumes if v.get("type") == "pmtarray"), None)
    if mother_of_pmts_name:
        mother_of_pmts_vol = structure_elem.find(f".//volume[@name='vol_{mother_of_pmts_name}']")
        if mother_of_pmts_vol is not None and pmt_pos_table:
            for i in range(len(pmt_pos_table["x"])):
                pmt_physvol = ET.SubElement(mother_of_pmts_vol, "physvol", name=f"phys_pmt_{i}")
                ET.SubElement(pmt_physvol, "volumeref", ref=R11065_INFO["volume_name"])
                ET.SubElement(pmt_physvol, "position", name=f"pos_pmt_{i}", unit="mm",
                              x=str(pmt_pos_table["x"][i]), y=str(pmt_pos_table["y"][i]), z=str(pmt_pos_table["z"][i]))
                rot = get_rotation(pmt_pos_table["dir_x"][i], pmt_pos_table["dir_y"][i], pmt_pos_table["dir_z"][i])
                ET.SubElement(pmt_physvol, "rotation", name=f"rot_pmt_{i}", unit="deg", x=rot["x"], y=rot["y"], z=rot["z"])


def define_setup(gdml_root, world_volume_name):
    """Creates the <setup> section."""
    setup_elem = ET.SubElement(gdml_root, "setup", name="Default", version="1.0")
    ET.SubElement(setup_elem, "world", ref=f"vol_{world_volume_name}")


def prettify(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = xml.dom.minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

# --- 3. Main Orchestrator ---
def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <geo_file> <pmtinfo_file> <pmt_table_name>", file=sys.stderr)
        sys.exit(1)

    geo_file, pmtinfo_file, pmt_table_name = sys.argv[1], sys.argv[2], sys.argv[3]
    
    geo_data = parse_ratdb_file(geo_file)
    pmtinfo_data = parse_ratdb_file(pmtinfo_file)
    
    geo_volumes = [b for b in geo_data if "type" in b]
    
    pmt_pos_table = next((b for b in pmtinfo_data if b.get("name") == pmt_table_name), None)
    if not pmt_pos_table:
        print(f"Error: PMT table '{pmt_table_name}' not found in '{pmtinfo_file}'", file=sys.stderr)
        sys.exit(1)

    world_volume = next((v for v in geo_volumes if v.get("index") == "world"), None)
    if not world_volume:
        print("Error: No 'world' volume with index: \"world\" found in geometry file.", file=sys.stderr)
        sys.exit(1)

    gdml_root = ET.Element("gdml")
    define_materials(gdml_root)
    define_solids(gdml_root, geo_volumes)
    define_structure(gdml_root, geo_volumes, pmt_pos_table)
    define_setup(gdml_root, world_volume["index"])

    output_filename = "out.gdml"
    with open(output_filename, "w") as f:
        f.write(prettify(gdml_root))
    print(f"GDML file generated successfully: {output_filename}")

if __name__ == "__main__":
    main()
