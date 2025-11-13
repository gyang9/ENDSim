import re

def parse_pmt_data(file_content):
    """
    Parses the PMTINFO.txt file content to extract configurations.

    Args:
        file_content: The string content of the PMTINFO.txt file.

    Returns:
        A list of dictionaries, where each dictionary represents a configuration.
    """
    configurations = []
    # Use regex to find all configuration blocks
    config_blocks = re.findall(r'{\s*name:\s*"(.*?)"(.*?)}(?=\s*{|$)', file_content, re.DOTALL)
    for name, data in config_blocks:
        config = {'name': name}
        # Extract numerical data for x, y, z coordinates
        for coord in ['x', 'y', 'z']:
            match = re.search(r'' + coord + r':\s*\[(.*?)\]', data, re.DOTALL)
            if match:
                # Convert the comma-separated string of numbers to a list of floats
                config[coord] = [float(f) for f in match.group(1).split(',') if f.strip()]
        configurations.append(config)
    return configurations

def write_gdml(config):
    """
    Writes a single configuration to a GDML file.

    Args:
        config: A dictionary representing a single PMT configuration.
    """
    # Sanitize the filename
    filename = "".join(c for c in config['name'] if c.isalnum() or c in (' ', '_')).rstrip() + "2nd.gdml"
    
    with open(filename, 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">\n')
        
        # Define materials
        f.write('<materials>\n')
        f.write('    <material name="Air" Z="7.6">\n')
        f.write('        <D value="1.29" unit="mg/cm3"/>\n')
        f.write('    </material>\n')
        f.write('</materials>\n')

        # Define solids
        f.write('<solids>\n')
        f.write('    <box name="WorldBox" x="10000" y="10000" z="40000" lunit="mm"/>\n')
        f.write('    <sphere name="PMT" rmin="0" rmax="100" lunit="mm"/>\n')
        f.write('</solids>\n')

        # Define structure
        f.write('<structure>\n')
        f.write('    <volume name="lvPMT">\n')
        f.write('        <materialref ref="Air"/>\n')
        f.write('        <solidref ref="PMT"/>\n')
        f.write('    </volume>\n')
        f.write('    <volume name="World">\n')
        f.write('        <materialref ref="Air"/>\n')
        f.write('        <solidref ref="WorldBox"/>\n')

        # Place PMTs
        for i in range(len(config['x'])):
            f.write(f'        <physvol name="pvPMT_{i}">\n')
            f.write(f'            <volumeref ref="lvPMT"/>\n')
            f.write(f'            <position name="posPMT_{i}" unit="mm" x="{config["x"][i]}" y="{config["y"][i]}" z="{config["z"][i]}"/>\n')
            f.write('        </physvol>\n')
        
        f.write('    </volume>\n')
        f.write('</structure>\n')

        # Setup
        f.write('<setup name="Default" version="1.0">\n')
        f.write('    <world ref="World"/>\n')
        f.write('</setup>\n')
        
        f.write('</gdml>\n')
    print(f"Generated {filename}")

# Read the content of the uploaded file
with open('PMTINFO.txt', 'r') as f:
    file_content = f.read()

# Parse the data and generate GDML files
configurations = parse_pmt_data(file_content)
for config in configurations:
    write_gdml(config)
