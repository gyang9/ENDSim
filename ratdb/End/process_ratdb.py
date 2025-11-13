import re
import math

def parse_ratdb_entry(entry_text):
    name_match = re.search(r'name:\s*"(.*?)"', entry_text)
    if not name_match:
        return None
    name = name_match.group(1)

    def get_float_array(field_name):
        match = re.search(fr'{field_name}:\s*\[(.*?)\]', entry_text, re.DOTALL)
        if not match:
            return []
        return [float(v) for v in match.group(1).replace('\n', '').split(',') if v.strip()]

    x = get_float_array('x')
    y = get_float_array('y')
    z = get_float_array('z')
    dir_x = get_float_array('dir_x')
    dir_y = get_float_array('dir_y')
    dir_z = get_float_array('dir_z')

    if not (len(x) == len(y) == len(z) == len(dir_x) == len(dir_y) == len(dir_z)):
        print(f"Warning: Mismatched array lengths in {name}, skipping.")
        return None

    pmts = []
    for i in range(len(x)):
        pmts.append({
            'position': {'x': x[i], 'y': y[i], 'z': z[i]},
            'direction': {'x': dir_x[i], 'y': dir_y[i], 'z': dir_z[i]}
        })
    return {'name': name, 'pmts': pmts}

def create_gdml(config_data):
    name = config_data['name']
    pmts = config_data['pmts']

    gdml = '<?xml version="1.0" encoding="UTF-8"?>\n'
    gdml += '<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">\n'

    gdml += '  <define>\n'
    for i, pmt in enumerate(pmts):
        pos = pmt['position']
        gdml += f'    <position name="pos_{i}" x="{pos["x"]}" y="{pos["y"]}" z="{pos["z"]}" unit="cm"/>\n'
        
        direction = pmt['direction']
        dx, dy, dz = direction['x'], direction['y'], direction['z']
        
        mag = math.sqrt(dx**2 + dy**2 + dz**2)
        if mag > 0:
            dx /= mag
            dy /= mag
            dz /= mag

        if (dx*dx + dy*dy) > 1e-9:
            alpha_z = math.degrees(math.atan2(dy, dx))
            beta_y = math.degrees(math.atan2(math.sqrt(dx*dx+dy*dy), dz))
            gamma_x = 0
        else:
            alpha_z = 0
            if dz > 0:
                beta_y = 0
            else:
                beta_y = 180
            gamma_x = 0
            
        gdml += f'    <rotation name="rot_{i}" x="{gamma_x}" y="{beta_y}" z="{alpha_z}" unit="deg"/>\n'

    gdml += '  </define>\n'

    gdml += '  <materials>\n'
    gdml += '    <material name="Vacuum" Z="1.0" density="1.e-25" unit="g/cm3">\n'
    gdml += '      <atom ref="H"/>\n'
    gdml += '    </material>\n'
    gdml += '  </materials>\n'

    gdml += '  <solids>\n'
    gdml += '    <box name="WorldBox" x="20000" y="20000" z="40000" unit="cm"/>\n'
    gdml += '    <sphere name="PMT" rmax="25.4" startphi="0" deltaphi="360" starttheta="0" deltatheta="180" unit="cm"/>\n'
    gdml += '  </solids>\n'

    gdml += '  <structure>\n'
    gdml += '    <volume name="PMTVolume">\n'
    gdml += '      <material_ref ref="Vacuum"/>\n'
    gdml += '      <solid_ref ref="PMT"/>\n'
    gdml += '    </volume>\n'
    gdml += '    <volume name="World">\n'
    gdml += '      <material_ref ref="Vacuum"/>\n'
    gdml += '      <solid_ref ref="WorldBox"/>\n'
    for i in range(len(pmts)):
        gdml += f'      <physvol>\n'
        gdml += f'        <volumeref ref="PMTVolume"/>\n'
        gdml += f'        <positionref ref="pos_{i}"/>
'        gdml += f'        <rotationref ref="rot_{i}"/>
'
        gdml += f'      </physvol>\n'
    gdml += '    </volume>\n'
    gdml += '  </structure>\n'

    gdml += '  <setup name="Default" version="1.0">\n'
    gdml += '    <world ref="World"/>\n'
    gdml += '  </setup>\n'

    gdml += '</gdml>\n'
    return gdml

with open("/home/guang/work/END/ENDSim/ratdb/End/PMTINFO.ratdb", "r") as f:
    content = f.read()

entries = content.split('}')
for entry in entries:
    if not entry.strip():
        continue
    if not entry.strip().startswith('{'):
        entry = '{' + entry
    config_data = parse_ratdb_entry(entry + '}')
    if config_data:
        gdml_content = create_gdml(config_data)
        file_name = f"{config_data['name']}.gdml"
        with open(file_name, "w") as f:
            f.write(gdml_content)
        print(f"Generated {file_name}")
