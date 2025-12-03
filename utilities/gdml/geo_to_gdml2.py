import json
import numpy as np
import pyg4ometry
import sys
import re

# Attempt to import matplotlib for plotting
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# Attempt to import scipy for robust rotation calculations
try:
    from scipy.spatial.transform import Rotation as R
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def plot_pmts(pmtinfo):
    """
    Creates a 3D plot of PMT positions and orientations.
    """
    if not MATPLOTLIB_AVAILABLE:
        print("\nMatplotlib is not installed. Skipping plot.")
        print("To install it, run: pip install matplotlib")
        return

    print("\nGenerating 3D plot of PMT geometry...")
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Extract data for plotting
    x = pmtinfo["x"]
    y = pmtinfo["y"]
    z = pmtinfo["z"]
    u = pmtinfo["dir_x"]
    v = pmtinfo["dir_y"]
    w = pmtinfo["dir_z"]

    # Create a quiver plot to show PMTs as arrows
    # The length is scaled for better visualization.
    ax.quiver(x, y, z, u, v, w, length=100, normalize=True, color='b')

    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Z (mm)')
    ax.set_title('PMT Positions and Orientations')
    
    # Set aspect ratio to be equal to avoid distortion
    # Handle case where all points are the same
    if max(x) == min(x) and max(y) == min(y) and max(z) == min(z):
        ax.set_xlim(max(x)-1, max(x)+1)
        ax.set_ylim(max(y)-1, max(y)+1)
        ax.set_zlim(max(z)-1, max(z)+1)
    else:
        max_range = np.array([max(x)-min(x), max(y)-min(y), max(z)-min(z)]).max() / 2.0
        mid_x = (max(x)+min(x)) * 0.5
        mid_y = (max(y)+min(y)) * 0.5
        mid_z = (max(z)+min(z)) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)


    plt.show()

def create_gdml(ratdb_file, config_name, output_file):
    """
    Parses a specific ratdb file format, finds the specified configuration,
    and creates a GDML file from it.

    Args:
        ratdb_file (str): The path to the ratdb file.
        config_name (str): The name of the configuration to use.
        output_file (str): The path to the output GDML file.
    """
    if not SCIPY_AVAILABLE:
        print("\nScipy is not installed. This script requires it for robust rotation calculations.")
        print("To install it, run: pip install scipy")
        return

    # 1. New parsing strategy: Bypass JSON decoding entirely.
    try:
        with open(ratdb_file, 'r') as f:
            text = f.read()

        # Step 1: Find the start of the desired configuration block.
        start_pattern = r'{\s*name:\s*"' + re.escape(config_name) + r'"'
        start_match = re.search(start_pattern, text)
        
        if not start_match:
            print(f"Error: Could not find the start of configuration block for '{config_name}' in {ratdb_file}.")
            return
            
        # Step 2: Extract the text from the start of our block to the end of the file.
        text_from_start = text[start_match.start():]
        
        # Step 3: Find the matching closing brace for this block.
        open_braces = 0
        end_index = -1
        for i, char in enumerate(text_from_start):
            if char == '{':
                open_braces += 1
            elif char == '}':
                open_braces -= 1
                if open_braces == 0:
                    end_index = i + 1
                    break
        
        if end_index == -1:
            print(f"Error: Could not find the closing brace for configuration block '{config_name}'.")
            return

        config_block = text_from_start[:end_index]
        
        # Step 4: Extract each array individually using regex from the isolated block.
        pmtinfo = {"name": config_name}
        data_keys = ["x", "y", "z", "dir_x", "dir_y", "dir_z"]
        
        for key in data_keys:
            # This pattern finds the key, and captures everything between the following '[' and ']'
            array_pattern = r'' + key + r'\s*:\s*\[([\s\S]*?)\]'
            array_match = re.search(array_pattern, config_block)
            if array_match:
                # Get the content, split by comma, and convert to float
                # This handles trailing commas by filtering out empty strings after split.
                content_str = array_match.group(1)
                values = [float(v) for v in content_str.split(',') if v.strip()]
                pmtinfo[key] = values
            else:
                raise ValueError(f"Could not find data array for key '{key}' in the config block.")

    except FileNotFoundError:
        print(f"Error: The file {ratdb_file} was not found.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during file parsing: {e}")
        return

    # 2. Create a pyg4ometry registry
    reg = pyg4ometry.geant4.Registry()

    # 3. Define materials
    world_material = pyg4ometry.geant4.MaterialPredefined("G4_WATER")
    el_B = pyg4ometry.geant4.ElementSimple("Boron", "B", 5, 10.81)
    el_O = pyg4ometry.geant4.ElementSimple("Oxygen", "O", 8, 16.00)
    el_Na = pyg4ometry.geant4.ElementSimple("Sodium", "Na", 11, 22.99)
    el_Al = pyg4ometry.geant4.ElementSimple("Aluminum", "Al", 13, 26.98)
    el_Si = pyg4ometry.geant4.ElementSimple("Silicon", "Si", 14, 28.09)
    pmt_material = pyg4ometry.geant4.MaterialCompound("Pyrex", 2.23, 5, reg)
    pmt_material.add_element_massfraction(el_Si, 0.3768)
    pmt_material.add_element_massfraction(el_B, 0.0403)
    pmt_material.add_element_massfraction(el_Na, 0.0297)
    pmt_material.add_element_massfraction(el_Al, 0.0122)
    pmt_material.add_element_massfraction(el_O, 0.5410)
    reg.addMaterial(world_material)

    # 4. Define the world volume
    world_solid = pyg4ometry.geant4.solid.Box("world_solid", 4000, 4000, 4000, reg, lunit="mm")
    world_logical = pyg4ometry.geant4.LogicalVolume(world_solid, world_material, "world_logical", reg)
    reg.setWorld(world_logical.name)

    # 5. Define the PMT solid
    pmt_solid = pyg4ometry.geant4.solid.Cons("pmt_solid", 0, 25.4, 0, 0, 25.4, 0, 2 * np.pi, reg, "mm", "rad")
    pmt_logical = pyg4ometry.geant4.LogicalVolume(pmt_solid, pmt_material, "pmt_logical", reg)

    # 7. Place the PMTs in the world volume
    print("\n--- PMT Placement Information ---")
    for i in range(len(pmtinfo["x"])):
        position = [pmtinfo["x"][i], pmtinfo["y"][i], pmtinfo["z"][i]]
        direction = np.array([pmtinfo["dir_x"][i], pmtinfo["dir_y"][i], pmtinfo["dir_z"][i]])
        
        print(f"PMT {i}: Pos = {position}, Dir = {direction.tolist()}")

        norm = np.linalg.norm(direction)
        if norm == 0:
            direction = np.array([0., 0., 1.])
        else:
            direction = direction / norm

        # The cone solid's base faces its local -Z axis.
        # We need to rotate this local -Z vector to align with the 'direction' from the file.
        base_axis = np.array([0., 0., -1.])
        
        # Use Scipy to find the rotation that aligns the base_axis with the target direction.
        rotation_object, _ = R.align_vectors([direction], [base_axis])
        
        # Get the rotation as ZYX extrinsic Euler angles in radians. This is a common
        # convention that corresponds to a sequence of rotations about the fixed axes.
        # The GDML 'rotation' tag with x, y, z attributes corresponds to this sequence.
        euler_angles = rotation_object.as_euler('ZYX', degrees=False)
        rot_z, rot_y, rot_x = euler_angles[0], euler_angles[1], euler_angles[2]

        # Create the GDML rotation object
        rotation = pyg4ometry.gdml.Rotation("rotation_{}".format(i), rot_x, rot_y, rot_z, 'rad', reg)
        pmt_physical = pyg4ometry.geant4.PhysicalVolume(rotation, position, pmt_logical, "pmt_{}".format(i), world_logical, reg)

    # 8. Write the GDML file
    writer = pyg4ometry.gdml.Writer()
    writer.addDetector(reg)
    writer.write(output_file)
    print(f"\nSuccessfully created GDML file: {output_file}")
    
    # After creating the GDML, plot the geometry
    plot_pmts(pmtinfo)


if __name__ == '__main__':
    # Define the file and the specific configuration name you want to use.
    ratdb_file_to_use = "PMTINFO.txt"
    config_to_build = "DOMINFO_pen_small_short" 
    
    create_gdml(ratdb_file_to_use, config_to_build, "detector.gdml")

