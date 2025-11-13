import sys
import subprocess
import os

def display_gdml_in_root(gdml_file_path):
    """
    Displays a GDML file in ROOT.

    Args:
        gdml_file_path (str): The absolute path to the GDML file.
    """
    if not os.path.exists(gdml_file_path):
        print(f"Error: GDML file not found at {gdml_file_path}")
        return

    # ROOT commands to import and draw the GDML file
    root_commands = f"""
    TGeoManager::Import("{gdml_file_path}");
    if (gGeoManager) {{
        gGeoManager->GetTopVolume()->Draw();
        gApplication->Run(); // Keep the canvas open
    }} else {{
        std::cerr << "Error: TGeoManager not initialized." << std::endl;
    }}
    .q
    """

    print(f"Attempting to display {gdml_file_path} in ROOT...")
    try:
        # Run ROOT in batch mode and execute the commands
        process = subprocess.run(
            ["root", "-l", "-b"],  # -l: no ROOT splash screen, -b: batch mode
            input=root_commands,
            text=True,
            check=True,
            capture_output=True
        )
        print("ROOT Stdout:")
        print(process.stdout)
        print("ROOT Stderr:")
        print(process.stderr)
        print(f"Successfully attempted to display {gdml_file_path} in ROOT.")
    except subprocess.CalledProcessError as e:
        print(f"Error running ROOT: {e}")
        print(f"ROOT Stdout: {e.stdout}")
        print(f"ROOT Stderr: {e.stderr}")
    except FileNotFoundError:
        print("Error: 'root' command not found. Please ensure ROOT is installed and in your PATH.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 display_gdml.py <path_to_gdml_file>")
        sys.exit(1)

    gdml_file = sys.argv[1]
    display_gdml_in_root(gdml_file)
