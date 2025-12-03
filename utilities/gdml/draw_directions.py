import matplotlib.pyplot as plt
import numpy as np

def visualize_pmts_with_indices():
    """
    This function creates a 3D plot of PMT positions and their direction vectors.
    The positions and directions are hardcoded from the user's data.
    Each PMT is annotated with its index number.
    The plot is interactive, allowing for rotation and zooming.
    """

    # --- PMT Positions (in mm) ---
    # Copied from the GDML file provided earlier.
    positions = np.array([
        [0.0, -195.95494889, 0.0],
        [0.014519849, -165.571658, 108.535225],
        [-93.98700187, -165.5716577, 54.28018691],
        [-94.00152172, -165.5716577, -54.25503779],
        [-0.014519849, -165.571658, -108.535225],
        [93.98700187, -165.5716577, -54.28018691],
        [94.00152172, -165.5716577, 54.25503779],
        [0.0187642369, -54.5766634, 190.299738],
        [-166.975329, -106.348057, -0.0100552998],
        [-164.7950254, -54.57666342, 95.16611935],
        [-83.4789564, -106.3480569, -144.60990444],
        [-164.81378963, -54.57666342, -95.13361874],
        [83.49637269, -106.3480569, -144.59984914],
        [-0.0187642369, -54.5766634, -190.299738],
        [166.975329, -106.348057, 0.0100552998],
        [164.7950254, -54.57666342, -95.16611935],
        [83.4789564, -106.3480569, 144.60990444],
        [164.81378963, -54.57666342, 95.13361874],
        [-83.49637269, -106.3480569, 144.59984914],
        [0.00758869613, 106.365682, 166.975496],
        [-190.309996, 54.5766634, 0.0357977582],
        [-144.60122658, 106.36568243, 83.49431976],
        [-95.18599956, 54.57666342, -164.79539189],
        [-144.60881528, 106.36568243, -83.48117575],
        [95.12399602, 54.57666342, -164.83118965],
        [-0.00758869613, 106.365682, -166.975496],
        [190.309996, 54.5766634, -0.0357977582],
        [144.60122658, 106.36568243, -83.49431976],
        [95.18599956, 54.57666342, 164.79539189],
        [144.60881528, 106.36568243, 83.48117575],
        [-95.12399602, 54.57666342, 164.83118965]
    ])

    # --- PMT Direction Vectors ---
    # These are the normalized vectors indicating the facing of each PMT.
    directions = np.array([
        [0., -1., 0.],
        [7.33420476e-05, -8.36328561e-01, 5.48228541e-01],
        [-0.47474317, -0.83632856, 0.27417779],
        [-0.47481651, -0.83632856, -0.27405075],
        [-7.33420476e-05, -8.36328561e-01, -5.48228541e-01],
        [0.47474317, -0.83632856, -0.27417779],
        [0.47481651, -0.83632856, 0.27405075],
        [9.47826516e-05, -2.75679789e-01, 9.61249523e-01],
        [-8.43453151e-01, -5.37202736e-01, -5.07929782e-05],
        [-0.83241911, -0.27567979, 0.48070685],
        [-0.42168259, -0.53720274, -0.73047725],
        [-0.8325139, -0.27567979, -0.48054268],
        [0.42177056, -0.53720274, -0.73042646],
        [-9.47826516e-05, -2.75679789e-01, -9.61249523e-01],
        [8.43453151e-01, -5.37202736e-01, 5.07929782e-05],
        [0.83241911, -0.27567979, -0.48070685],
        [0.42168259, -0.53720274, 0.73047725],
        [0.8325139, -0.27567979, 0.48054268],
        [-0.42177056, -0.53720274, 0.73042646],
        [3.83314051e-05, 5.37265690e-01, 8.43413052e-01],
        [-9.61253449e-01, 2.75666056e-01, 1.80814037e-04],
        [-0.73039796, 0.53726569, 0.42173972],
        [-0.48078331, 0.27566606, -0.8323795],
        [-0.73043629, 0.53726569, -0.42167333],
        [0.48047013, 0.27566606, -0.83256031],
        [-3.83314051e-05, 5.37265690e-01, -8.43413052e-01],
        [9.61253449e-01, 2.75666056e-01, -1.80814037e-04],
        [0.73039796, 0.53726569, -0.42173972],
        [0.48078331, 0.27566606, 0.8323795],
        [0.73043629, 0.53726569, 0.42167333],
        [-0.48047013, 0.27566606, 0.83256031]
    ])

    # --- Create the 3D Plot ---
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Unpack the position and direction coordinates
    x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
    u, v, w = directions[:, 0], directions[:, 1], directions[:, 2]

    # Plot the PMT positions as blue dots
    ax.scatter(x, y, z, color='blue', s=50, label='PMT Positions')

    # Plot the direction vectors as red arrows (quivers)
    ax.quiver(x, y, z, u, v, w, length=50, normalize=True, color='red', label='Direction')
    
    # --- ADDED: Annotate each PMT with its index ---
    for i, pos in enumerate(positions):
        ax.text(pos[0], pos[1], pos[2], f'  {i}', color='black', fontsize=9)

    # --- Formatting the Plot ---
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Z (mm)')
    ax.set_title('3D Visualization of PMT Positions, Directions, and Indices')
    ax.legend()

    # Set axis limits to be equal to avoid distortion
    max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
    mid_x = (x.max()+x.min()) * 0.5
    mid_y = (y.max()+y.min()) * 0.5
    mid_z = (z.max()+z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Show the plot
    plt.show()


if __name__ == '__main__':
    # To run this script, you need to have matplotlib and numpy installed.
    # You can install them using pip:
    # pip install matplotlib numpy
    visualize_pmts_with_indices()
