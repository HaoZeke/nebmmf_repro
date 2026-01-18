import os
import glob
import re
import numpy as np
from ase.io import read, write


def organize_with_ase(delta=0.01):
    config_content = """[Main]
job = saddle_search
temperature = 300
random_seed = 706253457

# [Potential]
# potential = metatomic
# [Metatomic]
# model_path = pet-mad-v1.1.0.pt

[Potential]
potential = xtb
[XTBPot]
paramset = GFN2xTB

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 1000
max_move = 0.05
convergence_metric = norm

[LBFGS]
lbfgs_memory = 25
lbfgs_inverse_curvature = 0.01
lbfgs_auto_scale = true

[Saddle Search]
displace_least_coordinated_weight = 1.0
displace_radius = 3.3
displace_magnitude = 0.01
min_mode_method = dimer
max_energy = 10.0

[Debug]
write_movies = True

[Dimer]
improved = true
opt_method = cg
remove_rotation = False
converged_angle = 0.5
"""

    # Identify peak files
    mode_files = glob.glob("peak*_mode.dat")

    for mode_file in mode_files:
        match = re.search(r"peak(\d+)_mode\.dat", mode_file)
        if not match:
            continue

        idx = match.group(1)
        dir_name = f"run_{idx}"
        orig_pos = f"peak{idx}_pos.con"

        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        # Define destination paths
        dest_pos = os.path.join(dir_name, "pos.con")
        dest_dir = os.path.join(dir_name, "direction.dat")
        dest_disp = os.path.join(dir_name, "displacement.con")
        dest_config = os.path.join(dir_name, "config.ini")

        # Handle position and direction files
        if os.path.exists(orig_pos):
            atoms = read(orig_pos)
            write(dest_pos, atoms)
            os.remove(orig_pos)  # Clean up original

        if os.path.exists(mode_file):
            # Read the eigenvector direction from the mode file
            # We assume a simple Nx3 array format in peakXX_mode.dat
            direction_vector = np.loadtxt(mode_file)
            np.savetxt(dest_dir, direction_vector, fmt="%12.8f")
            os.remove(mode_file)  # Clean up original

        # Generate displacement.con using finite difference
        if "atoms" in locals() and "direction_vector" in locals():
            # Create a copy to displace
            displaced_atoms = atoms.copy()

            # Apply finite difference: R_new = R_old + delta * direction
            # ASE Atoms.positions is a numpy array (N, 3)
            displaced_atoms.positions += delta * direction_vector

            # Save the new configuration
            write(dest_disp, displaced_atoms)

        # Write the config.ini
        with open(dest_config, "w") as f:
            f.write(config_content)

        print(f"Success: Organized {dir_name} using ASE.")


if __name__ == "__main__":
    organize_with_ase()
