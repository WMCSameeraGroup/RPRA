import math

def read_xyz_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    return lines


def write_xyz_file(filename, xyz_data):
    with open(filename, 'w') as file:
        file.writelines(xyz_data)


def extract_xyz_files(xyz_movie):
    # Find the indices of the start of each XYZ file
    indices = [i for i, line in enumerate(xyz_movie) if line.strip().isdigit()]

    # Read the first XYZ file
    first_xyz = xyz_movie[indices[0]:indices[1]]

    # Read the last XYZ file
    last_xyz = xyz_movie[indices[-1]:]

    return first_xyz, last_xyz


def calculate_distance(atom1, atom2):
    return math.sqrt((atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2 + (atom1[3] - atom2[3])**2)


def find_closest_h_atoms(xyz_data, cutoff_distance=1.2):
    atoms = []
    for line in xyz_data[2:]:  # Skip the first two lines (header)
        parts = line.split()
        atom_type = parts[0]
        x, y, z = map(float, parts[1:])
        atoms.append((atom_type, x, y, z))

    o_atoms = [atom for atom in atoms if atom[0] == 'O']
    h_atoms = [atom for atom in atoms if atom[0] == 'H']

    closest_h_atoms = []
    for o_atom in o_atoms:
        distances = [(calculate_distance(o_atom, h_atom), h_atom) for h_atom in h_atoms]
        distances = [d for d in distances if d[0] <= cutoff_distance]
        distances.sort(key=lambda x: x[0])
        if len(distances) >= 2:
            closest_h_atoms.extend([distances[0][1], distances[1][1]])
            h_atoms.remove(distances[0][1])
            h_atoms.remove(distances[1][1])

    return len(o_atoms), len(closest_h_atoms)


def find_closest_atoms(xyz_data, cutoff_distance, atom1_type, atom2_type, num_bonds):
    atoms = []
    for line in xyz_data[2:]:  # Skip the first two lines (header)
        parts = line.split()
        atom_type = parts[0]
        x, y, z = map(float, parts[1:])
        atoms.append((atom_type, x, y, z))

    atom1_list = [atom for atom in atoms if atom[0] == atom1_type]
    atom2_list = [atom for atom in atoms if atom[0] == atom2_type]

    closest_atoms = []
    for atom1 in atom1_list:
        distances = [(calculate_distance(atom1, atom2), atom2) for atom2 in atom2_list]
        distances = [d for d in distances if d[0] <= cutoff_distance]
        distances.sort(key=lambda x: x[0])
        if len(distances) >= num_bonds:
            closest_atoms.extend([distances[i][1] for i in range(num_bonds)])
            for i in range(num_bonds):
                atom2_list.remove(distances[i][1])

    return len(atom1_list), len(closest_atoms)


def find_hcn_molecules(xyz_data, h_c_cutoff, c_n_cutoff):
    atoms = []
    for line in xyz_data[2:]:  # Skip the first two lines (header)
        parts = line.split()
        atom_type = parts[0]
        x, y, z = map(float, parts[1:])
        atoms.append((atom_type, x, y, z))

    h_atoms = [atom for atom in atoms if atom[0] == 'H']
    c_atoms = [atom for atom in atoms if atom[0] == 'C']
    n_atoms = [atom for atom in atoms if atom[0] == 'N']

    hcn_molecules = 0
    used_atoms = set()

    for h_atom in h_atoms:
        for c_atom in c_atoms:
            if calculate_distance(h_atom, c_atom) <= h_c_cutoff:
                for n_atom in n_atoms:
                    if calculate_distance(c_atom, n_atom) <= c_n_cutoff:
                        if h_atom not in used_atoms and c_atom not in used_atoms and n_atom not in used_atoms:
                            hcn_molecules += 1
                            used_atoms.update([h_atom, c_atom, n_atom])
                        break
                break
    return hcn_molecules