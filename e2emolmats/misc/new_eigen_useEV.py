import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

# Number of atoms
N = 92

# Initialize a 3N x 3N matrix with zeros

datafile = open("FORCE_CONSTANTS.txt", "r")
lines = datafile.readlines()
print(lines[0])
lines.pop(0)
print(lines[0])


def read_dynamical_matrix(input_lines, N):
    dyn_mat = np.zeros((3 * N, 3 * N))

    for i in range(0, len(input_lines), 4):
        # Read atom indices and convert them to 0-based index
        indices = input_lines[i].split()
        atom_i = int(indices[0]) - 1
        atom_j = int(indices[1]) - 1
        # Read the 3x3 block for this pair of atoms
        block = np.array([
            list(map(float, input_lines[i + 1].split())),
            list(map(float, input_lines[i + 2].split())),
            list(map(float, input_lines[i + 3].split()))
        ])

        # Fill the corresponding block in the dynamical matrix
        dyn_mat[3 * atom_i:3 * atom_i + 3, 3 * atom_j:3 * atom_j + 3] = block
        # Since the dynamical matrix is symmetric, fill the transpose block
        if atom_i != atom_j:
            dyn_mat[3 * atom_j:3 * atom_j + 3, 3 * atom_i:3 * atom_i + 3] = block.T

    return dyn_mat


# Read the matrix data from a file (replace 'dynamical_matrix.txt' with your file path)
dynamical_matrix = read_dynamical_matrix(lines, N)
# Calculate eigenvalues and eigenvectors
eigenvalues_real, eigenvectors_real = np.linalg.eigh(dynamical_matrix)
# Print results
#print("Eigenvalues_real_unit:", eigenvalues_real)
#print("Eigenvectors_real_unit:", eigenvectors_real)
#print(len(eigenvalues_real)) This was 3 * N
kcal_to_J = 4184  # kcal/mol to J/mol
g_mol_to_kg = 1e-3  # g/mol to kg/mol
angstrom_to_m = 1e-10  # Å to m
speed_of_light = 29979245800  #cm/s
PlanckConstant = 6.62607015e-34  # [J s]
THzToJ = PlanckConstant * 1e12  #[J]
Kb = 1.3806504e-23
temp = 300
Na = 6.02214179e23

JTokJmol = 1e-3 * Na

# Convert eigenvalues to (rad/s)²
# Conversion factor: sqrt((J/(kg·m²)))
conversion_factor = np.sqrt(kcal_to_J / (g_mol_to_kg)) / (angstrom_to_m) / (2 * np.pi) / 1e12

#conversion_factor = np.sqrt(kcal_to_J / (g_mol_to_kg * (angstrom_to_m ** 2)))
frequencies = np.sqrt(np.abs(eigenvalues_real.real)) * np.sign(eigenvalues_real.real) * conversion_factor  #[THz]
#This is phonopy's way for negative eigenvalue
#angular_frequencies = np.sqrt(np.abs(eigenvalues_si)) * np.sign(eigenvalues_si)
# Calculate vibrational frequencies (Hz)
#frequencies = angular_frequencies / (2 * np.pi)
#k = inverse of wavelength, freq / speed of light
wavenumber = frequencies / speed_of_light * 1e12  #[cm^-1]
p_frequencies = []
VEnergy = []
for i in range(0, len(frequencies)):
    if frequencies[i] > 0.0:
        p_frequencies.append(frequencies[i] * THzToJ)  #[J]
#VEnergy = frequencies * THzToJ / 2 + Kb * temp * np.log(1.0-np.exp((-frequencies * THzToJ) / (Kb * temp)))
for i in range(0, len(p_frequencies)):
    VEnergy.append((p_frequencies[i] / 2 + Kb * temp * np.log(1.0 - np.exp((-p_frequencies[i]) / (Kb * temp)))) * JTokJmol)  #[kJ/mol]
print(sum(VEnergy))

sns.kdeplot(wavenumber, bw_adjust=0.2)
plt.hist(wavenumber, bins=50, density=True, alpha=0.3)

plt.xlabel('Wavenumber (cm$^{-1}$)')
plt.ylabel('Density of States')
plt.title('Phonon Density of States of form 2 acridine')
plt.grid(True)
plt.show()

# Define the number of bins or the range of wavenumbers
#num_bins = 100
# Compute the histogram
#hist, bin_edges = np.histogram(wavenumber, bins=num_bins, density=True)
#hist, bin_edges = np.histogram(frequencies, bins=num_bins, density=True)

# Plot the DOS
#plt.figure(figsize=(8, 6))
#plt.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), edgecolor="black", align="edge")
#plt.xlabel('Wavenumber (cm$^{-1}$)')
#plt.ylabel('Density of States')
#plt.title('Phonon Density of States of form 2 acridine')
#plt.grid(True)
#plt.show()

#print(frequencies)
