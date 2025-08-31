#!/bin/python3

from schrodinger import structure
from schrodinger.structutils import analyze, build
from schrodinger.forcefield import minimizer
from schrodinger.application.bioluminate import protein
import numpy as np
import pandas as pd


# file naming section

# naming all files
# read in structure and read out structure
fname = '/.../bachelor-thesis/starting-structures/KRAS-G12D-pdb5xco_1_A_prepared_minimized.mae'  # name of file to load
modified_fname = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_mutation_output.mae'  # name of exported file

# sasa files
sasa_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-abs1.txt'
sasa_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-rel1.txt'
sasa_res_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-abs1.txt'
sasa_res_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-rel1.txt'
sasa_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-abs1.h5'
sasa_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-rel1.h5'
sasa_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-abs1.h5'
sasa_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-rel1.h5'

# energy files
# potential energy files
energy_pot_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-abs1.txt'
energy_pot_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-rel1.txt'
energy_pot_res_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-abs1.txt'
energy_pot_res_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-rel1.txt'
energy_pot_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-abs1.h5'
energy_pot_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-rel1.h5'
energy_pot_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-abs1.h5'
energy_pot_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-rel1.h5'
# prime energy files
energy_prime_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-abs1.txt'
energy_prime_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-rel1.txt'
energy_prime_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-abs1.h5'
energy_prime_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-rel1.h5'
# interaction energy files
energy_int_res_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-abs1.txt'
energy_int_res_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-rel1.txt'
energy_int_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-abs1.h5'
energy_int_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-rel1.h5'



# calculate reference values section

# creating reader and writer
st = structure.StructureReader.read(fname)
writer = structure.StructureWriter(modified_fname)
st.title = str("KRAS-G12D-pdb5xco_reference-not-minimized")
writer.append(st)  # write starting structure to file

# minimize wildtype and make starting reference calculation for SASA, energy and stability
minimizer.minimize_structure(st)  # minimize wildtype
st.title = str("KRAS-G12D-pdb5xco_reference")
writer.append(st)  # write minimized reference structure to file
sasa_ref = analyze.calculate_sasa(st)
sasa_res_ref = analyze.calculate_sasa_by_residue(st)  # make list of SASA from wildtype as reference
energy_pot_ref = protein.PropertyCalculator(st, "energy_pot_calc").getTotalPotentialEnergy()
energy_prime_ref = protein.PropertyCalculator(st, "energy_prime_calc").getTotalPrimeEnergy()
energy_pot_res_ref = []
energy_int_res_ref = []
for res in st.residue:
    energy_pot_res_ref.append(protein.PropertyCalculator(st, "energy_pot_res_calc").getResiduePotentialEnergy(res))
    energy_int_res_ref.append(protein.PropertyCalculator(st, "energy_int_res_calc").getResidueInteractionEnergy(res))

if len(sasa_res_ref) > 167:  # remove data in list after aminoacid, like GDP and Mg2+
    len_dif = len(sasa_res_ref) - 167
    for l in range(0, len_dif):
        sasa_res_ref.pop()
        energy_pot_res_ref.pop()
        energy_int_res_ref.pop()


i = int(0)  # for loop counting

# mutation calculation section

# looping through the protein by residue and then mutate, minimize, calculate SASA and
# Energy (potential, prime, interaction)
for res in st.residue:  # iterates over all residues of the structure
    i += 1  # for numbering residues and stopping while testing

    if i == 168:  # i == 168, to stop after all aminoacids
        break

    # structure reset, mutation, minimization
    work_st = structure.StructureReader.read(fname)  # reset structure
    if i != 17:  # AS17 (serin) has interaction with magnesium, so do not mutate it
        atom_in_res = res.getAtomIndices()[0]
        build.mutate(work_st, atom_in_res, 'CYS')  # mutate structure
        minimizer.minimize_structure(work_st)  # minimize structure
        work_st.title = str("KRAS-G12D-pdb5xco-" + str(i) + "C")
    else:
        minimizer.minimize_structure(work_st)  # minimize structure
        work_st.title = str("KRAS-G12D-pdb5xco-" + str(i) + "-not-mutated")
    writer.append(work_st)

    if i == 1:
        # sasa calculations and creating arrays
        sasa = analyze.calculate_sasa(work_st)
        sasa_abs_array = np.array([sasa])
        sasa_rel_array = np.array([sasa - sasa_ref])
        sasa_res = analyze.calculate_sasa_by_residue(work_st)  # sasa by residue
        if len(sasa_res) > 167:
            len_dif = len(sasa_res) - 167
            for l in range(0, len_dif):
                sasa_res.pop()
        sasa_res_abs_array = np.array([sasa_res])  # absolute sasa values by residues

        sasa_buffer = []
        for x in range(0, 167):  # loops through list, calculate relative changes in sasa by residue
            sasa_buffer.append(sasa_res[x] - sasa_res_ref[x])
        sasa_res_rel_array = np.array([sasa_buffer])

        # energy calculations and creating arrays
        energy_pot = protein.PropertyCalculator(work_st, "energy_calc").getTotalPotentialEnergy()
        energy_pot_abs_array = np.array([energy_pot])
        energy_pot_rel_array = np.array([energy_pot - energy_pot_ref])
        energy_prime = protein.PropertyCalculator(work_st, "energy_calc").getTotalPrimeEnergy()
        energy_prime_abs_array = np.array([energy_prime])
        energy_prime_rel_array = np.array([energy_prime - energy_prime_ref])
        energy_pot_res = []
        energy_int_res = []
        for res1 in work_st.residue:  # calculate energy by residue
            energy_pot_res.append(protein.PropertyCalculator(work_st, "energy_pot_res_calc").getResiduePotentialEnergy(res1))
            energy_int_res.append(protein.PropertyCalculator(work_st, "energy_int_res_calc").getResidueInteractionEnergy(res1))

        if len(energy_pot_res) > 167:
            len_dif = len(energy_pot_res) - 167
            for l in range(0, len_dif):
                energy_pot_res.pop()
                energy_int_res.pop()
        energy_pot_res_abs_array = np.array([energy_pot_res])  # absolute energy values by residue
        energy_int_res_abs_array = np.array([energy_int_res])

        energy_pot_buffer = []
        energy_int_buffer = []
        for x in range(0, 167):  # loops through list, calculate relative changes in energy by residue
            energy_pot_buffer.append(energy_pot_res[x] - energy_pot_res_ref[x])
            energy_int_buffer.append(energy_int_res[x] - energy_int_res_ref[x])
        energy_pot_res_rel_array = np.array([energy_pot_buffer])
        energy_int_res_rel_array = np.array([energy_int_buffer])

    else:
        # sasa calculations and writing to arrays
        sasa = analyze.calculate_sasa(work_st)
        sasa_abs_array = np.append(sasa_abs_array, [sasa], axis=0)  # absolute sasa
        sasa_rel_array = np.append(sasa_rel_array, [sasa - sasa_ref], axis=0)  # relative sasa
        sasa_res = analyze.calculate_sasa_by_residue(work_st)  # sasa by residue
        if len(sasa_res) > 167:
            len_dif = len(sasa_res) - 167
            for l in range(0, len_dif):
                sasa_res.pop()
        sasa_res_abs_array = np.append(sasa_res_abs_array, [sasa_res], axis=0)

        sasa_buffer = []
        for x in range(0, 167):  # loops through list, calculate relative changes in sasa by residue
            sasa_buffer.append(sasa_res[x] - sasa_res_ref[x])
        sasa_res_rel_array = np.append(sasa_res_rel_array, [sasa_buffer], axis=0)

        # energy calculations and writing to arrays
        energy_pot = protein.PropertyCalculator(work_st, "energy_calc").getTotalPotentialEnergy()
        energy_pot_abs_array = np.append(energy_pot_abs_array, [energy_pot], axis=0)
        energy_pot_rel_array = np.append(energy_pot_rel_array, [energy_pot - energy_pot_ref], axis=0)
        energy_prime = protein.PropertyCalculator(work_st, "energy_calc").getTotalPrimeEnergy()
        energy_prime_abs_array = np.append(energy_prime_abs_array, [energy_prime], axis=0)
        energy_prime_rel_array = np.append(energy_prime_rel_array, [energy_prime - energy_prime_ref], axis=0)
        energy_pot_res = []
        energy_int_res = []
        for res1 in work_st.residue:  # calculate energy by residue
            energy_pot_res.append(protein.PropertyCalculator(work_st, "energy_pot_res_calc").getResiduePotentialEnergy(res1))
            energy_int_res.append(protein.PropertyCalculator(work_st, "energy_int_res_calc").getResidueInteractionEnergy(res1))

        if len(energy_pot_res) > 167:
            len_dif = len(energy_pot_res) - 167
            for l in range(0, len_dif):
                energy_pot_res.pop()
                energy_int_res.pop()
        energy_pot_res_abs_array = np.append(energy_pot_res_abs_array, [energy_pot_res], axis=0)
        energy_int_res_abs_array = np.append(energy_int_res_abs_array, [energy_int_res], axis=0)

        energy_pot_buffer = []
        energy_int_buffer = []
        for x in range(0, 167):  # loops through list, calculate relative changes in energy by residue
            energy_pot_buffer.append(energy_pot_res[x] - energy_pot_res_ref[x])
            energy_int_buffer.append(energy_int_res[x] - energy_int_res_ref[x])
        energy_pot_res_rel_array = np.append(energy_pot_res_rel_array, [energy_pot_buffer], axis=0)  # write energy to array
        energy_int_res_rel_array = np.append(energy_int_res_rel_array, [energy_int_buffer], axis=0)


# writing to file section

# writing SASA to txt files
with open(sasa_abs_file, 'a') as f:
    f.write("SASA_mut-absolut:\n" + str(sasa_abs_array) + "\n")
with open(sasa_rel_file, 'a') as f:
    f.write("SASA_mut-relative:\n" + str(sasa_rel_array) + "\n")
with open(sasa_res_abs_file, 'a') as f:
    f.write("SASA_mut-residue-absolute:\n" + str(sasa_res_abs_array) + "\n")
with open(sasa_res_rel_file, 'a') as f:
    f.write("SASA_mut-residue-relative:\n" + str(sasa_res_rel_array) + "\n")

# creating pandas array of SASA then writing to file
sasa_abs_pd_array = pd.DataFrame(sasa_abs_array).transpose()
sasa_rel_pd_array = pd.DataFrame(sasa_rel_array).transpose()
sasa_res_abs_pd_array = pd.DataFrame(sasa_res_abs_array)
sasa_res_rel_pd_array = pd.DataFrame(sasa_res_rel_array)
with pd.HDFStore(sasa_abs_pd_file, 'w') as f:
    f.append('sasa_abs_pd_array', sasa_abs_pd_array)  # f.append('key', value)
with pd.HDFStore(sasa_rel_pd_file, 'w') as f:
    f.append('sasa_rel_pd_array', sasa_rel_pd_array)
with pd.HDFStore(sasa_res_abs_pd_file, 'w') as f:
    f.append('sasa_res_abs_pd_array', sasa_res_abs_pd_array)
with pd.HDFStore(sasa_res_rel_pd_file, 'w') as f:
    f.append('sasa_res_rel_pd_array', sasa_res_rel_pd_array)

# writing energy to txt files
with open(energy_pot_abs_file, 'a') as f:
    f.write("Energy_mut_pot_absolut:\n" + str(energy_pot_abs_array) + "\n")
with open(energy_pot_rel_file, 'a') as f:
    f.write("Energy_mut_pot_relative:\n" + str(energy_pot_rel_array) + "\n")
with open(energy_pot_res_abs_file, 'a') as f:
    f.write("Energy_mut_pot_res_absolute:\n" + str(energy_pot_res_abs_array) + "\n")
with open(energy_pot_res_rel_file, 'a') as f:
    f.write("Energy_mut_pot_res_relative:\n" + str(energy_pot_res_rel_array) + "\n")
with open(energy_prime_abs_file, 'a') as f:
    f.write("Energy_mut_prime_absolut mit Ref als 1:\n" + str(energy_prime_abs_array) + "\n")
with open(energy_prime_rel_file, 'a') as f:
    f.write("Energy_mut_prime_relative:\n" + str(energy_prime_rel_array) + "\n")
with open(energy_int_res_abs_file, 'a') as f:
    f.write("Energy_mut_int_res_absolute:\n" + str(energy_int_res_abs_array) + "\n")
with open(energy_int_res_rel_file, 'a') as f:
    f.write("Energy_mut_int_res_relative:\n" + str(energy_int_res_rel_array) + "\n")

# creating pandas array of Energy then writing to file
# potential energy
energy_pot_abs_pd_array = pd.DataFrame(energy_pot_abs_array).transpose()
energy_pot_rel_pd_array = pd.DataFrame(energy_pot_rel_array).transpose()
energy_pot_res_abs_pd_array = pd.DataFrame(energy_pot_res_abs_array)
energy_pot_res_rel_pd_array = pd.DataFrame(energy_pot_res_rel_array)
with pd.HDFStore(energy_pot_abs_pd_file, 'w') as f:
    f.append('energy_pot_abs_pd_array', energy_pot_abs_pd_array)  # f.append('key', value)
with pd.HDFStore(energy_pot_rel_pd_file, 'w') as f:
    f.append('energy_pot_rel_pd_array', energy_pot_rel_pd_array)
with pd.HDFStore(energy_pot_res_abs_pd_file, 'w') as f:
    f.append('energy_pot_res_abs_pd_array', energy_pot_res_abs_pd_array)
with pd.HDFStore(energy_pot_res_rel_pd_file, 'w') as f:
    f.append('energy_pot_res_rel_pd_array', energy_pot_res_rel_pd_array)
# prime energy
energy_prime_abs_pd_array = pd.DataFrame(energy_prime_abs_array).transpose()
energy_prime_rel_pd_array = pd.DataFrame(energy_prime_rel_array).transpose()
with pd.HDFStore(energy_prime_abs_pd_file, 'w') as f:
    f.append('energy_prime_abs_pd_array', energy_prime_abs_pd_array)  # f.append('key', value)
with pd.HDFStore(energy_prime_rel_pd_file, 'w') as f:
    f.append('energy_prime_rel_pd_array', energy_prime_rel_pd_array)
# interaction energy
energy_int_res_abs_pd_array = pd.DataFrame(energy_int_res_abs_array)
energy_int_res_rel_pd_array = pd.DataFrame(energy_int_res_rel_array)
with pd.HDFStore(energy_int_res_abs_pd_file, 'w') as f:
    f.append('energy_int_res_abs_pd_array', energy_int_res_abs_pd_array)
with pd.HDFStore(energy_int_res_rel_pd_file, 'w') as f:
    f.append('energy_int_res_rel_pd_array', energy_int_res_rel_pd_array)
