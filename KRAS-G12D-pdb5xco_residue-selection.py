#!/bin/python3

import numpy as np
import pandas as pd


# file naming section

# naming all files
# sasa files
sasa_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-abs1.h5'
sasa_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-rel1.h5'
sasa_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-abs1.h5'
sasa_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-rel1.h5'

# energy files
# potential energy files
energy_pot_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-abs1.h5'
energy_pot_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-rel1.h5'
energy_pot_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-abs1.h5'
energy_pot_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-rel1.h5'
# prime energy files
energy_prime_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-abs1.h5'
energy_prime_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-rel1.h5'
# interaction energy files
energy_int_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-abs1.h5'
energy_int_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-rel1.h5'

# stability files
stab_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_residue-scanning_stability.csv'

# residue selection files
res_sel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_res-sel1.txt'



# can read in the SASA pandas file, gets the arrays, if more arrays in a file, can include the 'key' from writing
sasa_abs_array_readin = pd.read_hdf(sasa_abs_pd_file)
sasa_rel_array_readin = pd.read_hdf(sasa_rel_pd_file)
sasa_res_abs_array_readin = pd.read_hdf(sasa_res_abs_pd_file)
sasa_res_rel_array_readin = pd.read_hdf(sasa_res_rel_pd_file)

sasa_res_abs_diag_array_readin_buffer = np.array(sasa_res_abs_array_readin)
sasa_res_rel_diag_array_readin_buffer = np.array(sasa_res_rel_array_readin)
sasa_res_abs_diag_array = []
sasa_res_rel_diag_array = []
for x in range(0, 167):
    sasa_res_abs_diag_array.append(sasa_res_abs_diag_array_readin_buffer[x, x])
    sasa_res_rel_diag_array.append(sasa_res_rel_diag_array_readin_buffer[x, x])


# can read in the energy pandas file, gets the arrays, if more arrays in a file, can include the 'key' from writing
energy_pot_abs_array_readin = pd.read_hdf(energy_pot_abs_pd_file)
energy_pot_rel_array_readin = pd.read_hdf(energy_pot_rel_pd_file)
energy_pot_res_abs_array_readin = pd.read_hdf(energy_pot_res_abs_pd_file)
energy_pot_res_rel_array_readin = pd.read_hdf(energy_pot_res_rel_pd_file)
energy_prime_abs_array_readin = pd.read_hdf(energy_prime_abs_pd_file)
energy_prime_rel_array_readin = pd.read_hdf(energy_prime_rel_pd_file)
energy_int_res_abs_array_readin = pd.read_hdf(energy_int_res_abs_pd_file)
energy_int_res_rel_array_readin = pd.read_hdf(energy_int_res_rel_pd_file)

energy_pot_res_abs_diag_array_readin_buffer = np.array(energy_pot_res_abs_array_readin)
energy_pot_res_rel_diag_array_readin_buffer = np.array(energy_pot_res_rel_array_readin)
energy_int_res_abs_diag_array_readin_buffer = np.array(energy_int_res_abs_array_readin)
energy_int_res_rel_diag_array_readin_buffer = np.array(energy_int_res_rel_array_readin)
energy_pot_res_abs_diag_array = []
energy_pot_res_rel_diag_array = []
energy_int_res_abs_diag_array = []
energy_int_res_rel_diag_array = []
for x in range(0, 167):
    energy_pot_res_abs_diag_array.append(energy_pot_res_abs_diag_array_readin_buffer[x, x])
    energy_pot_res_rel_diag_array.append(energy_pot_res_rel_diag_array_readin_buffer[x, x])
    energy_int_res_abs_diag_array.append(energy_int_res_abs_diag_array_readin_buffer[x, x])
    energy_int_res_rel_diag_array.append(energy_int_res_rel_diag_array_readin_buffer[x, x])


# read in stability

stab_file_data = pd.read_csv(stab_file, index_col=0)
stab_list = stab_file_data['d Stability (solvated)'].tolist()

# selection section
selected_res = []
sasa_out = []
energy_pot_out = []
energy_int_out = []
stab_out = []

for r in range(0, 167):
    selected_res_count = int(1)

    if sasa_res_abs_diag_array[r] < 29.6:
        sasa_out.append(r+1)
        selected_res_count = int(0)
    else:
        if energy_pot_res_rel_diag_array[r] > 100:
            energy_pot_out.append(r+1)
            selected_res_count = int(0)
        else:
            if energy_int_res_rel_diag_array[r] > 100:
                energy_int_out.append(r+1)
                selected_res_count = int(0)
            else:
                if stab_list[r] > 20:
                    stab_out.append(r+1)
                    selected_res_count = int(0)

    if selected_res_count == 1:
        selected_res.append(r+1)

with open(res_sel_file, 'a') as f:
    f.write("SASA kicked out:\n" + "Nr: " + str(len(sasa_out)) + "\nkicked out residues: " + str(sasa_out) + "\n\n")
    f.write("Energy pot kicked out:\n" + "Nr: " + str(len(energy_pot_out)) + "\nkicked out residues: " + str(energy_pot_out) + "\n\n")
    f.write("Energy int kicked out:\n" + "Nr: " + str(len(energy_int_out)) + "\nkicked out residues: " + str(energy_int_out) + "\n\n")
    f.write("Stability kicked out:\n" + "Nr: " + str(len(stab_out)) + "\nkicked out residues: " + str(stab_out) + "\n\n")
    f.write("Residues remaining:\n" + "Nr: " + str(len(selected_res)) + "\nremaining residues: " + str(selected_res) + "\n\n")


