#!/bin/python3

from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.forcefield import minimizer
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv


# file naming section

# read in structure and read out structure
fname = '/.../bachelor-thesis/starting-structures/KRAS-G12D-pdb5xco_1_A_prepared_minimized.mae'  # name of file to load

# naming all files

# sasa input files
sasa_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-abs1.h5'
empirical_sasa_values = '/.../bachelor-thesis/data/empirical_sasa_values_for_rsa.csv'
sasa_empirical_ref_list = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa_empirical_ref_list.txt'

# rsa files
rsa_res_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-abs1.txt'
rsa_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-abs1.h5'
rsa_res_abs_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-abs-diag1.csv'
rsa_res_abs_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-abs1_heatmap.png'
rsa_res_abs_diag_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-abs-diag1_barplot.png'
rsa_res_ref_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-ref-abs1.txt'
rsa_res_ref_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-ref-abs1.h5'
rsa_res_ref_abs_csv_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-ref-abs1.csv'
rsa_res_ref_abs_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-ref-abs1_barplot.png'
rsa_res_rel_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/RSA/KRAS-G12D-pdb5xco_rsa-res-rel1_barplot.png'

pbd_sequence = "MTEYKLVVVGADGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHK"


csv_file_data = pd.read_csv(empirical_sasa_values, sep=';')  # values from this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3836772/
sasa_amino_code = csv_file_data['AA Code'].tolist()
theo_sasa_values = csv_file_data['Empirical SASA'].tolist()
empirical_sasa_of_sequence_list = []


for AA in pbd_sequence:
    index = 0
    for code in sasa_amino_code:
        if AA == code:
            empirical_sasa_of_sequence_list.append(theo_sasa_values[index])
            break
        index += 1


# writing empirical SASA values for reference structure to txt files
with open(sasa_empirical_ref_list, 'a') as f:
    f.write("empirical_SASA_of_each_AA_in_ref_sequence:\n" + str(empirical_sasa_of_sequence_list) + "\n")


# calculate reference values section
# creating reader and writer
st = structure.StructureReader.read(fname)

# minimize wildtype and make starting reference calculation for SASA, energy and stability
minimizer.minimize_structure(st)  # minimize wildtype

sasa_res_ref = analyze.calculate_sasa_by_residue(st)  # make list of SASA from wildtype as reference
if len(sasa_res_ref) > 167:  # remove data in list after aminoacid, like GDP and Mg2+
    len_dif = len(sasa_res_ref) - 167
    for l in range(0, len_dif):
        sasa_res_ref.pop()


# RSA calculation from SASA arrays

# read in SASA array from file
sasa_res_abs_array_readin = pd.read_hdf(sasa_res_abs_pd_file)
sasa_res_abs_array_readin_buffer = np.array(sasa_res_abs_array_readin)

cys_value = 148  # empirical value of Cys (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3836772/)

# calculate the SASA res array to RSA res array
for res_nr_in_array in range(0, 167):
    rsa_list_buffer = []
    for mut_res_nr_in_array in range(0, 167):
        sasa_value = sasa_res_abs_array_readin_buffer[res_nr_in_array, mut_res_nr_in_array]
        if mut_res_nr_in_array == res_nr_in_array:  # if its the diagonal value in the array it is the mutated AA -> Cys
            rsa_list_buffer.append(sasa_value / cys_value)
        else:
            rsa_list_buffer.append(sasa_value / empirical_sasa_of_sequence_list[mut_res_nr_in_array])
    if res_nr_in_array == 0:
        rsa_res_abs_array = np.array([rsa_list_buffer])
    else:
        rsa_res_abs_array = np.append(rsa_res_abs_array, [rsa_list_buffer], axis=0)  # absolute sasa
    
rsa_res_ref_array = []
   
# calculate the SASA res ref array to RSA res ref array
for mut_res_nr_in_array in range(0, 167):
    sasa_value = sasa_res_ref[mut_res_nr_in_array]
    rsa_res_ref_array.append([sasa_value / empirical_sasa_of_sequence_list[mut_res_nr_in_array]])



# writing to file section

# writing SASA to txt files
with open(rsa_res_abs_file, 'a') as f:
    f.write("RSA_res_mut-absolut:\n" + str(rsa_res_abs_array) + "\n")
with open(rsa_res_ref_abs_file, 'a') as f:
    f.write("RSA_res_ref-absolut:\n" + str(rsa_res_ref_array) + "\n")


# creating pandas array of SASA then writing to file
rsa_res_abs_pd_array = pd.DataFrame(rsa_res_abs_array)
with pd.HDFStore(rsa_res_abs_pd_file, 'w') as f:
    f.append('rsa_res_abs_pd_array', rsa_res_abs_pd_array)  # f.append('key', value)
rsa_res_ref_pd_array = pd.DataFrame(rsa_res_ref_array)
with pd.HDFStore(rsa_res_ref_abs_pd_file, 'w') as f:
    f.append('rsa_res_ref_abs_pd_array', rsa_res_ref_pd_array)  # f.append('key', value)


rsa_res_abs_diag_array = []
for x in range(0, 167):
    rsa_res_abs_diag_array.append([rsa_res_abs_array[x, x]])

# writing SASA to csv files
with open(rsa_res_abs_diag_file, 'a') as f:
    csv.writer(f).writerows(rsa_res_abs_diag_array)
with open(rsa_res_ref_abs_csv_file, 'a') as f:
    csv.writer(f).writerows(rsa_res_ref_array)


# figure section

# can read in the rsa pandas file, gets the arrays, if more arrays in a file, can include the 'key' from writing
rsa_res_abs_array_readin = pd.read_hdf(rsa_res_abs_pd_file)
rsa_res_ref_array_readin_buffer = pd.read_hdf(rsa_res_ref_abs_pd_file)

rsa_res_ref_array_readin = pd.DataFrame(rsa_res_ref_array_readin_buffer).transpose()
rsa_res_abs_diag_array_readin_buffer = np.array(rsa_res_abs_array_readin)
rsa_res_abs_diag_array = []
rsa_res_rel_array = []
for x in range(0, 167):
    rsa_res_abs_diag_array.append([rsa_res_abs_diag_array_readin_buffer[x, x]])
    rsa_res_rel_array.append([rsa_res_ref_array_readin[x] - rsa_res_abs_diag_array_readin_buffer[x, x]])


# plots a heatmap of absolute rsa by residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(rsa_res_abs_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'absolute RSA []'})
ax.set_title("Heatmap of absolute RSA by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA to Cys")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(rsa_res_abs_heatmap_fig)



# plots a barplot of diagonal of absolute rsa by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=rsa_res_abs_diag_array)
ax.set_title("Barplot of diagonal absolute RSA by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("absolute RSA []")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rsa_res_abs_diag_barplot_fig)


# plots a barplot of absolute reference rsa by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=rsa_res_ref_array)
ax.set_title("Barplot of absolute reference RSA by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("absolute RSA []")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rsa_res_ref_abs_barplot_fig)

# plots a barplot of relative rsa by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=rsa_res_rel_array)
ax.set_title("Barplot of relative RSA by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("relative RSA []")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rsa_res_rel_barplot_fig)
