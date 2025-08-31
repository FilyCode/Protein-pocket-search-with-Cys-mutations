#!/bin/python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv


# file naming section

# naming all files
# read in structure and read out structure
fname = '/.../bachelor-thesis/starting-structures/KRAS-G12D-pdb5xco_1_A_prepared_minimized.mae'  # name of file to load
modified_fname = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_mutation_output.mae'  # name of exported file

# sasa files
sasa_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-abs1.h5'
sasa_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-rel1.h5'
sasa_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-abs1.h5'
sasa_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-rel1.h5'
sasa_res_abs_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-abs-diag1.csv'
sasa_res_rel_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-rel-diag1.csv'
sasa_res_abs_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-abs1_heatmap.png'
sasa_res_rel_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-rel1_heatmap.png'
sasa_abs_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-abs1_barplot.png'
sasa_rel_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-rel1_barplot.png'
sasa_res_abs_diag_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-abs-diag1_barplot.png'
sasa_res_rel_diag_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_sasa-res-rel-diag1_barplot.png'

# energy files
# potential energy files
energy_pot_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-abs1.h5'
energy_pot_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-rel1.h5'
energy_pot_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-abs1.h5'
energy_pot_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-rel1.h5'
energy_pot_res_abs_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-abs-diag1.csv'
energy_pot_res_rel_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-rel-diag1.csv'
# prime energy files
energy_prime_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-abs1.h5'
energy_prime_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-rel1.h5'
# interaction energy files
energy_int_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-abs1.h5'
energy_int_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-rel1.h5'
energy_int_res_abs_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-abs-diag1.csv'
energy_int_res_rel_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-rel-diag1.csv'
# energy figures
energy_pot_abs_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-abs1_barplot.png'
energy_pot_rel_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-rel1_barplot.png'
energy_pot_res_abs_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-abs1_heatmap.png'
energy_pot_res_rel_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-rel1_heatmap.png'
energy_pot_res_abs_diag_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-abs-diag1_barplot.png'
energy_pot_res_rel_diag_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-pot-res-rel-diag1_barplot.png'
energy_prime_abs_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-abs1_barplot.png'
energy_prime_rel_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-prime-rel1_barplot.png'
energy_int_res_abs_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-abs1_heatmap.png'
energy_int_res_rel_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-rel1_heatmap.png'
energy_int_res_abs_diag_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-abs-diag1_barplot.png'
energy_int_res_rel_diag_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco_energy-int-res-rel-diag1_barplot.png'



# figure section

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
    sasa_res_abs_diag_array.append([sasa_res_abs_diag_array_readin_buffer[x, x]])
    sasa_res_rel_diag_array.append([sasa_res_rel_diag_array_readin_buffer[x, x]])



# plots a heatmap of absolute SASA by residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(sasa_res_abs_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'absolute SASA [$\AA^2$]'})
ax.set_title("Heatmap of absolute SASA by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA to Cys")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(sasa_res_abs_heatmap_fig)

# plots a heatmap of relative SASA by residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(sasa_res_rel_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'relative SASA [$\AA^2$]'})
ax.set_title("Heatmap of relative SASA by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA to Cys")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(sasa_res_rel_heatmap_fig)

# plots a barplot of absolute SASA
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=sasa_abs_array_readin)
ax.set_title("Barplot of absolute total SASA")
ax.set_xlabel('mutated AA')
ax.set_ylabel("absolute total SASA [$\AA^2$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(sasa_abs_barplot_fig)

# plots a barplot of relative SASA
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=sasa_rel_array_readin)
ax.set_title("Barplot of relative total SASA")
ax.set_xlabel('mutated AA')
ax.set_ylabel("relative total SASA [$\AA^2$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(sasa_rel_barplot_fig)

# plots a barplot of diagonal of absolute SASA by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=sasa_res_abs_diag_array)
ax.set_title("Barplot of diagonal absolute SASA by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("absolute SASA [$\AA^2$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(sasa_res_abs_diag_barplot_fig)

# plots a barplot of diagonal of relative SASA by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=sasa_res_rel_diag_array)
ax.set_title("Barplot of diagonal relative SASA by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("relative SASA [$\AA^2$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(sasa_res_rel_diag_barplot_fig)


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
    energy_pot_res_abs_diag_array.append([energy_pot_res_abs_diag_array_readin_buffer[x, x]])
    energy_pot_res_rel_diag_array.append([energy_pot_res_rel_diag_array_readin_buffer[x, x]])
    energy_int_res_abs_diag_array.append([energy_int_res_abs_diag_array_readin_buffer[x, x]])
    energy_int_res_rel_diag_array.append([energy_int_res_rel_diag_array_readin_buffer[x, x]])



# writing SASA and energy diag to txt files
with open(sasa_res_abs_diag_file, 'a') as f:
    csv.writer(f).writerows(sasa_res_abs_diag_array)
with open(sasa_res_rel_diag_file, 'a') as f:
    csv.writer(f).writerows(sasa_res_rel_diag_array)
with open(energy_pot_res_abs_diag_file, 'a') as f:
    csv.writer(f).writerows(energy_pot_res_abs_diag_array)
with open(energy_pot_res_rel_diag_file, 'a') as f:
    csv.writer(f).writerows(energy_pot_res_rel_diag_array)
with open(energy_int_res_abs_diag_file, 'a') as f:
    csv.writer(f).writerows(energy_int_res_abs_diag_array)
with open(energy_int_res_rel_diag_file, 'a') as f:
    csv.writer(f).writerows(energy_int_res_rel_diag_array)



# plots a heatmap of absolute potential Energy by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(energy_pot_res_abs_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'Potential Energy [kcal/mol]'})
ax.set_title("Heatmap of absolute Potential Energy by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(energy_pot_res_abs_heatmap_fig)

# plots a heatmap of relative potential Energy by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(energy_pot_res_rel_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'Potential Energy [kcal/mol]'})
ax.set_title("Heatmap of relative Potential Energy by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(energy_pot_res_rel_heatmap_fig)

# plots a barplot of absolute potential energy
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=energy_pot_abs_array_readin)
ax.set_title("Barplot of absolute total Potential Energy")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Potential Energy [kcal/mol]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(energy_pot_abs_barplot_fig)

# plots a barplot of relative potential energy
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=energy_pot_rel_array_readin)
ax.set_title("Barplot of relative total Potential Energy")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Potential Energy [kcal/mol]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(energy_pot_rel_barplot_fig)



# plots a barplot of absolute prime energy
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=energy_prime_abs_array_readin)
ax.set_title("Barplot of absolute total Prime Energy")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Prime Energy [kcal/mol]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(energy_prime_abs_barplot_fig)

# plots a barplot of relative prime energy
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=energy_prime_rel_array_readin)
ax.set_title("Barplot of relative total Prime Energy")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Prime Energy [kcal/mol]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(energy_prime_rel_barplot_fig)



# plots a heatmap of absolute interaction Energy by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(energy_int_res_abs_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'Interaction Energy [kcal/mol]'})
ax.set_title("Heatmap of absolute Interaction Energy by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(energy_int_res_abs_heatmap_fig)


# plots a heatmap of relative interaction Energy by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(energy_int_res_rel_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'Interaction Energy [kcal/mol]'})
ax.set_title("Heatmap of relative Interaction Energy by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(energy_int_res_rel_heatmap_fig)

# plots a barplot of diagonal absolute potential energy by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=energy_pot_res_abs_diag_array)
ax.set_title("Barplot of diagonal absolute Potential Energy by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Potential Energy [kcal/mol]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(energy_pot_res_abs_diag_fig)

# plots a barplot of diagonal relative potential energy by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=energy_pot_res_rel_diag_array)
ax.set_title("Barplot of diagonal relative Potential Energy by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Potential Energy [kcal/mol]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(energy_pot_res_rel_diag_fig)

# plots a barplot of diagonal absolute interaction energy by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=energy_int_res_abs_diag_array)
ax.set_title("Barplot of diagonal absolute Interaction Energy by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Interaction Energy [kcal/mol]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(energy_int_res_abs_diag_fig)

# plots a barplot of diagonal relative interaction energy by residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=energy_int_res_rel_diag_array)
ax.set_title("Barplot of diagonal relative Interaction Energy by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Interaction Energy [kcal/mol]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(energy_int_res_rel_diag_fig)
