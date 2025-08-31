#!/bin/python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv

# file naming section

# naming all files
# pKa files
pKa_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-abs1.h5'
pKa_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-rel1.h5'
pKa_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-abs1.h5'
pKa_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-rel1.h5'
pKa_ref_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-ref1.h5'
pKa_res_ref_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-ref1.h5'
pKa_res_ref_csv_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-ref1.csv'
pKa_res_abs_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-abs-diag1.csv'
pKa_res_rel_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-rel-diag1.csv'
# pKa figures
pKa_res_abs_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-abs1_heatmap.png'
pKa_res_rel_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-rel1_heatmap.png'
pKa_abs_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-abs1_barplot.png'
pKa_rel_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-rel1_barplot.png'
pKa_res_abs_diag_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-abs-diag1_barplot.png'
pKa_res_rel_diag_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-rel-diag1_barplot.png'
pKa_ref_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-ref1_barplot.png'
pKa_res_ref_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-ref1_barplot.png'


# figure section

# can read in the pKa pandas file, gets the arrays, if more arrays in a file, can include the 'key' from writing
pKa_abs_array_readin = pd.read_hdf(pKa_abs_pd_file)
pKa_rel_array_readin = pd.read_hdf(pKa_rel_pd_file)
pKa_res_abs_array_readin = pd.read_hdf(pKa_res_abs_pd_file)
pKa_res_rel_array_readin = pd.read_hdf(pKa_res_rel_pd_file)
pKa_ref_readin = pd.read_hdf(pKa_ref_pd_file)
pKa_res_ref_readin = pd.read_hdf(pKa_res_ref_pd_file)


pKa_res_ref_readin_buffer = np.array(pKa_res_ref_readin)
pKa_res_abs_diag_array_readin_buffer = np.array(pKa_res_abs_array_readin)
pKa_res_rel_diag_array_readin_buffer = np.array(pKa_res_rel_array_readin)
pKa_res_ref_array = []
pKa_res_abs_diag_array = []
pKa_res_rel_diag_array = []
for x in range(0, 167):
    pKa_res_ref_array.append([pKa_res_ref_readin_buffer[0, x]])
    pKa_res_abs_diag_array.append([pKa_res_abs_diag_array_readin_buffer[x, x]])
    pKa_res_rel_diag_array.append([pKa_res_rel_diag_array_readin_buffer[x, x]])


# writing pKa ref and diag to file
# writing SASA and energy diag to txt files
with open(pKa_res_ref_csv_file, 'a') as f:
    csv.writer(f).writerows(pKa_res_ref_array)
with open(pKa_res_abs_diag_file, 'a') as f:
    csv.writer(f).writerows(pKa_res_abs_diag_array)
with open(pKa_res_rel_diag_file, 'a') as f:
    csv.writer(f).writerows(pKa_res_rel_diag_array)



# plots a heatmap of absolute pKa by Residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(pKa_res_abs_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'absolute p$K_a$'})
ax.set_title("Heatmap of absolute p$K_a$ by Residue")
ax.set_xlabel("Amino acids")
ax.set_ylabel("mutated AA to Cys")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=360)
plt.savefig(pKa_res_abs_heatmap_fig)

# plots a heatmap of relative pKa by Residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(pKa_res_rel_array_readin, center=0.0, cmap="PiYG", cbar_kws={'label': 'relative p$K_a$'})
ax.set_title("Heatmap of relative p$K_a$ by Residue")
ax.set_xlabel("Amino acids")
ax.set_ylabel("mutated AA to Cys")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=360)
plt.savefig(pKa_res_rel_heatmap_fig)

# plots a barplot of absolute pKa
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=pKa_abs_array_readin)
ax.set_title("Barplot of absolute total p$K_a$")
ax.set_xlabel('mutated AA')
ax.set_ylabel("p$K_a$")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(pKa_abs_barplot_fig)

# plots a barplot of relative pKa
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=pKa_rel_array_readin)
ax.set_title("Barplot of relative total p$K_a$")
ax.set_xlabel('mutated AA')
ax.set_ylabel("p$K_a$")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(pKa_rel_barplot_fig)


# plots a barplot of diagonal of absolute pKa by Residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=pKa_res_abs_diag_array)
ax.set_title("Barplot of diagonal absolute p$K_a$ by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("p$K_a$")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(pKa_res_abs_diag_barplot_fig)

# plots a barplot of diagonal of relative pKa by Residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=pKa_res_rel_diag_array)
ax.set_title("Barplot of diagonal relative p$K_a$ by Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("p$K_a$")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(pKa_res_rel_diag_barplot_fig)

# plots a barplot of reference total pKa
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=pKa_ref_readin)
ax.set_title("Barplot of total reference p$K_a$")
ax.set_xlabel('Amino acids')
ax.set_ylabel("p$K_a$")
plt.savefig(pKa_ref_barplot_fig)

# plots a barplot of reference pKa by Residue
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=pKa_res_ref_readin)
ax.set_title("Barplot of reference p$K_a$ by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("p$K_a$")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(pKa_res_ref_barplot_fig)