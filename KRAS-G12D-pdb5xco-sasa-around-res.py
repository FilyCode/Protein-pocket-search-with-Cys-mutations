#!/bin/python3

from schrodinger import structure
from schrodinger.structutils import analyze, build, measure
from schrodinger.forcefield import minimizer
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
sasa_around_res_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-abs1.txt'
sasa_around_res_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-rel1.txt'
sasa_around_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-abs1.h5'
sasa_around_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-rel1.h5'
sasa_around_res_abs_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-abs-diag1.csv'
sasa_around_res_rel_diag_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-rel-diag1.csv'
sasa_around_res_ref_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-ref-abs1.txt'
sasa_around_res_ref_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-ref-abs1.h5'
sasa_around_res_ref_abs_csv_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-ref-abs1.csv'
# sasa figures
sasa_around_res_abs_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa_around-res-abs1_heatmap.png'
sasa_around_res_rel_heatmap_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa_around-res-rel1_heatmap.png'
sasa_around_res_diag_abs_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-diag-abs1_barplot.png'
sasa_around_res_diag_rel_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-diag-rel1_barplot.png'
sasa_around_res_ref_abs_barplot_fig = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/sasa-around-res/KRAS-G12D-pdb5xco_sasa-around-res-ref-abs1_barplot.png'



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
sasa_around_res_ref = []
for x in st.residue:
    res_atom_list = x.getAtomIndices()
    atom_list = measure.get_atoms_close_to_subset(st, res_atom_list, 8.0)
    sasa_ref_calc = analyze.calculate_sasa_by_atom(st, atoms=atom_list)
    sasa_ref = sum(sasa_ref_calc)
    sasa_around_res_ref.append(sasa_ref)


if len(sasa_around_res_ref) > 167:  # remove data in list after aminoacid, like GDP and Mg2+
    len_dif = len(sasa_around_res_ref) - 167
    for l in range(0, len_dif):
        sasa_around_res_ref.pop()

i = int(0)  # for loop counting

# mutation calculation section

# looping through the protein by residue and then mutate, minimize, calculate SASA and
# Energy (potential, prime, interaction)
for res in st.residue:  # iterates over all residues of the structure
    i += 1  # for numbering residues and stopping while testing

    if i == 168:  # i == 168, to stop after all Amino acids
        break

    # structure reset, mutation, minimization
    work_st = structure.StructureReader.read(fname)  # reset structure
    if i != 17:  # AS17 (serin) has interaction with magnesium, so don´t mutate it
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

        sasa_abs_buffer = []
        sasa_rel_buffer = []
        t = 0
        for x in work_st.residue:
            res_atom_list = x.getAtomIndices()
            atom_list = measure.get_atoms_close_to_subset(work_st, res_atom_list, 8.0)
            sasa_calc = analyze.calculate_sasa_by_atom(work_st, atoms=atom_list)
            sasa = sum(sasa_calc)
            sasa_abs_buffer.append(sasa)
            sasa_rel_buffer.append(sasa - sasa_around_res_ref[t])
            t += 1
            if t == 167:
                break

        if len(sasa_abs_buffer) > 167:
            len_dif = len(sasa_abs_buffer) - 167
            for l in range(0, len_dif):
                sasa_abs_buffer.pop()
                sasa_rel_buffer.pop()
        sasa_around_res_abs_array = np.array([sasa_abs_buffer])
        sasa_around_res_rel_array = np.array([sasa_rel_buffer])

    else:
        # sasa calculations and writing to arrays
        sasa_abs_buffer = []
        sasa_rel_buffer = []
        t = 0
        for x in work_st.residue:
            res_atom_list = x.getAtomIndices()
            atom_list = measure.get_atoms_close_to_subset(work_st, res_atom_list, 8.0)
            sasa_calc = analyze.calculate_sasa_by_atom(work_st, atoms=atom_list)
            sasa = sum(sasa_calc)
            sasa_abs_buffer.append(sasa)
            sasa_rel_buffer.append(sasa - sasa_around_res_ref[t])
            t += 1
            if t == 167:
                break

        if len(sasa_abs_buffer) > 167:
            len_dif = len(sasa_abs_buffer) - 167
            for l in range(0, len_dif):
                sasa_abs_buffer.pop()
                sasa_rel_buffer.pop()
        sasa_around_res_abs_array = np.append(sasa_around_res_abs_array, [sasa_abs_buffer], axis=0)
        sasa_around_res_rel_array = np.append(sasa_around_res_rel_array, [sasa_rel_buffer], axis=0)



# writing to file section

# writing SASA to txt files
with open(sasa_around_res_abs_file, 'a') as f:
    f.write("SASA_around_res_mut-absolut:\n" + str(sasa_around_res_abs_array) + "\n")
with open(sasa_around_res_rel_file, 'a') as f:
    f.write("SASA_around_res_mut-relative:\n" + str(sasa_around_res_rel_array) + "\n")
with open(sasa_around_res_ref_abs_file, 'a') as f:
    f.write("SASA_around_res_mut-absolut:\n" + str(sasa_around_res_ref) + "\n")


# creating pandas array of SASA then writing to file
sasa_around_res_abs_pd_array = pd.DataFrame(sasa_around_res_abs_array)
sasa_around_res_rel_pd_array = pd.DataFrame(sasa_around_res_rel_array)
sasa_around_res_ref_abs_pd_array = pd.DataFrame(sasa_around_res_ref)
with pd.HDFStore(sasa_around_res_abs_pd_file, 'w') as f:
    f.append('sasa_around_res_abs_pd_array', sasa_around_res_abs_pd_array)  # f.append('key', value)
with pd.HDFStore(sasa_around_res_rel_pd_file, 'w') as f:
    f.append('sasa_around_res_rel_pd_array', sasa_around_res_rel_pd_array)
with pd.HDFStore(sasa_around_res_ref_abs_pd_file, 'w') as f:
    f.append('sasa_around_res_abs_pd_array', sasa_around_res_ref_abs_pd_array)  # f.append('key', value)

# figure section

# can read in the SASA pandas file, gets the arrays, if more arrays in a file, can include the 'key' from writing
sasa_around_res_abs_array_readin = pd.read_hdf(sasa_around_res_abs_pd_file)
sasa_around_res_rel_array_readin = pd.read_hdf(sasa_around_res_rel_pd_file)
sasa_around_res_ref_abs_array_readin = pd.read_hdf(sasa_around_res_ref_abs_pd_file)


sasa_around_res_abs_diag_array_readin_buffer = np.array(sasa_around_res_abs_array_readin)
sasa_around_res_rel_diag_array_readin_buffer = np.array(sasa_around_res_rel_array_readin)
sasa_around_res_abs_diag_array = []
sasa_around_res_rel_diag_array = []
for x in range(0, 167):
    sasa_around_res_abs_diag_array.append([sasa_around_res_abs_diag_array_readin_buffer[x, x]])
    sasa_around_res_rel_diag_array.append([sasa_around_res_rel_diag_array_readin_buffer[x, x]])

# writing SASA to csv files
with open(sasa_around_res_abs_diag_file, 'a') as f:
    csv.writer(f).writerows(sasa_around_res_abs_diag_array)
with open(sasa_around_res_rel_diag_file, 'a') as f:
    csv.writer(f).writerows(sasa_around_res_rel_diag_array)



# plots a barplot of reference absolute SASA
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=sasa_around_res_ref_abs_array_readin)
ax.set_title("Barplot of absolute reference SASA 8$\AA$ around Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("absolute SASA [$\AA^2$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(sasa_around_res_diag_abs_barplot_fig)


# plots a heatmap of absolute SASA around res by residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(sasa_around_res_abs_array_readin, center=0.0,cmap="PiYG", cbar_kws={'label': 'absolute SASA [$\AA^2$]'})
ax.set_title("Heatmap of absolute SASA 8$\AA$ around Residue by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA to Cys")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(sasa_around_res_abs_heatmap_fig)

# plots a heatmap of relative SASA around res by residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(sasa_around_res_rel_array_readin, center=0.0,cmap="PiYG", cbar_kws={'label': 'relative SASA [$\AA^2$]'})
ax.set_title("Heatmap of relative SASA 8$\AA$ around Residue by Residue")
ax.set_xlabel('Amino acids')
ax.set_ylabel("mutated AA to Cys")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.yticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5))
plt.savefig(sasa_around_res_rel_heatmap_fig)

# plots a barplot of diagonal absolute SASA
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=sasa_around_res_abs_diag_array)
ax.set_title("Barplot of absolute SASA 8$\AA$ around Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("absolute SASA [$\AA^2$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(sasa_around_res_diag_abs_barplot_fig)

# plots a barplot of diagonal relative SASA
fig = plt.figure(dpi=600)
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.barplot(data=sasa_around_res_rel_diag_array)
ax.set_title("Barplot of relative SASA 8$\AA$ around Residue")
ax.set_xlabel('mutated AA')
ax.set_ylabel("relative SASA [$\AA^2$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(sasa_around_res_diag_rel_barplot_fig)
