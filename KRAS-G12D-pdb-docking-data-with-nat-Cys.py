#!/bin/python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv

sel_AA = [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154]
docking_data_csv_file = "/.../bachelor-thesis/Docking/Docking_scoring-with-natural-cys118.csv"
docking_data_csv_transposed_file = "/.../bachelor-thesis/Docking/Docking_scoring-for-boxplot-with-natural-cys118.csv"
docking_data_txt_file = "/.../bachelor-thesis/Docking/Docking_scoring-with-natural-cys118.txt"
docking_data_pd_file = "/.../bachelor-thesis/Docking/Docking_scoring-with-natural-cys118.h5"
docking_heatmap_file = "/.../bachelor-thesis/Docking/Docking_heatmap-with-natural-cys118.png"
sel_AA_data_file = "/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/KRAS-G12D-pdb5xco-data-for-selected-residues-with-natural-cys118.csv"
docking_sasa_figure_file = "/.../bachelor-thesis/Docking/Docking-sasa_compare_figure-with-natural-cys118.png"
docking_sasa_around_res_figure_file = "/.../bachelor-thesis/Docking/Docking-sasa-around-res_compare_figure-with-natural-cys118.png"
docking_both_sasa_figure_file = "/.../bachelor-thesis/Docking/Docking-sasa2-res_compare_figure-with-natural-cys118.png"
ligand_scores_boxplot_fig = "/.../bachelor-thesis/Docking/ligands-score-boxplot-with-natural-cys118.png"
residue_scores_boxplot_fig = "/.../bachelor-thesis/Docking/residue-score-boxplot-with-natural-cys118.png"
all_data_merged_figure_file = "/.../bachelor-thesis/Docking/Docking-all-data_compare_figure-with-natural-cys118.png"

# corrected files
heavy_atoms_file = "/.../bachelor-thesis/Docking/ligands/ligands_heavy_atoms.csv"
docking_data_corr_csv_file = "/.../bachelor-thesis/Docking/Docking_scoring_corrected-with-natural-cys118.csv"
docking_data_corr_csv_transposed_file = "/.../bachelor-thesis/Docking/Docking_scoring_corrected-for-boxplot-with-natural-cys118.csv"
docking_data_corr_txt_file = "/.../bachelor-thesis/Docking/Docking_scoring_corrected-with-natural-cys118.txt"
docking_data_corr_pd_file = "/.../bachelor-thesis/Docking/Docking_scoring_corrected-with-natural-cys118.h5"
docking_corr_heatmap_file = "/.../bachelor-thesis/Docking/Docking_heatmap_corrected-with-natural-cys118.png"
docking_sasa_corr_figure_file = "/.../bachelor-thesis/Docking/Docking-sasa_compare_corrected_figure-with-natural-cys118.png"
docking_sasa_around_res_corr_figure_file = "/.../bachelor-thesis/Docking/Docking-sasa-around-res_compare_corrected_figure-with-natural-cys118.png"
docking_both_sasa_corr_figure_file = "/.../bachelor-thesis/Docking/Docking-sasa2-res_compare_corrected_figure-with-natural-cys118.png"
ligand_scores_corr_boxplot_fig = "/.../bachelor-thesis/Docking/ligands-score_corrected-boxplot-with-natural-cys118.png"
residue_scores_corr_boxplot_fig = "/.../bachelor-thesis/Docking/residue-score_corrected-boxplot-with-natural-cys118.png"


data = np.zeros((56, 18))  # make scoring array for heatmap

# list which is the overall number of the ligand in the broads library
alpha = [7, 13, 17, 18, 20, 23, 25, 26, 28, 29, 31, 36, 38, 40, 41, 47, 50, 51]
beta = [1, 2, 3, 4, 5, 6, 11, 12, 15, 16, 21, 22, 32, 34, 37, 39, 42, 43, 44, 45, 48, 49, 53, 54, 55, 56]
disulfid = [52]
nitril = [8, 9, 10, 14, 19, 24, 27, 30, 33, 35]
nitril_plus = [46]


t = 0
# loop through all directories
for i in sel_AA:
    # file directories

    # alpha Ketone
    # get data out of file
    src = "/.../bachelor-thesis/Docking/AS" + str(i) + "/alphaKeton/output/bestranking.lst"
    names = np.loadtxt(src, dtype=str, usecols=10)
    score_all = np.loadtxt(src, usecols=0)

    # look which ligand has which score and write to big array
    for x in range(1, 19):
        sub = str("m" + str(x) + "_")
        index = -1
        row = 0  # indicates the row
        for s in names:
            if sub in s:
                index = row  # if found then write row to index
                break
            row += 1  # if not found increases row by one
        if index >= 0:
            score = score_all[index]
        else:
            score = -100  # used to test if mx_ is found or not
        index_all = alpha[x - 1]
        data[index_all - 1, t] = score

    # beta Ketone
    # get data out of file
    src = "/.../bachelor-thesis/Docking/AS" + str(i) + "/betaKeton/output/bestranking.lst"
    names = np.loadtxt(src, dtype=str, usecols=10)
    score_all = np.loadtxt(src, usecols=0)

    # look which ligand has which score and write to big array
    for x in range(1, 28):
        sub = str("m" + str(x) + "_")
        index = -1
        row = 0
        for s in names:
            if sub in s:
                index = row
                break
            row += 1
        if index >= 0:
            score = score_all[index]
        else:
            score = -100  # used to test if mx_ is found or not

        if x < 14:
            index_all = beta[x - 1]
            data[index_all - 1, t] = score

        if x > 14:
            index_all = beta[x - 2]
            data[index_all - 1, t] = score

    # nitril
    # get data out of file
    src = "/.../bachelor-thesis/Docking/AS" + str(i) + "/nitril/output/bestranking.lst"
    names = np.loadtxt(src, dtype=str, usecols=10)
    score_all = np.loadtxt(src, usecols=0)

    # look which ligand has which score and write to big array
    for x in range(1, 11):
        sub = str("m" + str(x) + "_")
        index = -1
        row = 0
        for s in names:
            if sub in s:
                index = row
                break
            row += 1
        if index >= 0:
            score = score_all[index]
        else:
            score = -100  # used to test if mx_ is found or not
        index_all = nitril[x - 1]
        data[index_all - 1, t] = score

    # nitril+
    # get data out of file
    src = "/.../bachelor-thesis/Docking/AS" + str(i) + "/nitril+/output/bestranking.lst"
    names = np.loadtxt(src, dtype=str, usecols=10)
    score_all = np.loadtxt(src, usecols=0)

    # look which ligand has which score and write to big array
    score = score_all
    index_all = nitril_plus[0]
    data[index_all - 1, t] = score

    # disulfid
    # get data out of file
    src = "/.../bachelor-thesis/Docking/AS" + str(i) + "/disulfid/output/bestranking.lst"
    names = np.loadtxt(src, dtype=str, usecols=10)
    score_all = np.loadtxt(src, usecols=0)

    # look which ligand has which score and write to big array
    score = score_all
    index_all = disulfid[0]
    data[index_all - 1, t] = score

    t += 1

# correct the scoring depending on the heavy atoms
data_corrected = np.zeros((56, 18))

csv_file_data = pd.read_csv(heavy_atoms_file, sep=';')
heavy_atoms_list = csv_file_data['heavy atoms'].tolist()

for col in range(0, 56):
    for value in range(0, 18):
        data_corrected[col, value] = data[col, value] / heavy_atoms_list[col]

# write scores to file
with open(docking_data_csv_file, 'a') as f:
    csv.writer(f).writerows(data)
data_trans = np.transpose(data)
with open(docking_data_csv_transposed_file, 'a') as f:
    csv.writer(f).writerows(data_trans)
with open(docking_data_txt_file, 'a') as f:
    f.write(str(data))
docking_pd_array = pd.DataFrame(data)  # .transpose
with pd.HDFStore(docking_data_pd_file, 'w') as f:
    f.append('docking_pd_array', docking_pd_array)  # f.append('key', value)

# write corrected scores to file
with open(docking_data_corr_csv_file, 'a') as f:
    csv.writer(f).writerows(data_corrected)
data_corr_trans = np.transpose(data_corrected)
with open(docking_data_corr_csv_transposed_file, 'a') as f:
    csv.writer(f).writerows(data_corr_trans)
with open(docking_data_corr_txt_file, 'a') as f:
    f.write(str(data_corrected))
docking_corr_pd_array = pd.DataFrame(data_corrected)  # .transpose
with pd.HDFStore(docking_data_corr_pd_file, 'w') as f:
    f.append('docking_pd_array', docking_corr_pd_array)  # f.append('key', value)

# get data for the residues out of a csv file and put it into a figure later
csv_file_data = pd.read_csv(sel_AA_data_file, sep=';')
sasa_list = csv_file_data['abs SASA by residue after mut'].tolist()
sasa_around_res_list = csv_file_data['abs SASA 8A around res by residue after mut'].tolist()
energy_pot_rel_list = csv_file_data['rel pot Energy by residue after mutation'].tolist()
stab_list = csv_file_data['delta Stability'].tolist()
pka_list_buffer = csv_file_data['pka of residue after mutation'].tolist()
pka_list = []

for p in pka_list_buffer:
    pka_list.append(p - 7.5)

# plots a heatmap of scoring of the ligands by residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(data=data, center=0.0, cmap="PiYG", cbar_kws={'label': 'Docking Score'})
ax.set_title("Heatmap of Docking Scores with natural Cys118")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Ligands")
plt.xticks(np.arange(0, 18, step=1) + 0.5, [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154],
           rotation=90)
plt.yticks(np.arange(0, 57, step=3) + 0.5, np.arange(1, 60, step=3), rotation=360)
plt.savefig(docking_heatmap_file)

# plots a heatmap with a barplot with sasa underneath
sns.set('paper', 'white', rc={'font.size': 20, 'axes.labelsize': 14, 'legend.fontsize': 12, 'axes.titlesize': 24,
                              'xtick.labelsize': 14,'ytick.labelsize': 14, "pgf.rcfonts": False})  # can also be sns.set()
fig, ((ax1, cbar_ax), (ax2, dummy_ax)) = plt.subplots(nrows=2, ncols=2, sharex='col', 
							dpi=600, figsize=(16, 12),
                                                      	gridspec_kw={'height_ratios': [5, 1], 
							'width_ratios': [20, 1]})

sns.heatmap(data=data, center=0.0, cmap="PiYG", cbar_ax=cbar_ax, xticklabels=False, ax=ax1,
            cbar_kws={'label': 'Docking Score'})

ax1.set_title("Heatmap of Docking Scores with natural Cys118 and absolute SASA of mutated Residues")
ax1.set_ylabel("Ligands")
ax2.set_ylabel("absolute SASA [$\AA^2$]")
ax2.set_xlabel('selected AA')

y_ticks = np.arange(1, 60, step=3)
y_tick_pos = [i + 0.3 for i in np.arange(0, 57, step=3)]
ax1.set_yticks(y_tick_pos)
ax1.set_yticklabels(y_ticks, rotation=360)
# ax1.set_yticks(np.arange(0, 57, step=3)+0.5, np.arange(1, 60, step=3))
x_ticks = [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154]
x_tick_pos = [i + 0.5 for i in range(len(x_ticks))]
ax2.bar(x_tick_pos, sasa_list, align='center')
ax2.set_xticks(x_tick_pos)  # maybe ,rotation=90)
ax2.set_xticklabels(x_ticks, rotation=90)
dummy_ax.axis('off')

plt.tight_layout()
plt.savefig(docking_sasa_figure_file)

# plots a heatmap with a barplot with sasa around res underneath
sns.set('paper', 'white', rc={'font.size': 20, 'axes.labelsize': 14, 'legend.fontsize': 12, 'axes.titlesize': 24,
                              'xtick.labelsize': 14,
                              'ytick.labelsize': 14, "pgf.rcfonts": False})  # can also be sns.set()
fig, ((ax1, cbar_ax), (ax2, dummy_ax)) = plt.subplots(nrows=2, ncols=2, sharex='col', 
							dpi=600, figsize=(16, 12),
                                                     	gridspec_kw={'height_ratios': [5, 1],
							'width_ratios': [20, 1]})

sns.heatmap(data=data, center=0.0, cmap="PiYG", cbar_ax=cbar_ax, xticklabels=False, ax=ax1,
            cbar_kws={'label': 'Docking Score'})

ax1.set_title("Heatmap of Docking Scores with natural Cys118 and absolute SASA 8$\AA$ around Residue of mutated Residues")
ax1.set_ylabel("Ligands")
ax2.set_ylabel("abs. SASA 8$\AA$ around Residue [$\AA^2$]")
ax2.set_xlabel('selected AA')

y_ticks = np.arange(1, 60, step=3)
y_tick_pos = [i + 0.3 for i in np.arange(0, 57, step=3)]
ax1.set_yticks(y_tick_pos)
ax1.set_yticklabels(y_ticks, rotation=360)
# ax1.set_yticks(np.arange(0, 57, step=3)+0.5, np.arange(1, 60, step=3))
x_ticks = [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154]
x_tick_pos = [i + 0.5 for i in range(len(x_ticks))]
ax2.bar(x_tick_pos, sasa_around_res_list, align='center')
ax2.set_xticks(x_tick_pos)  # maybe ,rotation=90)
ax2.set_xticklabels(x_ticks, rotation=90)
dummy_ax.axis('off')

plt.tight_layout()
plt.savefig(docking_sasa_around_res_figure_file)

# plots a heatmap with a barplot of abs sasa and abs sasa around res underneath
sns.set('paper', 'white', rc={'font.size': 20, 'axes.labelsize': 16, 'legend.fontsize': 18, 'axes.titlesize': 28,
                              'xtick.labelsize': 16,
                              'ytick.labelsize': 16, "pgf.rcfonts": False})  # can also be sns.set()
fig, ((ax1, cbar_ax), (ax2, dummy_ax2),	
	(ax3, dummy_ax3)) = plt.subplots(nrows=3, ncols=2, sharex='col', dpi=600,
					figsize=(20, 16), gridspec_kw={'height_ratios': [5, 1.5, 1.5],
									'width_ratios': [20, 1]})

sns.heatmap(data=data, center=0.0, cmap="PiYG", cbar_ax=cbar_ax, xticklabels=False, ax=ax1,
            cbar_kws={'label': 'Docking Score'})

ax1.set_title("Heatmap of Docking Scores with natural Cys118 and absolute SASA of mutated Residues")
ax1.set_ylabel("Ligands")
ax2.set_ylabel("absolute SASA [$\AA^2$]")
ax3.set_ylabel("abs. SASA 8$\AA$ around Residue [$\AA^2$]")
ax3.set_xlabel('selected AA')

y_ticks = np.arange(1, 60, step=3)
y_tick_pos = [i + 0.3 for i in np.arange(0, 57, step=3)]
ax1.set_yticks(y_tick_pos)
ax1.set_yticklabels(y_ticks, rotation=360)
x_ticks = [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154]
x_tick_pos = [i + 0.5 for i in range(len(x_ticks))]
y_ticks_pka = [7.5, 8, 8.5, 9, 9.5, 10]
y_tick_pos_pka = [0, 0.5, 1, 1.5, 2, 2.5]
ax2.bar(x_tick_pos, sasa_list, align='center')
ax3.bar(x_tick_pos, sasa_around_res_list, align='center')
ax3.set_xticks(x_tick_pos)  # maybe ,rotation=90)
ax3.set_xticklabels(x_ticks, rotation=90)
dummy_ax2.axis('off')
dummy_ax3.axis('off')

plt.tight_layout()
plt.savefig(docking_both_sasa_figure_file)

# boxplot of docking scores of ligands
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.boxplot(data=data_trans)  # , cmap="PiYG")
ax.set_title("Boxplot of Docking Scores with natural Cys118 of the Ligands")
ax.set_xlabel('Ligands')
ax.set_ylabel("Docking Score")
plt.xticks(np.arange(0, 56, step=2), np.arange(1, 57, step=2), rotation=90)
# plt.yticks(np.arange(0, 17, step=3)+0.5, np.arange(1, 60, step=3), rotation=360)
plt.savefig(ligand_scores_boxplot_fig)

# boxplot of docking scores of selectes residues
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.boxplot(data=data)  # , cmap="PiYG")
ax.set_title("Boxplot of Docking Scores with natural Cys118 of the selected Residues")
ax.set_xlabel('Amino acids')
ax.set_ylabel("Docking Score")
plt.xticks(np.arange(0, 18, step=1), [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154],
           rotation=90)
# plt.yticks(np.arange(0, 17, step=3)+0.5, np.arange(1, 60, step=3), rotation=360)
plt.savefig(residue_scores_boxplot_fig)

# plots a heatmap with a barplot of abs sasa, pot rel energy, stability and pka underneath
sns.set('paper', 'white', rc={'font.size': 20, 'axes.labelsize': 16, 'legend.fontsize': 18, 'axes.titlesize': 28,
                              'xtick.labelsize': 16,
                              'ytick.labelsize': 16, "pgf.rcfonts": False})  # can also be sns.set()
fig, ((ax1, cbar_ax), (ax2, dummy_ax2), (ax3, dummy_ax3), 
	(ax4, dummy_ax4), (ax5, dummy_ax5)) = plt.subplots(nrows=5, ncols=2, sharex='col', 
							dpi=400, figsize=(20,16), 
							gridspec_kw={'height_ratios': [5, 1.5, 1.5, 1.5, 1.5],
									'width_ratios': [ 20, 1]})

sns.heatmap(data=data, center=0.0, cmap="PiYG", cbar_ax=cbar_ax, xticklabels=False, ax=ax1,
            cbar_kws={'label': 'Docking Score'})

ax1.set_title("Heatmap of Docking Scores with natural Cys118 and absolute SASA of mutated Residues")
ax1.set_ylabel("Ligands")
ax2.set_ylabel("absolute SASA [$\AA^2$]")
ax3.set_ylabel("rel pot Energy [kcal/mol]")
ax4.set_ylabel("$\Delta$Stability")
ax5.set_ylabel("pKa")
ax5.set_xlabel('selected AA')

y_ticks = np.arange(1, 60, step=3)
y_tick_pos = [i + 0.3 for i in np.arange(0, 57, step=3)]
ax1.set_yticks(y_tick_pos)
ax1.set_yticklabels(y_ticks, rotation=360)
x_ticks = [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154]
x_tick_pos = [i + 0.5 for i in range(len(x_ticks))]
y_ticks_pka = [7.5, 8, 8.5, 9, 9.5, 10]
y_tick_pos_pka = [0, 0.5, 1, 1.5, 2, 2.5]
ax2.bar(x_tick_pos, sasa_list, align='center')
ax3.bar(x_tick_pos, energy_pot_rel_list, align='center')
ax4.bar(x_tick_pos, stab_list, align='center')
ax5.bar(x_tick_pos, pka_list)
ax5.set_xticks(x_tick_pos)  # maybe ,rotation=90)
ax5.set_xticklabels(x_ticks, rotation=90)
ax5.set_yticks(y_tick_pos_pka)
ax5.set_yticklabels(y_ticks_pka, rotation=360)
dummy_ax2.axis('off')
dummy_ax3.axis('off')
dummy_ax4.axis('off')
dummy_ax5.axis('off')

plt.tight_layout()
plt.savefig(all_data_merged_figure_file)

# figures of corrected value

# plots a heatmap of corrected scoring of the ligands by residue
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.heatmap(data=data_corrected, center=0.0, cmap="PiYG", cbar_kws={'label': 'Docking Score'})
ax.set_title("Heatmap of corrected Docking Scores with natural Cys118")
ax.set_xlabel('mutated AA')
ax.set_ylabel("Ligands")
plt.xticks(np.arange(0, 18, step=1) + 0.5, [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154],
           rotation=90)
plt.yticks(np.arange(0, 57, step=3) + 0.5, np.arange(1, 60, step=3), rotation=360)
plt.savefig(docking_corr_heatmap_file)

# plots a corrected heatmap with a barplot with sasa underneath
sns.set('paper', 'white', rc={'font.size': 20, 'axes.labelsize': 14, 'legend.fontsize': 12, 'axes.titlesize': 24,
                              'xtick.labelsize': 14,
                              'ytick.labelsize': 14, "pgf.rcfonts": False})  # can also be sns.set()
fig, ((ax1, cbar_ax), (ax2, dummy_ax)) = plt.subplots(nrows=2, ncols=2, sharex='col', 
							dpi=600, figsize=(16, 12),
                                                      	gridspec_kw={'height_ratios': [5, 1], 
									'width_ratios': [20, 1]})

sns.heatmap(data=data_corrected, center=0.0, cmap="PiYG", cbar_ax=cbar_ax, xticklabels=False, ax=ax1,
            cbar_kws={'label': 'Docking Score'})

ax1.set_title("Heatmap of corrected Docking Scores with natural Cys118 and absolute SASA of mutated Residues")
ax1.set_ylabel("Ligands")
ax2.set_ylabel("absolute SASA [$\AA^2$]")
ax2.set_xlabel('selected AA')

y_ticks = np.arange(1, 60, step=3)
y_tick_pos = [i + 0.3 for i in np.arange(0, 57, step=3)]
ax1.set_yticks(y_tick_pos)
ax1.set_yticklabels(y_ticks, rotation=360)
# ax1.set_yticks(np.arange(0, 57, step=3)+0.5, np.arange(1, 60, step=3))
x_ticks = [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154]
x_tick_pos = [i + 0.5 for i in range(len(x_ticks))]
ax2.bar(x_tick_pos, sasa_list, align='center')
ax2.set_xticks(x_tick_pos)  # maybe ,rotation=90)
ax2.set_xticklabels(x_ticks, rotation=90)
dummy_ax.axis('off')

plt.tight_layout()
plt.savefig(docking_sasa_corr_figure_file)

# plots a corrected heatmap with a barplot with sasa around residue underneath
sns.set('paper', 'white', rc={'font.size': 20, 'axes.labelsize': 14, 'legend.fontsize': 12, 'axes.titlesize': 24,
                              'xtick.labelsize': 14,
                              'ytick.labelsize': 14, "pgf.rcfonts": False})  # can also be sns.set()
fig, ((ax1, cbar_ax), (ax2, dummy_ax)) = plt.subplots(nrows=2, ncols=2, sharex='col', 
							dpi=600, figsize=(16, 12),
                                                      	gridspec_kw={'height_ratios': [5, 1], 
									'width_ratios': [20, 1]})

sns.heatmap(data=data_corrected, center=0.0, cmap="PiYG", cbar_ax=cbar_ax, xticklabels=False, ax=ax1,
            cbar_kws={'label': 'Docking Score'})

ax1.set_title("Heatmap of corrected Docking Scores with natural Cys118 and absolute SASA 8$\AA$ around Residue of mutated Residues")
ax1.set_ylabel("Ligands")
ax2.set_ylabel("abs. SASA 8$\AA$ around Residue [$\AA^2$]")
ax2.set_xlabel('selected AA')

y_ticks = np.arange(1, 60, step=3)
y_tick_pos = [i + 0.3 for i in np.arange(0, 57, step=3)]
ax1.set_yticks(y_tick_pos)
ax1.set_yticklabels(y_ticks, rotation=360)
# ax1.set_yticks(np.arange(0, 57, step=3)+0.5, np.arange(1, 60, step=3))
x_ticks = [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154]
x_tick_pos = [i + 0.5 for i in range(len(x_ticks))]
ax2.bar(x_tick_pos, sasa_around_res_list, align='center')
ax2.set_xticks(x_tick_pos)  # maybe ,rotation=90)
ax2.set_xticklabels(x_ticks, rotation=90)
dummy_ax.axis('off')

plt.tight_layout()
plt.savefig(docking_sasa_around_res_corr_figure_file)

# plots a corrected heatmap with a barplot of abs sasa and abs sasa around res underneath
sns.set('paper', 'white', rc={'font.size': 20, 'axes.labelsize': 16, 'legend.fontsize': 18, 'axes.titlesize': 28,
                              'xtick.labelsize': 16,
                              'ytick.labelsize': 16, "pgf.rcfonts": False})  # can also be sns.set()
fig, ((ax1, cbar_ax), (ax2, dummy_ax2), (ax3, dummy_ax3)) = plt.subplots(nrows=3, ncols=2, sharex='col', 
									dpi=600, figsize=(20, 16),
									gridspec_kw={'height_ratios': [5, 1.5, 1.5],
                                                                                      'width_ratios': [20, 1]})

sns.heatmap(data=data_corrected, center=0.0, cmap="PiYG", cbar_ax=cbar_ax, xticklabels=False, ax=ax1,
            cbar_kws={'label': 'Docking Score'})

ax1.set_title("Heatmap of corrected Docking Scores with natural Cys118 and absolute SASA of mutated Residues")
ax1.set_ylabel("Ligands")
ax2.set_ylabel("absolute SASA [$\AA^2$]")
ax3.set_ylabel("abs. SASA 8$\AA$ around Residue [$\AA^2$]")
ax3.set_xlabel('selected AA')

y_ticks = np.arange(1, 60, step=3)
y_tick_pos = [i + 0.3 for i in np.arange(0, 57, step=3)]
ax1.set_yticks(y_tick_pos)
ax1.set_yticklabels(y_ticks, rotation=360)
x_ticks = [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154]
x_tick_pos = [i + 0.5 for i in range(len(x_ticks))]
y_ticks_pka = [7.5, 8, 8.5, 9, 9.5, 10]
y_tick_pos_pka = [0, 0.5, 1, 1.5, 2, 2.5]
ax2.bar(x_tick_pos, sasa_list, align='center')
ax3.bar(x_tick_pos, sasa_around_res_list, align='center')
ax3.set_xticks(x_tick_pos)  # maybe ,rotation=90)
ax3.set_xticklabels(x_ticks, rotation=90)
dummy_ax2.axis('off')
dummy_ax3.axis('off')

plt.tight_layout()
plt.savefig(docking_both_sasa_corr_figure_file)

# boxplot of corrected docking scores of ligands
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.boxplot(data=data_corr_trans)  # , cmap="PiYG")
ax.set_title("Boxplot of corrected Docking Scores with natural Cys118 of the Ligands")
ax.set_xlabel('Ligands')
ax.set_ylabel("Docking Score")
plt.xticks(np.arange(0, 56, step=2), np.arange(1, 57, step=2), rotation=90)
# plt.yticks(np.arange(0, 17, step=3)+0.5, np.arange(1, 60, step=3), rotation=360)
plt.savefig(ligand_scores_corr_boxplot_fig)

# boxplot of corrected docking scores of selectes residues
fig = plt.figure(dpi=600)  # necessary to reset figure, otherwise will have some parts of figures before in new figure
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.boxplot(data=data_corrected)  # , cmap="PiYG")
ax.set_title("Boxplot of corrected Docking Scores with natural Cys118 of the selected Residues")
ax.set_xlabel('Amino acids')
ax.set_ylabel("Docking Score")
plt.xticks(np.arange(0, 18, step=1), [2, 12, 25, 38, 39, 45, 62, 64, 74, 86, 98, 103, 108, 118, 129, 135, 148, 154],
           rotation=90)
# plt.yticks(np.arange(0, 17, step=3)+0.5, np.arange(1, 60, step=3), rotation=360)
plt.savefig(residue_scores_corr_boxplot_fig)
