#!/bin/python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# file naming section
# rms-rmsf files
rmsd_apo_wt1_file = '/.../bachelor-thesis/MD/runs/apo/wildtype-without-ligand/production-run1/rmsd.csv'
rmsf_apo_wt1_file = '/.../bachelor-thesis/MD/runs/apo/wildtype-without-ligand/production-run1/atomic-fluctuations.csv'
rmsfB_apo_wt1_file = '/.../bachelor-thesis/MD/runs/apo/wildtype-without-ligand/production-run1/atomic-fluctuations-bFactor.csv'
rmsd_apo_wt2_file = '/.../bachelor-thesis/MD/runs/apo/wildtype-without-ligand/production-run2/rmsd.csv'
rmsf_apo_wt2_file = '/.../bachelor-thesis/MD/runs/apo/wildtype-without-ligand/production-run2/atomic-fluctuations.csv'
rmsfB_apo_wt2_file = '/.../bachelor-thesis/MD/runs/apo/wildtype-without-ligand/production-run2/atomic-fluctuations-bFactor.csv'
rmsd_apo_lig1_file = '/.../bachelor-thesis/MD/runs/apo/39C-ligand38/production-run1/rmsd.csv'
rmsf_apo_lig1_file = '/.../bachelor-thesis/MD/runs/apo/39C-ligand38/production-run1/atomic-fluctuations.csv'
rmsfB_apo_lig1_file = '/.../bachelor-thesis/MD/runs/apo/39C-ligand38/production-run1/atomic-fluctuations-bFactor.csv'
rmsd_apo_lig2_file = '/.../bachelor-thesis/MD/runs/apo/39C-ligand38/production-run2/rmsd.csv'
rmsf_apo_lig2_file = '/.../bachelor-thesis/MD/runs/apo/39C-ligand38/production-run2/atomic-fluctuations.csv'
rmsfB_apo_lig2_file = '/.../bachelor-thesis/MD/runs/apo/39C-ligand38/production-run2/atomic-fluctuations-bFactor.csv'
rmsd_holo_wt1_file = '/.../bachelor-thesis/MD/runs/holo/wildtype-without-ligand/production-run1/rmsd.csv'
rmsf_holo_wt1_file = '/.../bachelor-thesis/MD/runs/holo/wildtype-without-ligand/production-run1/atomic-fluctuations.csv'
rmsfB_holo_wt1_file = '/.../bachelor-thesis/MD/runs/holo/wildtype-without-ligand/production-run1/atomic-fluctuations-bFactor.csv'
rmsd_holo_wt2_file = '/.../bachelor-thesis/MD/runs/holo/wildtype-without-ligand/production-run2/rmsd.csv'
rmsf_holo_wt2_file = '/.../bachelor-thesis/MD/runs/holo/wildtype-without-ligand/production-run2/atomic-fluctuations.csv'
rmsfB_holo_wt2_file = '/.../bachelor-thesis/MD/runs/holo/wildtype-without-ligand/production-run2/atomic-fluctuations-bFactor.csv'


# figure files
rmsd_apo_wt1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsd_apo_wt1.png'
rmsd_apo_wt_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsd_apo_wt-combined.png'
rmsf_apo_wt1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_apo_wt1.png'
rmsf_apo_wt_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_apo_wt-combined.png'
rmsfB_apo_wt1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_weighted_apo_wt1.png'
rmsfB_apo_wt_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_weighted_apo_wt-combined.png'
rmsd_apo_lig1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsd_apo_lig1.png'
rmsd_apo_lig_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsd_apo_lig-combined.png'
rmsf_apo_lig1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_apo_lig1.png'
rmsf_apo_lig_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_apo_lig-combined.png'
rmsfB_apo_lig1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_weighted_apo_lig1.png'
rmsfB_apo_lig_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_weighted_apo_lig-combined.png'
rmsd_holo_wt1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsd_holo_wt1.png'
rmsd_holo_wt_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsd_holo_wt-combined.png'
rmsf_holo_wt1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_holo_wt1.png'
rmsf_holo_wt_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_holo_wt-combined.png'
rmsfB_holo_wt1_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_weighted_holo_wt1.png'
rmsfB_holo_wt_combined_fig = '/.../bachelor-thesis/MD/runs/rmsd-rmsf/rmsf_weighted_holo_wt-combined.png'


# figure section

# read in rmsd and rmsf files
rmsd_apo_wt1_list = pd.read_csv(rmsd_apo_wt1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsf_apo_wt1_list = pd.read_csv(rmsf_apo_wt1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsfB_apo_wt1_list = pd.read_csv(rmsfB_apo_wt1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)

rmsd_apo_wt2_list = pd.read_csv(rmsd_apo_wt2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsf_apo_wt2_list = pd.read_csv(rmsf_apo_wt2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsfB_apo_wt2_list = pd.read_csv(rmsfB_apo_wt2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)

rmsd_apo_lig1_list = pd.read_csv(rmsd_apo_lig1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsf_apo_lig1_list = pd.read_csv(rmsf_apo_lig1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsfB_apo_lig1_list = pd.read_csv(rmsfB_apo_lig1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)

rmsd_apo_lig2_list = pd.read_csv(rmsd_apo_lig2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsf_apo_lig2_list = pd.read_csv(rmsf_apo_lig2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsfB_apo_lig2_list = pd.read_csv(rmsfB_apo_lig2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)

rmsd_holo_wt1_list = pd.read_csv(rmsd_holo_wt1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsf_holo_wt1_list = pd.read_csv(rmsf_holo_wt1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsfB_holo_wt1_list = pd.read_csv(rmsfB_holo_wt1_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)

rmsd_holo_wt2_list = pd.read_csv(rmsd_holo_wt2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsf_holo_wt2_list = pd.read_csv(rmsf_holo_wt2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)
rmsfB_holo_wt2_list = pd.read_csv(rmsfB_holo_wt2_file, sep=';', skipinitialspace=1, decimal='.', dtype=float)




frames_x_axes = [1, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500, 25000]

# figures


# apo structure without ligand production-run1

# plots a lineplot of rmsd of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsd_apo_wt1_list, legend=False, lw=0.5)
plt.title("Lineplot of RMSD of non-mutated apo structure")
plt.xlabel('Frames')
plt.ylabel("RMSD [$\AA$]")
plt.xticks(np.arange(0, 25001, step=2500), frames_x_axes , rotation=45)
plt.savefig(rmsd_apo_wt1_fig)

# plots a lineplot of rmsf of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsf_apo_wt1_list, legend=False)
plt.title("Lineplot of RMSF of non-mutated apo structure")
plt.xlabel('Amino acids')
plt.ylabel("RMSF [$\AA$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsf_apo_wt1_fig)

# plots a lineplot of weighted rmsf of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsfB_apo_wt1_list, legend=False)
plt.title("Lineplot of weighted RMSF of non-mutated apo structure")
plt.xlabel('Amino acids')
plt.ylabel("weighted RMSF ["r"$\AA$ $\cdot$ $\frac{3}{8}$$\pi$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsfB_apo_wt1_fig)



# apo structure without ligand combined production-runs

# plots a lineplot of rmsd of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsd_apo_wt1_list, lw=0.5, legend=False)
sns.lineplot(data=rmsd_apo_wt2_list, palette=['red'], lw=0.5, legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of RMSD of non-mutated apo structure")
plt.xlabel('Frames')
plt.ylabel("RMSD [$\AA$]")
plt.xticks(np.arange(0, 25001, step=2500), frames_x_axes , rotation=45)
plt.savefig(rmsd_apo_wt_combined_fig)

# plots a lineplot of rmsf of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsf_apo_wt1_list, legend=False)
sns.lineplot(data=rmsf_apo_wt2_list, palette=['red'], legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of RMSF of non-mutated apo structure")
plt.xlabel('Amino acids')
plt.ylabel("RMSF [$\AA$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsf_apo_wt_combined_fig)

# plots a lineplot of weighted rmsf of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsfB_apo_wt1_list, legend=False)
sns.lineplot(data=rmsfB_apo_wt2_list, palette=['red'], legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of weighted RMSF of non-mutated apo structure")
plt.xlabel('Amino acids')
plt.ylabel("weighted RMSF ["r"$\AA$ $\cdot$ $\frac{3}{8}$$\pi$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsfB_apo_wt_combined_fig)




# apo structure 39C with ligand 38 production-run1

# plots a lineplot of rmsd of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsd_apo_lig1_list, legend=False, lw=0.5)
plt.title("Lineplot of RMSD of non-mutated apo structure")
plt.xlabel('Frames')
plt.ylabel("RMSD [$\AA$]")
plt.xticks(np.arange(0, 25001, step=2500), frames_x_axes , rotation=45)
plt.savefig(rmsd_apo_lig1_fig)

# plots a lineplot of rmsf of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsf_apo_lig1_list, legend=False)
plt.title("Lineplot of RMSF of non-mutated apo structure")
plt.xlabel('Amino acids')
plt.ylabel("RMSF [$\AA$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsf_apo_lig1_fig)

# plots a lineplot of weighted rmsf of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsfB_apo_lig1_list, legend=False)
plt.title("Lineplot of weighted RMSF of non-mutated apo structure")
plt.xlabel('Amino acids')
plt.ylabel("weighted RMSF ["r"$\AA$ $\cdot$ $\frac{3}{8}$$\pi$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsfB_apo_lig1_fig)


# apo structure 39C with ligand 38 combined production-runs

# plots a lineplot of rmsd of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsd_apo_lig1_list, lw=0.5, legend=False)
sns.lineplot(data=rmsd_apo_lig2_list, palette=['red'], lw=0.5, legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of RMSD of non-mutated apo structure")
plt.xlabel('Frames')
plt.ylabel("RMSD [$\AA$]")
plt.xticks(np.arange(0, 25001, step=2500), frames_x_axes , rotation=45)
plt.savefig(rmsd_apo_lig_combined_fig)

# plots a lineplot of rmsf of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsf_apo_lig1_list, legend=False)
sns.lineplot(data=rmsf_apo_lig2_list, palette=['red'], legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of RMSF of non-mutated apo structure")
plt.xlabel('Amino acids')
plt.ylabel("RMSF [$\AA$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsf_apo_lig_combined_fig)

# plots a lineplot of weighted rmsf of apo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsfB_apo_lig1_list, legend=False)
sns.lineplot(data=rmsfB_apo_lig2_list, palette=['red'], legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of weighted RMSF of non-mutated apo structure")
plt.xlabel('Amino acids')
plt.ylabel("weighted RMSF ["r"$\AA$ $\cdot$ $\frac{3}{8}$$\pi$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsfB_apo_lig_combined_fig)




# holo structure without ligand production-run1

# plots a lineplot of rmsd of holo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsd_holo_wt1_list, legend=False, lw=0.5)
plt.title("Lineplot of RMSD of non-mutated holo structure")
plt.xlabel('Frames')
plt.ylabel("RMSD [$\AA$]")
plt.xticks(np.arange(0, 25001, step=2500), frames_x_axes , rotation=45)
plt.savefig(rmsd_holo_wt1_fig)

# plots a lineplot of rmsf of holo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsf_holo_wt1_list, legend=False)
plt.title("Lineplot of RMSF of non-mutated holo structure")
plt.xlabel('Amino acids')
plt.ylabel("RMSF [$\AA$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsf_holo_wt1_fig)

# plots a lineplot of weighted rmsf of holo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsfB_holo_wt1_list, legend=False)
plt.title("Lineplot of weighted RMSF of non-mutated holo structure")
plt.xlabel('Amino acids')
plt.ylabel("weighted RMSF ["r"$\AA$ $\cdot$ $\frac{3}{8}$$\pi$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsfB_holo_wt1_fig)



# holo structure without ligand combined production-runs

# plots a lineplot of rmsd of holo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
sns.lineplot(data=rmsd_holo_wt1_list, lw=0.5, legend=False)
sns.lineplot(data=rmsd_holo_wt2_list, palette=['red'], lw=0.5, legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of RMSD of non-mutated holo structure")
plt.xlabel('Frames')
plt.ylabel("RMSD [$\AA$]")
plt.xticks(np.arange(0, 25001, step=2500), frames_x_axes , rotation=45)
plt.savefig(rmsd_holo_wt_combined_fig)

# plots a lineplot of rmsf of holo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsf_holo_wt1_list, legend=False)
sns.lineplot(data=rmsf_holo_wt2_list, palette=['red'], legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of RMSF of non-mutated holo structure")
plt.xlabel('Amino acids')
plt.ylabel("RMSF [$\AA$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsf_holo_wt_combined_fig)

# plots a lineplot of weighted rmsf of holo structure, without ligand
fig = plt.figure(dpi=600, figsize=(6, 6))
sns.set('paper', 'white', rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.titlesize': 10,
                              'xtick.labelsize': 8,
                              'ytick.labelsize': 8, "pgf.rcfonts": False})  # can also be sns.set()
plt.rc('text', usetex=True)
ax = sns.lineplot(data=rmsfB_holo_wt1_list, legend=False)
sns.lineplot(data=rmsfB_holo_wt2_list, palette=['red'], legend=False)
plt.legend(labels=['production-run1','production-run2'], facecolor='white')
plt.title("Lineplot of weighted RMSF of non-mutated holo structure")
plt.xlabel('Amino acids')
plt.ylabel("weighted RMSF ["r"$\AA$ $\cdot$ $\frac{3}{8}$$\pi$]")
plt.xticks(np.arange(0, 167, step=5), np.arange(1, 170, step=5), rotation=90)
plt.savefig(rmsfB_holo_wt_combined_fig)

