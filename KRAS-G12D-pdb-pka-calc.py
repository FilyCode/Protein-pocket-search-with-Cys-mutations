#!/bin/python3

from schrodinger import structure
from schrodinger.structutils import build
from schrodinger.forcefield import minimizer
from schrodinger.application.bioluminate import protein
import numpy as np
import pandas as pd


# file naming section

# naming all files
# read in structure and read out structure
fname = '/.../bachelor-thesis/starting-structures/KRAS-G12D-pdb5xco_1_A_prepared_minimized.mae'  # name of file to load
modified_fname = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_mutation_output.mae'  # name of exported file

# pka files
pka_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-abs1.txt'
pka_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-rel1.txt'
pka_res_abs_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-abs1.txt'
pka_res_rel_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-rel1.txt'
pka_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-abs1.h5'
pka_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-rel1.h5'
pka_res_abs_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-abs1.h5'
pka_res_rel_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-rel1.h5'
pka_ref_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-ref1.txt'
pka_res_ref_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-ref1.txt'
pka_ref_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-ref1.h5'
pka_res_ref_pd_file = '/.../bachelor-thesis/data/KRAS-G12D-pdb5xco_1_A/pka/KRAS-G12D-pdb5xco_pka-res-ref1.h5'



# calculate reference values section

# creating reader and writer
st = structure.StructureReader.read(fname)
writer = structure.StructureWriter(modified_fname)
st.title = str("KRAS-G12D-pdb5xco_reference-not-minimized")
writer.append(st)  # write starting structure to file

# minimize wildtype and make starting reference calculation for pka, energy and stability
minimizer.minimize_structure(st)  # minimize wildtype
st.title = str("KRAS-G12D-pdb5xco_reference")
writer.append(st)  # write minimized reference structure to file
pka_ref = float(protein.PropertyCalculator(st, "pka_ref_calc").getTotalpKa() or 0)
pka_res_ref = []
for res in st.residue:
    pka_res_ref.append(float(protein.PropertyCalculator(st, "pka_res_ref_calc").getResiduepKa(res) or 0))

if len(pka_res_ref) > 167:  # remove data in list after aminoacid, like GDP and Mg2+
    len_dif = len(pka_res_ref) - 167
    for l in range(0, len_dif):
        pka_res_ref.pop()

i = int(0)  # for loop counting

# mutation calculation section

# looping through the protein by residue and then mutate, minimize, calculate pka and
# Energy (potential, prime, interaction)
for res in st.residue:  # iterates over all residues of the structure
    i += 1  # for numbering residues and stopping while testing

    if i == 3:  # i == 168, to stop after all aminoacids
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
        pka = float(protein.PropertyCalculator(work_st, "pka_ref_calc").getTotalpKa() or 0)
        pka_abs_array = np.array([pka])
        pka_rel_array = np.array([pka - pka_ref])
        pka_res = []
        for res in work_st.residue:
            pka_res.append(float(protein.PropertyCalculator(work_st, "pka_res_calc").getResiduepKa(res) or 0))
        if len(pka_res) > 167:
            len_dif = len(pka_res) - 167
            for l in range(0, len_dif):
                pka_res.pop()
        pka_res_abs_array = np.array([pka_res])  # absolute pka values by residues

        pka_buffer = []
        for x in range(0, 167):  # loops through list, calculate relative changes in pka by residue
            pka_buffer.append(pka_res[x] - pka_res_ref[x])
        pka_res_rel_array = np.array([pka_buffer])


    else:
        # pka calculations and writing to arrays
        pka = float(protein.PropertyCalculator(work_st, "pka_ref_calc").getTotalpKa() or 0)
        pka_abs_array = np.append(pka_abs_array, [pka], axis=0)  # absolute pka
        pka_rel_array = np.append(pka_rel_array, [pka - pka_ref], axis=0)  # relative pka
        pka_res = []
        for res in work_st.residue:
            pka_res.append(float(protein.PropertyCalculator(work_st, "pka_res_calc").getResiduepKa(res) or 0))
        if len(pka_res) > 167:
            len_dif = len(pka_res) - 167
            for l in range(0, len_dif):
                pka_res.pop()
        pka_res_abs_array = np.append(pka_res_abs_array, [pka_res], axis=0)

        pka_buffer = []
        for x in range(0, 167):  # loops through list, calculate relative changes in pka by residue
            pka_buffer.append(pka_res[x] - pka_res_ref[x])
        pka_res_rel_array = np.append(pka_res_rel_array, [pka_buffer], axis=0)



# writing to file section

# writing pka to txt files
with open(pka_abs_file, 'a') as f:
    f.write("pka_mut-absolut:\n" + str(pka_abs_array) + "\n")
with open(pka_rel_file, 'a') as f:
    f.write("pka_mut-relative:\n" + str(pka_rel_array) + "\n")
with open(pka_res_abs_file, 'a') as f:
    f.write("pka_mut-residue-absolute:\n" + str(pka_res_abs_array) + "\n")
with open(pka_res_rel_file, 'a') as f:
    f.write("pka_mut-residue-relative:\n" + str(pka_res_rel_array) + "\n")
with open(pka_ref_file, 'a') as f:
    f.write("pka_ref-absolute:\n" + str(pka_ref) + "\n")
with open(pka_res_ref_file, 'a') as f:
    f.write("pka_ref-residue-absolute:\n" + str(pka_res_ref) + "\n")

# creating pandas array of pka then writing to file
pka_abs_pd_array = pd.DataFrame(pka_abs_array).transpose()
pka_rel_pd_array = pd.DataFrame(pka_rel_array).transpose()
pka_res_abs_pd_array = pd.DataFrame(pka_res_abs_array)
pka_res_rel_pd_array = pd.DataFrame(pka_res_rel_array)
pka_ref_pd_array = pd.Series(pka_ref)
pka_res_ref_pd_array = pd.DataFrame(pka_res_ref)
with pd.HDFStore(pka_abs_pd_file, 'w') as f:
    f.append('pka_abs_pd_array', pka_abs_pd_array)  # f.append('key', value)
with pd.HDFStore(pka_rel_pd_file, 'w') as f:
    f.append('pka_rel_pd_array', pka_rel_pd_array)
with pd.HDFStore(pka_res_abs_pd_file, 'w') as f:
    f.append('pka_res_abs_pd_array', pka_res_abs_pd_array)
with pd.HDFStore(pka_res_rel_pd_file, 'w') as f:
    f.append('pka_res_rel_pd_array', pka_res_rel_pd_array)
with pd.HDFStore(pka_ref_pd_file, 'w') as f:
    f.append('pka_ref_pd_array', pka_ref_pd_array)
with pd.HDFStore(pka_res_ref_pd_file, 'w') as f:
    f.append('pka_res_ref_pd_array', pka_res_ref_pd_array)
