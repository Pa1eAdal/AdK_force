# -*- coding: utf-8 -*-
# Author : zzy
# Time : 2023/2/15 21:28
# File : calculate_mfpt_domain.py
# Software : pycharm

import numpy as np
import itertools
import os
import method

"""
Input: .data file got from ts2data
Output: .npz have mfpt "start-state to end-state" time list
"""

# itime,Rg,Etot,chi1,chi2,chi3,Qa,Qb,rmsda,rmsdb,           10 1-10
# ebind1,ebind2,ebind3,ebind4,ebind5,ebind6,ebind7,         7  11-17
# sasa1,sasa2,sasa3,sasa4,sasa5,sasa6,sasa7,                7  18-24
# istate1,istate2,istate3,istate4,istate5,istate6,istate7,  7  25-31
# iake,iconform,ichem,iake2,iake6,iatp,iamp                 7  32-38
# 10+7+7+7+7=38, istate: 24-30

# LID(ATPbd)    1   3   5   7
#               D   T   T*  M**
# NMP(AMPbd)    2   4   6
#               D   M   D*
# ilid = E T D M
# inmp = E M D

project_name = "mfpt_akesi_force_o2c_real"
task_name = "mfpt"
chi_cut = 0.9
dv_list = [0]
concentrate_list = [300, 3000]
position_list = [148]
force_list = [0, 5, 7, 10, 15, 20]

rmsd = 0.8
pale = 0

n_seed = 100
frame = 4e5

if "c2o" in project_name:
    start_state = "DD"
    p_s_start = "CC"
    end_state = ""  # any
    p_s_end = "OO"
if "o2c" in project_name:
    start_state = "EE"
    p_s_start = "OO"
    end_state = ""  # any
    p_s_end = "CC"

outdir = './domain_mfpt/'

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

outhead = outdir.strip("./")

tau = 112 / 2e5  # ms

state2label_chem = {}
label_chem2state = {}
state_chem = 'EE TE EM TM DD DE ED DM TD'.split()
state_chem_nod = 'EE TE EM'.split()
state_chem_1d = 'TD ED'.split()
list_nod = []
list_1d = []

for i in range(len(state_chem)):
    state2label_chem[state_chem[i]] = i
    label_chem2state[i] = state_chem[i]

for i in state_chem_nod:
    list_nod.append(state2label_chem[i])

for i in state_chem_1d:
    list_1d.append(state2label_chem[i])

state2label_phy = {}
label_phy2state = {}
state_phy = "MM CM MC CC OO OM MO OC CO".split()

for i in range(len(state_phy)):
    state2label_phy[state_phy[i]] = i
    label_phy2state[i] = state_phy[i]


def define_state_chem(istate1, istateD, istateT, istate4):
    cg_trj = np.full_like(istate1, 0, dtype=int)
    cg_trj_str = []
    '''
    if x is a np.ndarray as [1, 2, 1, 3],(x = 1) return [T, F, T, F]
    so cg_trj = [0, 0, ...0]first
    then it will be [2, 1, 1, 2, 3, 2, 1, ..., 2]-like
    '''
    f0 = (istate1 == 0) & (istateD == 0) & (istateT == 0) & (istate4 == 0)
    cg_trj[f0] = state2label_chem['EE']

    f1 = (istate1 == 0) & (istateD == 0) & (istateT > 0) & (istate4 == 0)
    cg_trj[f1] = state2label_chem['TE']

    f2 = (istate1 == 0) & (istateD == 0) & (istateT == 0) & (istate4 > 0)
    cg_trj[f2] = state2label_chem['EM']

    f3 = (istate1 == 0) & (istateD == 0) & (istateT > 0) & (istate4 > 0)
    cg_trj[f3] = state2label_chem['TM']

    f4 = (istate1 > 0) & (istateD > 0) & (istateT == 0) & (istate4 == 0)
    cg_trj[f4] = state2label_chem['DD']

    f5 = (istate1 > 0) & (istateD == 0) & (istateT == 0) & (istate4 == 0)
    cg_trj[f5] = state2label_chem['DE']

    f6 = (istate1 == 0) & (istateD > 0) & (istateT == 0) & (istate4 == 0)
    cg_trj[f6] = state2label_chem['ED']

    f7 = (istate1 > 0) & (istateD == 0) & (istateT == 0) & (istate4 > 0)
    cg_trj[f7] = state2label_chem['DM']

    f8 = (istate1 == 0) & (istateD > 0) & (istateT > 0) & (istate4 == 0)
    cg_trj[f8] = state2label_chem['TD']

    return cg_trj


def define_state_phy(chi1, chi2, chi_cut):
    """
    :param chi1: LID domain
    :param chi2: NMP domain
    :return: physical state [lid, nmp]
    """

    phy_trj = np.full_like(chi1, 3, dtype=int)

    f0 = (chi1 < chi_cut) & (chi1 > -chi_cut) & (chi2 < chi_cut) & (chi2 > -chi_cut)
    phy_trj[f0] = state2label_phy['MM']

    f1 = (chi1 > chi_cut) & (chi2 < chi_cut) & (chi2 > -chi_cut)
    phy_trj[f1] = state2label_phy['CM']

    f2 = (chi1 < chi_cut) & (chi1 > -chi_cut) & (chi2 > chi_cut)
    phy_trj[f2] = state2label_phy['MC']

    f3 = (chi1 > chi_cut) & (chi2 > chi_cut)
    phy_trj[f3] = state2label_phy['CC']

    f4 = (chi1 < -chi_cut) & (chi2 < -chi_cut)
    phy_trj[f4] = state2label_phy['OO']

    f5 = (chi1 < -chi_cut) & (chi2 < chi_cut) & (chi2 > -chi_cut)
    phy_trj[f5] = state2label_phy['OM']

    f6 = (chi1 < chi_cut) & (chi1 > -chi_cut) & (chi2 < -chi_cut)
    phy_trj[f6] = state2label_phy['MO']

    f7 = (chi1 < -chi_cut) & (chi2 > chi_cut)
    phy_trj[f7] = state2label_phy['OC']

    f8 = (chi1 > chi_cut) & (chi2 < -chi_cut)
    phy_trj[f8] = state2label_phy['CO']

    return phy_trj


def change_time_domain_chem(domain, chem_state_list, chem_state_start, phy_state_list, phy_state_start, phy_state_end):
    if domain == "nmp":
        posi = 1
    elif domain == "lid":
        posi = 0

    nframe = 0
    flag = 0

    phy_end = state2label_phy[phy_state_end]
    if phy_end == state2label_phy["CC"]:
        for i in range(len(chem_state_list)):
            cs = label_chem2state[chem_state_list[i]][posi]
            ps = label_phy2state[phy_state_list[i]][posi]
            if (flag == 0) & (ps == "O"):
                nframe = 0
                flag = 1
            elif (ps == "C") & (flag == 1) & (cs == "TM"[posi]):
                nframe += 1
                break
            else:
                nframe += 1
    if phy_end == state2label_phy["OO"]:
        for i in range(len(chem_state_list)):
            cs = label_chem2state[chem_state_list[i]][posi]
            ps = label_phy2state[phy_state_list[i]][posi]
            if (ps == "C") & (flag == 0):
                nframe = 0
                flag = 1
            elif (flag == 1) & (cs != "D") & (ps == "O"):
                nframe += 1
                break
            else:
                nframe += 1

    return nframe


def change_time_domain_phy(domain, phy_state_list, phy_state_start, phy_state_end):
    if domain == "nmp":
        posi = 1
    elif domain == "lid":
        posi = 0

    nframe = 0
    flag = 0

    phy_end = state2label_phy[phy_state_end]
    if phy_end == state2label_phy["CC"]:
        for i in range(len(phy_state_list)):
            ps = label_phy2state[phy_state_list[i]][posi]
            if (flag == 0) & (ps == "O"):
                nframe = 0
                flag = 1
            elif (ps == "C") & (flag == 1):
                nframe += 1
                break
            else:
                nframe += 1
    elif phy_end == state2label_phy["OO"]:
        for i in range(len(phy_state_list)):
            ps = label_phy2state[phy_state_list[i]][posi]
            if (ps == "C") & (flag == 0):
                nframe = 0
                flag = 1
            elif (flag == 1) & (ps == "O"):
                nframe += 1
                break
            else:
                nframe += 1

    return nframe


# main
for dv, co, s, f in itertools.product(dv_list, concentrate_list, position_list, force_list):
    # for dv, co, s, lod in itertools.product(dv_list, concentrate_list, position_list, lod_list):
    c = f"{co}".zfill(4)
    try:
        f
    except NameError:
        try:
            lod
        except NameError:
            print("no lod and f")
        else:
            outfn = f"{project_name}_cut{chi_cut}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}"
    else:
        outfn = f"{project_name}_cut{chi_cut}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}"

    print(outfn + '\n')

    nmp_tot_time_chem = []
    lid_tot_time_chem = []

    nmp_tot_time_phy = []
    lid_tot_time_phy = []

    for seed in range(n_seed):
        n = f"{seed}".zfill(3)
        try:
            f
        except NameError:
            try:
                lod
            except NameError:
                print("no lod and f")
            else:
                data_n = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/dp{dp}/l{lod}/io{io}/" \
                         f"d_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}_n{n}.data"
        else:
            if (f != 0) & ("long" in project_name):
                data_n = f"mfpt_without_cata/{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/f{f}/" \
                         f"d_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.data"
            else:
                data_n = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/f{f}/" \
                         f"d_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.data"

        # physical state
        chi1, chi2 = np.loadtxt(data_n, usecols=[3, 4], unpack=True, dtype=float)  # chi1-lid, chi2-nmp
        phy_state_list = define_state_phy(chi1, chi2, chi_cut)

        # chemical state
        iconf = np.loadtxt(data_n, usecols=0, unpack=True, dtype=int)
        istate1, istate2, istate3, istate4, istate5, istate6 = np.loadtxt(data_n, usecols=np.arange(24, 30),
                                                                          unpack=True,
                                                                          dtype=int)

        istateT = istate3 | istate5
        istateD = istate2 | istate6

        cg_trj_chem = define_state_chem(istate1, istateD, istateT, istate4)

        nmp_frame_chem = change_time_domain_chem("nmp", cg_trj_chem, "DD", phy_state_list, phy_state_start=p_s_start, phy_state_end=p_s_end)
        lid_frame_chem = change_time_domain_chem("lid", cg_trj_chem, "DD", phy_state_list, phy_state_start=p_s_start, phy_state_end=p_s_end)

        nmp_frame_phy = change_time_domain_phy("nmp", phy_state_list, phy_state_start=p_s_start, phy_state_end=p_s_end)
        lid_frame_phy = change_time_domain_phy("lid", phy_state_list, phy_state_start=p_s_start, phy_state_end=p_s_end)

        nmp_tot_time_chem.append(nmp_frame_chem)
        lid_tot_time_chem.append(lid_frame_chem)

        nmp_tot_time_phy.append(nmp_frame_phy)
        lid_tot_time_phy.append(lid_frame_phy)

    outname = outdir + outhead + "_" + outfn + "_" + f"{p_s_start}_{p_s_end}"
    np.savez(outname, nmp_time_chem=nmp_tot_time_chem, lid_time_chem=lid_tot_time_chem, nmp_time_phy=nmp_tot_time_phy,
             lid_time_phy=lid_tot_time_phy)
    print(outname)
