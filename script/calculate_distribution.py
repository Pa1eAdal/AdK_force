# -*- codeing = utf-8 -*-
# Time : 2023/2/1 16:26
# File : calculate_distribution.py
# Software : PyCharm
"""
input: .data file from ts2data
output:.npz file include o2c time and c2o time
"""

import mdtraj as md
import numpy as np
import os
import itertools

project_name = "akesi_force"
special = "confor_standard_test3_rmsd9"

site = [22, 42, 55, 73, 148, 177]
dv_list = [0, 14]
amp_list = [300, 3000]
site_list = [148]
# length_of_dna = [60]
# k_list = [0.01]
# rp_list = np.arange(-4, 2, 2)
# w_list = [10000]
n_seed = 20
rmsd = 0.8
pale = 0

rmsd_cut = 0.9 / 10

top_f = os.path.expanduser(f"~/data/pdb/4ake_5unit0_CA.pdb")
close_f = os.path.expanduser(f"~/data/pdb/1ake_5unit0_CA.pdb")

refer_close = md.load(close_f)

state2label_chem = {}
label_chem2state = {}
state_chem = 'EE TE EM TM DD DE ED DM TD'.split()
state_chem_nod = 'EE TE EM'.split()
state_chem_1d = 'ED DE'.split()
state_chem_1d_long = 'ED DE DM TD'.split()
state_chem_e = "TE EM".split()
state_chem_ee = "EE".split()
state_chem_e_long = "TE EM EE".split()
state_chem_tmdd = "TM DD".split()
state_chem_dd = "DD".split()
state_chem_xd = "ED TD".split()
state_chem_xy = "EE EM TE TM".split()
list_nod = []
list_e = []
list_ee = []
list_1d = []
list_tmdd = []
list_dd = []
list_xd = []
list_xy = []
list_1d_long = []
list_e_long = []

for i in range(len(state_chem)):
    state2label_chem[state_chem[i]] = i
    label_chem2state[i] = state_chem[i]

for i in state_chem_nod:
    list_nod.append(state2label_chem[i])

for i in state_chem_1d:
    list_1d.append(state2label_chem[i])

for i in state_chem_xy:
    list_xy.append(state2label_chem[i])

for i in state_chem_xd:
    list_xd.append(state2label_chem[i])

for i in state_chem_e:
    list_e.append(state2label_chem[i])

for i in state_chem_ee:
    list_ee.append(state2label_chem[i])

for i in state_chem_tmdd:
    list_tmdd.append(state2label_chem[i])

for i in state_chem_dd:
    list_dd.append(state2label_chem[i])

for i in state_chem_1d_long:
    list_1d_long.append(state2label_chem[i])

for i in state_chem_e_long:
    list_e_long.append(state2label_chem[i])

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


def define_state_phy(chi1, chi2):
    """
    :param chi1: LID domain
    :param chi2: NMP domain
    :return: physical state [lid, nmp]
    """

    phy_trj = np.full_like(chi1, 3, dtype=int)

    f0 = (chi1 < 0.5) & (chi1 > -0.5) & (chi2 < 0.5) & (chi2 > -0.5)
    phy_trj[f0] = state2label_phy['MM']

    f1 = (chi1 > 0.5) & (chi2 < 0.5) & (chi2 > -0.5)
    phy_trj[f1] = state2label_phy['CM']

    f2 = (chi1 < 0.5) & (chi1 > -0.5) & (chi2 > 0.5)
    phy_trj[f2] = state2label_phy['MC']

    f3 = (chi1 > 0.5) & (chi2 > 0.5)
    phy_trj[f3] = state2label_phy['CC']

    f4 = (chi1 < -0.5) & (chi2 < -0.5)
    phy_trj[f4] = state2label_phy['OO']

    f5 = (chi1 < -0.5) & (chi2 < 0.5) & (chi2 > -0.5)
    phy_trj[f5] = state2label_phy['OM']

    f6 = (chi1 < 0.5) & (chi1 > -0.5) & (chi2 < -0.5)
    phy_trj[f6] = state2label_phy['MO']

    f7 = (chi1 < -0.5) & (chi2 > 0.5)
    phy_trj[f7] = state2label_phy['OC']

    f8 = (chi1 > 0.5) & (chi2 < -0.5)
    phy_trj[f8] = state2label_phy['CO']

    return phy_trj


def change_time5(state_list, p_state_end, rmsd_list=[], posi=0):
    nframe = 0
    flag = 0
    frame_list = []

    if posi == 0:
        phy_end = state2label_phy[p_state_end]
        if "long" in special:
            list_e2d_start = list_e_long
        if "no_confor" in special:
            if "long" in special:
                list_e2d_start = list_e_long
                list_d2e_start = list_1d_long
            else:
                list_e2d_start = list_e
                list_d2e_start = list_1d
            if "standard" in special:
                list_e2d_start = list_ee
            if "test" in special:
                list_e2d_start = list_e_long
            # list_d2e_end = list_e2d_start
            if phy_end == state2label_phy["CC"]:
                # print(list_e2d_start)
                for i in range(len(state_list)):
                    cs = state_list[i]
                    if (cs in list_e2d_start) & (flag == 0):
                        #  TE EM
                        nframe = 0
                        flag = 1
                    elif (flag == 1) & (cs in list_tmdd):
                        #                        DD
                        nframe += 1
                        frame_list.append(nframe)
                        flag = 0
                    else:
                        nframe += 1
        else:
            if "standard" in special:
                list_e2d_start = list_ee
                list_d2e_start = list_xd
                list_d2e_end = list_xy
            else:
                list_d2e_start = list_1d
                list_d2e_end = list_e
            if phy_end == state2label_phy["CC"]:
                for i in range(len(state_list)):
                    cs = state_list[i]
                    rm = rmsd_list[i]
                    if (cs in list_e2d_start) & (flag == 0):
                        #  EE
                        nframe = 0
                        flag = 1
                    elif (rm < rmsd_cut) & (flag == 1) & (cs in list_tmdd):
                        #       0.8                              TM DD
                        nframe += 1
                        frame_list.append(nframe)
                        flag = 0
                    else:
                        nframe += 1

        if phy_end == state2label_phy["OO"]:
            for i in range(len(state_list)):
                cs = state_list[i]
                if (cs in list_d2e_start) & (flag == 0):
                    #       DE ED
                    nframe = 0
                    flag = 1
                elif (flag == 1) & (cs in list_d2e_end):
                    nframe += 1
                    frame_list.append(nframe)
                    flag = 0
                else:
                    nframe += 1


    return frame_list


def change_time_all(chem_state, rmsd_list):
    o2c = change_time5(chem_state[1:], rmsd_list=rmsd_list, p_state_end="CC")
    c2o = change_time5(chem_state[1:], p_state_end="OO")
    return o2c, c2o


work_dir = os.getcwd() + "/"

out_dir = './distribution/'

if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok=True)

out_head = out_dir.strip("/.")

for dv, co, s in itertools.product(dv_list, amp_list, site_list):
# for co, s, f, w in itertools.product(amp_list, site_list, force_list, w_list):
    # for co, k, rp, s in itertools.product(conc_list, k_list, rp_list, site):
    if dv == 14:
        force_list = [0, 5, 10, 15, 20, 25]
    elif dv == 0:
        force_list = [0, 5, 7, 10, 15, 20]
    c = str(co).zfill(4)
    for f in force_list:
        o2c_list = []
        c2o_list = []
        try:
            w
        except NameError:
            try:
                lod
            except NameError:
                try:
                    f
                except NameError:
                    print("no lod, f and w")
                else:
                    data_dir = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/f{f}/"
                    out_fn = f"{project_name}_{special}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}.npz"
                    if f == 0:
                        data_dir = f"rakesi8/dv{dv}/c{c}/s42/f{f}/"
                        xtc_dir = f"rakesi8/dv{dv}/c{c}/s42/f{f}/"
            else:
                data_dir = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/dp{dp}/l{lod}/io{io}/"
                out_fn = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}.npz"
        else:
            data_dir = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/f{f}/w{w}"
            # dirn = f"{project_name}/dv{dv}/c{c}/s{s}/f{f}/"
            out_fn = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_w{w}.npz"

        out_fn = out_dir + out_fn

        for ns in range(n_seed):
            n = str(ns).zfill(3)

            try:
                w
            except NameError:
                try:
                    f
                except NameError:
                    try:
                        lod
                    except NameError:
                        print("no w, f and lod")
                    else:
                        data_n = data_dir + f"d_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}_n{n}.data"
                else:
                    data_n = data_dir + f"d_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.data"
                    xtc_n = data_dir + f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.xtc"
                    if f == 0:
                        data_n = data_dir + f"d_dv{dv}_c{c}_s42_f{f}_n{n}.data"
                        xtc_n = xtc_dir + f"rakesi8_dv{dv}_c{c}_s42_f{f}_n{n}.xtc"
            else:
                data_n = data_dir + f"d_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_w{w}_n{n}.data"


            # chemical state
            iconf = np.loadtxt(data_n, usecols=0, unpack=True, dtype=int)
            istate1, istate2, istate3, istate4, istate5, istate6 = np.loadtxt(data_n, usecols=np.arange(24, 30),
                                                                              unpack=True,
                                                                              dtype=int)

            istateT = istate3 | istate5
            istateD = istate2 | istate6

            chem_state = define_state_chem(istate1, istateD, istateT, istate4)

            trj = md.load(xtc_n, top=top_f)
            rmsd_list = md.rmsd(trj, refer_close, parallel=False)
            # print(np.max(rmsd_list))

            o2c, c2o = change_time_all(chem_state[1:], rmsd_list)

            c2o_list.append(c2o)
            o2c_list.append(o2c)
        # print(c2o_list)
        c2o_list = [i for ii in c2o_list for i in ii]
        o2c_list = [i for ii in o2c_list for i in ii]
        if dv == 14:
            print(o2c_list)
        print(out_fn)
        np.savez(out_fn, c2o_list=c2o_list, o2c_list=o2c_list)
