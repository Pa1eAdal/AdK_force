# -*- codeing = utf-8 -*-
# Time : 2022/8/23 23:44
# File : ake_efficient.py
# Software : PyCharm
"""
Input: .data file from ts2data.F90
"""
import numpy as np
import itertools
import os
# from method import complement_data_name

control = "chem"
project_name = "akesi"
dv_list = [0, 14, 16]
amp_list = [10, 30, 100, 300, 1000, 3000]
atp_list = [1000]
adp_list = [0]
cut = "rmsd"
site_list = [148]
# lod_list= [40, 60, 70]
# io_list = [0.05]
# dp_list = [0]
# force_list = [5, 10, 15, 20, 25]
force_list = [0]
# w_list = [10000]
# k_list = [0.015]
# rp_list = np.asarray([0, 3, 12])
# ba_list = [[1, 1]]
# cf_list = [4]
# cb_list = [5]
# kof_list = [0.1]

n_seed = 20
frame = 2e5
rmsd = "0.8"
pale = 0

if cut == "chi":
    chi = 0.5
    cut_val = chi
    pale = 1
elif cut == "rmsd":
    cut_val = 0.8
    pale = 0

if pale == 0:
    p1 = 0.09
    p7 = 0.36
elif pale == 1:
    p1 = 0.11
    p7 = 0.41

outdir = "./txtdata/"

# for i in range(0, 9)
# state_list = "EE TE EM TM DD DE ED DM TD".split()
# ichem      =  0  1  2  3  4  5  6  7  8
# ia2        =  0  1  1  23 45 6  6  7  7
# iake2      =  0  1  1  2  3  4  4  5  5
# this is ia2, in other word, state_list[ia2] is the corresponding chemical state

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

dvn_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16]
p_close = [0.28, 0.52, 0.77, 0.15, 0.05, 0.02, 9e-3, 3e-3, 0.89, 0.94, 0.12, 0.10, 0.09, 0.986, 0.994, 0.008]
p_dict = dict(zip(dvn_list, p_close))

if control == "amend":  # for amendment
    for dv, amp, atp, adp, f in itertools.product(dv_list, amp_list, atp_list, adp_list, force_list):
        # c = format(conc, "04d")
        t = f"{atp:04d}"
        m = f"{amp:04d}"
        d = f"{adp:04d}"
        # work_dir = f"cefemol-ake-7/workdir/dv{dv}_g{p7}/t{t}_m{m}/"
        # work_dir = f"{project_name}/dv{dv}/chi{chi}/p1{p1}/p7{p7}/t{t}/m{m}/d{d}/"
        work_dir = f"{project_name}/dv{dv}/c{m}/s42/f{0}/"
        fn = work_dir + "data0filelist"
        # calculate turnover rate
        with open(fn, "r") as dlist:
            for datafile in dlist.readlines():
                catalyze3_old = -1
                catalyze3 = -1
                turnover3 = 0
                count = 0
                nframe = []
                data_n = work_dir + datafile.strip("\n")
                out_fn = outdir + f"akesi_force_dv{dv}_rm{rmsd}_pi{pale}_c{m}_s148_f{f}_n" + datafile.strip("\n.data")[-3:]
                # out_fn = out_fn.strip(".data")
                with open(data_n, 'r') as dfn:
                    for data_s in dfn.readlines():
                        data = data_s.split()
                        frame = int(data[0])
                        ia2 = int(data[34])  # data[34] not only data[-2]
                        if ia2 == 0:
                            iake2 = 0
                        if ia2 == 1:
                            iake2 = 1
                        if ia2 == 2:
                            iake2 = 2
                        if ia2 == 3:
                            iake2 = 2
                        if ia2 == 4:
                            iake2 = 3
                        if ia2 == 5:
                            iake2 = 3
                        if ia2 == 6:
                            iake2 = 4
                        if ia2 == 7:
                            iake2 = 5
                        if (iake2 == 0) or (iake2 == 1) or (iake2 == 2) or (iake2 == 3):
                            catalyze3 = -1
                        if (iake2 == 4) or (iake2 == 5):
                            catalyze3 = 1
                        if catalyze3 * catalyze3_old == -1:
                            turnover3 = turnover3 + 1
                            nframe.append(frame)
                            # p_eff.write(f"{frame}" + "\n")
                            catalyze3_old = catalyze3
                            # print(frame)
                        count += 1
                turnover3 = turnover3 / count
                np.savez(out_fn, nframe=nframe, turnover_rate=turnover3)
                print(out_fn)
else:
    for dv, amp, atp, adp, f, s in itertools.product(dv_list, amp_list, atp_list, adp_list, force_list, site_list):
        # c = format(conc, "04d")
        t = f"{atp:04d}"
        m = f"{amp:04d}"
        c = m
        try:
            adp
        except NameError:
            try:
                f
            except NameError:
                try:
                    lod
                except NameError:
                    print("no lod and f")
                else:
                    outfn = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}"
            else:
                work_dir = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/f{f}/"
        else:
            d = f"{adp:04d}"
            # work_dir = f"cefemol-ake-7/workdir/dv{dv}_g{p7}/t{t}_m{m}/"
            # work_dir = f"{project_name}/dv{dv}/chi{chi}/p1{p1}/p7{p7}/t{t}/m{m}/d{d}/"
        work_dir = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/f{f}/"
        fn = work_dir + "data0filelist"
        # calculate turnover rate
        with open(fn, "r") as dlist:
            for datafile in dlist.readlines():
                catalyze3_old = -1
                catalyze3 = -1
                turnover3 = 0
                count = 0
                nframe = []
                data_n = work_dir + datafile.strip("\n")
                out_fn = outdir + f"{project_name}_d" + datafile.strip("\nd_").strip(".data")
                # out_fn = out_fn.strip(".data")
                with open(data_n, 'r') as dfn:
                    for data_s in dfn.readlines():
                        data = data_s.split()
                        frame = int(data[0])
                        ia2 = int(data[34])  # data[34] not only data[-2]
                        if ia2 == 0:
                            iake2 = 0
                        if ia2 == 1:
                            iake2 = 1
                        if ia2 == 2:
                            iake2 = 2
                        if ia2 == 3:
                            iake2 = 2
                        if ia2 == 4:
                            iake2 = 3
                        if ia2 == 5:
                            iake2 = 3
                        if ia2 == 6:
                            iake2 = 4
                        if ia2 == 7:
                            iake2 = 5
                        if (iake2 == 0) or (iake2 == 1) or (iake2 == 2) or (iake2 == 3):
                            catalyze3 = -1
                        if (iake2 == 4) or (iake2 == 5):
                            catalyze3 = 1
                        if catalyze3 * catalyze3_old == -1:
                            turnover3 = turnover3 + 1
                            nframe.append(frame)
                            # p_eff.write(f"{frame}" + "\n")
                            catalyze3_old = catalyze3
                            # print(frame)
                        count += 1
                turnover3 = turnover3 / count
                np.savez(out_fn, nframe=nframe, turnover_rate=turnover3)
                print(out_fn)