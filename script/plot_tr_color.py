# -*- coding: utf-8 -*-
# Author : 
# Time : 2022/11/9 16:59
# File : plot_tr_color.py
# Software : pycharm
"""
Input: .npz file from ake_efficient.py and .txt file from ake_efficient.F90
"""
import itertools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
import os
from time import localtime

from scipy.optimize import curve_fit

import method

control = "force"

task_name = "tr_color" + "_" + control
dv_list = [0, 1, 10, 14, 15, 16]
conc_list = [300]
site = [22, 42, 55, 73, 148, 177]
site1 = [148]
# lod_list = [90]
force_list = [0, 1, 5, 10, 15, 25]
# k_list = [0.01]
# rp_list = np.arange(-4, 3, 2)
rmsd = "0.8"
pale = 0
n_seed = 20
special = f"pic_{localtime().tm_year % 100}_{localtime().tm_mon}_{localtime().tm_mday}"
data_dir = "txtdata/"

if "force" in control:
    project_name = "akesi_force"
elif "dna" in control:
    project_name = "akesi_dna"


def rdata(file, ctrl=0):
    # file is the filename
    # ctrl if how we get average and rmsd
    # rdata return aar, the list as [average, rmsd]
    aar = []
    tr = []
    with open(file, "r") as p:
        for i in p.readlines():
            line = i.strip()
            l = line.split()
            for tr0 in l:
                pass
            tr.append(float(tr0))
        tr = np.asarray(tr)
        tr = tr * 1000
    if ctrl == 0:
        # 20 data
        ave = np.mean(tr)
        aar.append(ave)
        rmsd = np.sqrt(np.mean(tr ** 2 - ave ** 2))
        aar.append(rmsd)
    if ctrl == 1:
        # random.choice 60 data
        p = np.random.choice(tr, 60)
        ave = np.mean(p)
        aar.append(ave)
        rmsd = np.sqrt(np.mean(p ** 2 - ave ** 2))
        aar.append(rmsd)
    if ctrl == 2:
        # random.choice 20 data * 3
        avelist = []
        p1 = np.random.choice(tr, 20)
        p2 = np.random.choice(tr, 20)
        p3 = np.random.choice(tr, 20)
        avelist.append(np.mean(p1))
        avelist.append(np.mean(p2))
        avelist.append(np.mean(p3))
        avelist = np.asarray(avelist)
        ave = np.mean(avelist)
        aar.append(ave)
        rmsd = np.sqrt(np.mean(avelist ** 2 - ave ** 2))
        aar.append(rmsd)

    return aar


# main
plot_setting = method.Plot_setting()
my_dpi = 350
rect = [0.135, 0.14, 0.877, 0.85]
if "both" in special:
    rect = [0.06, 0.11, 0.9, 0.8]
    y_dict = dict()

position_list = ["LID(148)-CORE(22)", "NMP(42)-LID(144)", "NMP(55)-CORE(166)", "NMP(73)-LID(142)", "LID(148)-CORE(177)", "NMP(55)-CORE(177)"]
position_dict = dict(zip(site, position_list))

dvlist = [0, 1, 2, 10, 14, 15, 16]
pcloselist = [0.283, 0.581, 0.775, 0.119, 0.994, 0.019, 0.008]
pclose_dict = dict(zip(dvlist, pcloselist))

dv_list.sort(key=lambda el: pclose_dict[el])

pclose_list = []
for dv in dv_list:
    pclose_list.append(format(pclose_dict[dv], ".2f"))

outdir = "./Figures/tr/pdf"

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

if control == "force":
    # Create output directory if it doesn't exist
    out_data_dir = "./value"
    if not os.path.exists(out_data_dir):
        os.makedirs(out_data_dir, exist_ok=True)
        
    out_data_name = f'{out_data_dir}/turnover_rate_pc.dat'
    with open(out_data_name, 'w') as data_file:
        z = np.zeros((len(force_list), len(dv_list)))
        x = np.arange(0, len(dv_list))
        y = np.arange(0, len(force_list))
        for s in site1:
            for co in conc_list:
                c = str(co).zfill(4)
                # picture setting
                picture_name = f"{project_name}_{task_name}_rm{rmsd}_pi{pale}_s{s}_c{c}"
                pn = special + "_" + picture_name
                fign = outdir.strip("pdf") + pn + ".png"
                pdfn = outdir + "/" + pn + ".pdf"
                fig = plt.figure(figsize=(2.25, 2), dpi=my_dpi)
                ax = fig.add_axes(rect)
                # ax.set_title(f"{position_dict[s]}, [sub] = {c}$\mu$M")
                ax.set_ylabel('Force(pN)')
                ax.set_xlabel('$P_{close}$')
                ax.set_xticks(x, pclose_list)
                ax.set_yticks(y, force_list)
                cmap = plt.colormaps['RdBu']
                # data loading
                for i_f in range(len(force_list)):
                    f = force_list[i_f]
                    for i_dv in range(len(dv_list)):
                        dv = dv_list[i_dv]
                        tr_list = []
                        if ((dv == 0) | (dv == 14) | (dv == 16)) & (f == 0):
                            filename = f"rakesi8{dv}_c{c}_s42_f{f}.txt"
                            fn = data_dir + filename
                            aar = rdata(fn)
                            ave_tr = aar[0]
                            # Write header for this condition
                            data_file.write(f"dv{dv}\tforce{f}\tsite{s}\t[AMP]{co}\n")
                            # Write single value for this case
                            data_file.write(f"n000\t{ave_tr}\n")
                        else:
                            # Write header for this condition
                            data_file.write(f"dv{dv}\tforce{f}\tsite{s}\t[AMP]{co}\n")
                            for seed in range(n_seed):
                                n = f"{seed}".zfill(3)
                                filename = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.npz"
                                fn = data_dir + filename
                                data = np.load(fn)
                                turnover_rate = np.squeeze(data["turnover_rate"]) * 2e5 / 112 / 2
                                tr_list.append(turnover_rate)
                                # Write data for each seed
                                data_file.write(f"n{n}\t{turnover_rate}\n")
                            ave_tr = np.mean(tr_list)
                        z[i_f][i_dv] = ave_tr
                z1 = np.zeros_like(z)
                # z1 = z - z[0]
                z1 = z
                # for i in range(len(z[0])):
                   # z1[i] = z[i] / z[0][i]
                zab = np.absolute(z1)
                zmax = np.max(zab)
                if "0" not in special:
                    cmap=plt.colormaps["viridis_r"]
                    levels = MaxNLocator(nbins=55).tick_values(0, zmax)
                else:
                    levels = MaxNLocator(nbins=55).tick_values(-zmax, zmax)
                norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                cf = ax.pcolormesh(x, y, z1, cmap=cmap, norm=norm, edgecolor='k', linewidth=0.5)
                cb = fig.colorbar(cf, ax=ax, pad=0.03)
                # cb.yaxis.set_major_locator(np.arange(0, 0.6, 0.1))
                # cb.axis.set_minor_locator([])
                cb.set_ticks(np.arange(0, 0.6, 0.1))
                cb.ax.minorticks_off()
                ax.tick_params(length=2.5)

                fig.savefig(fign, dpi=my_dpi)
                print(fign)