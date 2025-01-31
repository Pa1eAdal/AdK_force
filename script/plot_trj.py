# -*- codeing = utf-8 -*-
# Time : 2023/2/8 18:13
# File : plot_trj.py
# Software : PyCharm
"""
Input: .npz file from calculate_domain_distance.py
"""
import numpy as np
from matplotlib import pyplot as plt
import itertools
import enzyme_setting
import method
import os

plot_setting = method.Plot_setting()
adk = enzyme_setting.Enzyme_adk(cut="chi")
force_setting = enzyme_setting.Force_setting()

task_name = "trj"
project_name = "akesi_force"

adk.dv_list = [0]
adk.amp_list = [300]
adk.adp_list = [0]
force_setting.site_list = [148]
if project_name == "akesi_force":
    site = [148]
    force_setting.force_list = [10]
elif project_name == "rakesi8":
    site = [42]
    force_setting.force_list = [0]
n_seed = 1
# length_of_dna = [60]
# klist = [0.01]
# rplist = np.arange(0, 1)
special = "10ms"

my_dpi = plot_setting.dpi
rect = [0.06, 0.13, 0.9, 0.8]

data_dir = "domain_distance"
tr_dir = "txtdata"
outdir = "./Figures/trj/pdf"
method.mk_outdir(outdir)

if "dna" in project_name:
    io_list = [0.05]
    dp_list = [0]
    lod_list = [70]
    site = [148]
    for dv, amp, s, io, dp, lod in itertools.product(dv_list, amp_list, site, io_list, dp_list, lod_list):
        c = f"{amp}".zfill(4)
        if lod == 70:
            lim = [0, 112]
            # lim = [20, 30]
        lim_list = np.arange(lim[0], lim[1] + 1, 2)
        fdir = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/dp{dp}/l{lod}/io{io}"
        special = f"{lim[-1] - lim[0]}ms"

        for seed in range(n_seed):
            n = f"{seed:03d}"
            domdis_n = f"./domain_distance/domain_distance_{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}_n{n}.npz"
            domdis = np.load(domdis_n)
            rlc = domdis["rlc"] * 10
            rnc = domdis["rnc"] * 10
            # z = np.loadtxt(filename)  # file['frame'],file['pframe']
            # chi1, chi2 = [], []
            # for line in z:
            #     chi1.append(line[3])
            #     chi2.append(line[4])

            frame = np.arange(0, 20001)
            frame = frame * 112 / 2e4

            fig, (ax3, ax4) = plt.subplots(2, 1, figsize=(2.25, 2), dpi=my_dpi, sharex=True)
            plt.subplots_adjust(left=0.15, bottom=0.13, right=0.971, top=0.98, hspace=0.05)
            # fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(4.5, 4), dpi=my_dpi)

            fign = outdir.strip(
                "pdf") + f"{project_name}_{special}_{task_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}_n{n}.png"

            new_rlc = method.generate_new_array(rlc, 10)
            new_rnc = method.generate_new_array(rnc, 10)

            # ax1.plot(frame, chi1[1:], zorder=10)
            # ax1.set_ylabel("$\chi_{LID}$")
            # # ax1.set_xlabel("time(ms)")
            # # ax1.set_xlim(lim)
            # ax1.set_ylim(-1.1, 1.1)
            # ax1.set_yticks([-1, 0, 1])
            # ax1.set_xticks([])

            # ax2.plot(frame, chi2[1:], zorder=10)
            # ax2.set_ylabel("$\chi_{NMP}$")
            # # ax2.set_xlabel("time(ms)")
            # # ax2.set_xlim(lim)
            # ax2.set_ylim(-1.1, 1.1)
            # ax2.set_yticks([-1, 0, 1])
            # ax2.set_xticks([])

            ax3.plot(frame, new_rlc, lw=0.5, zorder=10, color='tab:blue')
            ax3.set_ylabel("$R_{LID-CORE}(\AA)$")
            # ax3.set_xlabel("time(ms)")
            ax3.set_ylim(19, 34)
            ax3.set_yticks([22, 29])
            ax3.set_xlim(lim)
            ax3.set_xticks([])
            ax3.grid(lw=0.5, axis='y')

            ax4.plot(frame, new_rnc, lw=0.5, zorder=10, color='tab:blue')
            ax4.set_ylabel("$R_{NMP-CORE}(\AA)$")
            ax4.set_xlabel("Time(ms)")
            ax4.set_ylim(17.5, 23.5)
            ax4.set_yticks([19, 22])
            ax4.set_xlim(lim)
            ax4.set_xticks(lim_list, np.arange(0, lim_list[-1] - lim_list[0] + 1, 2))
            ax4.grid(lw=0.5, axis='y')

            fig.savefig(fign, dpi=my_dpi)
            print(fign)
elif "force" in project_name:
    for dv, amp, s, f in itertools.product(adk.dv_list, adk.amp_list, force_setting.site_list, force_setting.force_list):
        c = f"{amp}".zfill(4)
        if f == 0:
            lim = [20, 30]
        elif f == 10:
            lim = [32, 42]
        lim_list = np.arange(lim[0], lim[1] + 1, 2)
        if project_name == "rakesi8":
            fdir = f"{project_name}/dv{dv}/c{c}/s{s}/f{f}"
        else:
            fdir = f"{project_name}/dv{dv}/rm{rmsd}/pi{pale}/c{c}/s{s}/f{f}"

        for seed in range(n_seed):
            n = f"{seed}".zfill(3)
            if project_name == "rakesi8":
                filename = f"{fdir}/d_dv{dv}_c{c}_s{s}_f{f}_n{n}.data"
            else:
                filename = f"{fdir}/d_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.data"
            domdis_n = f"./domain_distance/domain_distance_{project_name}_dv{dv}_c{c}_s{s}_f{f}_n{n}.npz"
            domdis = np.load(domdis_n)
            rlc = domdis["rlc"] * 10
            rnc = domdis["rnc"] * 10
            # z = np.loadtxt(filename)  # file['frame'],file['pframe']
            # chi1, chi2 = [], []
            # for line in z:
            #     chi1.append(line[3])
            #     chi2.append(line[4])

            frame = np.arange(0, 20001)
            frame = frame * 112 / 2e4

            fig, (ax3, ax4) = plt.subplots(2, 1, figsize=(2.25, 2), dpi=my_dpi, sharex=True)
            plt.subplots_adjust(left=0.15, bottom=0.13, right=0.975, top=0.98, hspace=0.05)
            # fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(4.5, 4), dpi=my_dpi)

            fign = outdir.strip("pdf") + f"{project_name}_{special}_{task_name}_dv{dv}_rm{rmsd}_pi{pale}_s{s}_f{f}_n{n}.png"

            new_rlc = method.generate_new_array(rlc, 10)
            new_rnc = method.generate_new_array(rnc, 10)

            # ax1.plot(frame, chi1[1:], zorder=10)
            # ax1.set_ylabel("$\chi_{LID}$")
            # # ax1.set_xlabel("time(ms)")
            # # ax1.set_xlim(lim)
            # ax1.set_ylim(-1.1, 1.1)
            # ax1.set_yticks([-1, 0, 1])
            # ax1.set_xticks([])

            # ax2.plot(frame, chi2[1:], zorder=10)
            # ax2.set_ylabel("$\chi_{NMP}$")
            # # ax2.set_xlabel("time(ms)")
            # # ax2.set_xlim(lim)
            # ax2.set_ylim(-1.1, 1.1)
            # ax2.set_yticks([-1, 0, 1])
            # ax2.set_xticks([])

            ax3.plot(frame, new_rlc, lw=0.5, zorder=10, color='tab:blue')
            ax3.set_ylabel("$R_{LID-CORE}(\AA)$")
            # ax3.set_xlabel("time(ms)")
            ax3.set_ylim(19, 34)
            ax3.set_yticks([22, 29])
            ax3.set_xlim(lim)
            ax3.set_xticks([])
            ax3.grid(lw=0.5, axis='y')

            ax4.plot(frame, new_rnc, lw=0.5, zorder=10, color='tab:blue')
            ax4.set_ylabel("$R_{NMP-CORE}(\AA)$")
            ax4.set_xlabel("Time(ms)")
            ax4.set_ylim(17.5, 23.5)
            ax4.set_yticks([19, 22])
            ax4.set_xlim(lim)
            ax4.set_xticks(lim_list, np.arange(0, lim_list[-1] - lim_list[0] + 1, 2))
            ax4.grid(lw=0.5, axis='y')

            fig.savefig(fign, dpi=my_dpi)
            print(fign)
