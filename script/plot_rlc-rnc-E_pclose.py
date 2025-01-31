#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 20:43:40 2021

@author: yyzhang

Modified on Fri Aug  5 23:50 2022

@developer: zzy
"""
"""
Input: .npz file from calc_rlc-rnc_pclose.py
"""
import numpy as np
import itertools
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import ticker, cm
import os
from time import localtime
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# import matplotlib as mpl
# mpl.rcParams['axes.unicode_minus'] = False

# plt.rcParams.items()
plt.rcParams.update({'text.usetex': False})
rc('font', **{'style': 'normal', 'size': 8, 'family': 'sans-serif', 'sans-serif': ['DejaVu Sans']})
rc('axes', **{'titlesize': 8})
rc('mathtext', **{'default': 'regular'})
rc(('xtick', 'ytick'), **{'labelsize': 8, 'direction': 'in', 'major.pad': 1})
rc(('axes'), **{'labelsize': 10, 'labelpad': 1, 'linewidth': 1})
rc(('legend'), **{'fontsize': 8, 'frameon': False})

control = "paper2"
special = f"{localtime().tm_year % 100}_{localtime().tm_mon}_{localtime().tm_mday}"

task_name = "rlc-rnc"
project_name = "akesi_dna"
dv_list = [0]
conc_list = [0]
site_list = [148]
force_list = [0]
lod_list = [70]
dp = 0
io = 0.05
# k_list = [0.005, 0.01]
# rp_list = np.arange(7) * 2
rmsd = 0.8
pale = 0

my_dpi = 350
rect = [0.15, 0.15, 0.75 * 2.5 / 2.25, 0.75]

position_list = ["LID(148)-CORE(22)", "NMP(42)-LID(144)", "NMP(55)-CORE(166)", "NMP(73)-LID(142)"]
position_dict = dict(zip(site_list, position_list))

dvnlist = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16]
pclose = [0.28, 0.52, 0.77, 0.15, 0.05, 0.02, 9e-3, 3e-3, 0.89, 0.94, 0.12, 0.10, 0.09, 0.986, 0.994, 0.008]
pdict = dict(zip(dvnlist, pclose))

outdir = './Figures/pdf'

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

if control == "paper2":
    if "dna" in project_name:
        for f in force_list:
            fig, ax1 = plt.subplots(1, 1, figsize=(2.25, 2.25), dpi=my_dpi)
            plt.subplots_adjust(left=0.15, bottom=0.13, right=0.971, top=0.98, wspace=0)
            for dv in dv_list:
                p = format(pdict[dv], '.2f')
                print(f"pclose={p}")
                for s, co, lod in itertools.product(site_list, conc_list, lod_list):
                    c = f"{co}".zfill(4)
                    prefix = f"{task_name}_{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}"
                    fn = prefix + ".npz"
                    fign = outdir.strip(
                        "pdf") + f"{task_name}_{special}_{project_name}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}_sa.png"
                    data = np.load(fn)
                    x, y = 10 * data['rlc'], 10 * data['rnc']
                    print(x.shape, y.shape)

                    H, xedges, yedges = np.histogram2d(x, y, bins=100)
                    xmid = 0.5 * (xedges[1:] + xedges[:-1])
                    ymid = 0.5 * (yedges[1:] + yedges[:-1])

                    X, Y = np.meshgrid(xmid, ymid)
                    print(H.max())

                    z = np.ma.masked_where(H <= 0, -np.log(H))
                    print(fign)
                    print(z.min())
                    z0 = z - z.min()

                    vmin = 0
                    vmax = z0.max()
                    levels = np.arange(vmin, vmax + 0.5, step=1)

                    plt.rcParams.update({'text.usetex': True})
                    ax1.set_xlabel('$R_{LID \u2010 CORE}(\AA)$')
                    ax1.set_ylabel('$R_{NMP \u2010 CORE}(\AA)$')
                    ax1.set_xlim([19, 37])
                    ax1.set_ylim([17, 25])
                    ax1.set_xticks([20, 25, 30, 35])
                    ax1.set_yticks([18, 20, 22, 24])
                    #cs = ax1.contour(X, Y, z0.T)
                    #cntr = ax1.contourf(X, Y, z0.T)
                    cs = ax1.contour(X, Y, z0.T, vmin=vmin, vmax=vmax, levels=levels, colors='k', linestyles='-',
                                     linewidths=0.1)
                    cntr = ax1.contourf(X, Y, z0.T, vmin=vmin, vmax=vmax, levels=levels, cmap=cm.coolwarm)
            # axins = inset_axes(
            #     ax1,
            #     width="5%",  # width: 5% of parent_bbox width
            #     height="30%",  # height: 50%
            #     loc="lower left",
            #     bbox_to_anchor=(0.9, 0., 1, 1),
            #     bbox_transform=ax1.transAxes,
            #     borderpad=0,
            # )
            # fig.colorbar(ax1, cax=axins)

            cb_ax = fig.add_axes([0.865, 0.16, 0.015, 0.4])
            cb_ax.spines['top'].set_linewidth(10)
            cb_ax.spines['top'].set_color('red')
            cbar = fig.colorbar(cntr, cax=cb_ax, ticks=np.arange(0, vmax+1, 3))
            # cbar.set_label('$k_B T$', loc='top')
            cbar.ax.tick_params(axis='y', direction='out', labelsize=8, length=2)

            fig.text(0.875, 0.595, '$\mathrm{k_{B}T}$', fontsize=8)
            fig.savefig(fign, dpi=my_dpi)
            plt.close("all")
    else:
        for f in force_list:
            fig, ax1 = plt.subplots(1, 1, figsize=(2.25, 2.25), dpi=my_dpi)
            plt.subplots_adjust(left=0.15, bottom=0.13, right=0.971, top=0.98, wspace=0)
            for dv in dv_list:
                p = format(pdict[dv], '.2f')
                print(f"pclose={p}")
                for s, co in itertools.product(site_list,conc_list):
                    c = f"{co}".zfill(4)
                    prefix = f"{task_name}_{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}"
                    fn = prefix + ".npz"
                    fign = outdir.strip("pdf") + f"{task_name}_{special}_{project_name}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_sa.png"
                    data = np.load(fn)
                    x, y = 10 * data['rlc'], 10 * data['rnc']
                    print(x.shape, y.shape)

                    H, xedges, yedges = np.histogram2d(x, y, bins=100)
                    xmid = 0.5 * (xedges[1:] + xedges[:-1])
                    ymid = 0.5 * (yedges[1:] + yedges[:-1])

                    X, Y = np.meshgrid(xmid, ymid)

                    z = np.ma.masked_where(H <= 0, -np.log(H))
                    print(fign)
                    print(z.min())
                    z0 = z - z.min()

                    vmin = 0
                    vmax = 12
                    levels = np.arange(vmin, vmax + 0.5, step=0.8)

                    plt.rcParams.update({'text.usetex': True})
                    ax1.set_xlabel('$R_{LID \u2010 CORE}(\AA)$')
                    ax1.set_ylabel('$R_{NMP \u2010 CORE}(\AA)$')
                    ax1.set_xlim([19, 37])
                    ax1.set_ylim([17, 25])
                    ax1.set_xticks([20, 25, 30, 35])
                    ax1.set_yticks([18, 20, 22, 24])
                    cs = ax1.contour(X, Y, z0.T, vmin=vmin, vmax=vmax, levels=levels, colors='k', linestyles='-',
                                     linewidths=0.1)
                    cntr = ax1.contourf(X, Y, z0.T, vmin=vmin, vmax=vmax, levels=levels, cmap=cm.coolwarm)
            # axins = inset_axes(
            #     ax1,
            #     width="5%",  # width: 5% of parent_bbox width
            #     height="30%",  # height: 50%
            #     loc="lower left",
            #     bbox_to_anchor=(0.9, 0., 1, 1),
            #     bbox_transform=ax1.transAxes,
            #     borderpad=0,
            # )
            # fig.colorbar(ax1, cax=axins)

            cb_ax = fig.add_axes([0.865, 0.16, 0.015, 0.4])
            cb_ax.spines['top'].set_linewidth(10)
            cb_ax.spines['top'].set_color('red')
            cbar = fig.colorbar(cntr, cax=cb_ax, ticks=np.arange(0, 11, 2))
            # cbar.set_label('$k_B T$', loc='top')
            cbar.ax.tick_params(axis='y', direction='out', labelsize=8, length=2)

            fig.text(0.875, 0.595, '$\mathrm{k_{B}T}$', fontsize=8)
            fig.savefig(fign, dpi=my_dpi)
            plt.close("all")
