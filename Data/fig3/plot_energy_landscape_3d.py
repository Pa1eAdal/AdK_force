#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 20:43:40 2021

@author: yyzhang

Modified on Fri Aug  5 23:50 2022

@developer: zhangzy
"""

import numpy as np
import itertools
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import ticker, cm
import os
from time import localtime
import plotly.graph_objects as go
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D
# import matplotlib as mpl
# mpl.rcParams['axes.unicode_minus'] = False

# plt.rcParams.items()
# plt.rcParams.update({'text.usetex': False})
rc('font', **{'style': 'normal', 'size': 8, 'family': 'sans-serif', 'sans-serif': ['DejaVu Sans']})
rc('axes', **{'titlesize': 8})
rc('mathtext', **{'default': 'regular'})
rc(('xtick', 'ytick'), **{'labelsize': 8, 'direction': 'in', 'major.pad': 1})
rc(('axes'), **{'labelsize': 10, 'labelpad': 1, 'linewidth': 1})
rc(('legend'), **{'fontsize': 8, 'frameon': False})

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.fontset'] = 'dejavusans'

special = f"{localtime().tm_year % 100}_{localtime().tm_mon}_{localtime().tm_mday}"

task_name = "rlc-rnc"
project_name = "akesi_force"
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
            fign = outdir.strip("pdf") + f"{task_name}_dv{dv}_{special}_{project_name}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_sa.png"
            data = np.load(fn)
            x, y = 10 * data['rlc'], 10 * data['rnc']
            print(x.shape, y.shape)
            H, xedges, yedges = np.histogram2d(x, y, bins=40)
            xmid = 0.5 * (xedges[1:] + xedges[:-1])
            ymid = 0.5 * (yedges[1:] + yedges[:-1])
            X, Y = np.meshgrid(xmid, ymid)
            z = np.ma.masked_where(H <= 0, -np.log(H))
            print(fign)
            print(z.min())
            z0 = z - z.min()
            # Find minima in lower left and upper right quadrants
            x_mid_idx = len(xmid) // 4
            y_mid_idx = len(ymid) // 2
            
            # Lower left quadrant
            z0_ll = z0[:y_mid_idx, :x_mid_idx]
            min_ll_idx = np.unravel_index(np.argmin(z0_ll), z0_ll.shape)
            min_ll_val = z0_ll[min_ll_idx]
            min_ll_x = xmid[min_ll_idx[1]]
            min_ll_y = ymid[min_ll_idx[0]]
            print(f"Lower left minimum: z0={min_ll_val:.2f} at x={min_ll_x:.2f}, y={min_ll_y:.2f}")
            
            # Upper right quadrant  
            z0_ur = z0[y_mid_idx:, x_mid_idx:]
            min_ur_idx = np.unravel_index(np.argmin(z0_ur), z0_ur.shape)
            min_ur_val = z0_ur[min_ur_idx]
            min_ur_x = xmid[min_ur_idx[1] + x_mid_idx]
            min_ur_y = ymid[min_ur_idx[0] + y_mid_idx]
            print(f"Upper right minimum: z0={min_ur_val:.2f} at x={min_ur_x:.2f}, y={min_ur_y:.2f}")
            #if f == 0:
            #    z0 += 1.21
            #if f == 10:
            #    z0 += 0.74
            space = 1.5
            vmin = space
            vmax = vmin + 6
            levels = np.arange(vmin-2, vmax + 0.5, step=0.5)
            # plt.rcParams.update({
            #     "text.usetex": True,
            #     "font.family": "sans",
            #     "font.sans": "Helvetica",
            # })
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
            
            # Create 3D plot
            fig3d = plt.figure(figsize=(4, 3))
            ax3d = fig3d.add_subplot(111, projection='3d')
            
            fig3d.subplots_adjust(left=0, right=1, bottom=0.05, top=1.1)
            
            ax3d.set_box_aspect((1.2, 1.2, 1.2))
            hist = H.copy()
            hist /= np.max(hist)
            Z = np.ma.masked_where(hist <= 0, -np.log(hist))
            zmin = np.min(Z)
            print(zmin)
            space = 1.5
            # Z *= 2
            Z_threshold = 8.5 + space
            Z_shifted = Z + space
            Z_masked = np.ma.masked_where(Z_shifted > Z_threshold, Z_shifted)
            surf = ax3d.plot_surface(X, Y, Z_masked, cmap='coolwarm', edgecolor='k',
                                   lw=0.1, antialiased=True, shade=True,
                                   vmin=space, vmax=10+space)
            ax3d.contour(X, Y, Z, zdir='z', offset=0, color='k')
            ax3d.contourf(X, Y, Z, zdir='z', offset=0, cmap='coolwarm')
            ax3d.xaxis.pane.set_facecolor('none')
            ax3d.yaxis.pane.set_facecolor('none')
            ax3d.zaxis.pane.set_facecolor('none')
            ax3d.xaxis.pane.set_edgecolor('none')
            ax3d.yaxis.pane.set_edgecolor('none')
            ax3d.zaxis.pane.set_edgecolor('none')
            ax3d.set_zticks([])
            ax3d.zaxis.line.set_visible(False)
            ax3d.zaxis.label.set_visible(False)
            ax3d.grid(False)
            ax3d.set_xlim([19, 37])
            ax3d.set_ylim([17, 25])
            #ax3d.set_zlim([0, 9])
            ax3d.set_xticks([20, 25, 30, 35])
            ax3d.set_yticks([18, 20, 22, 24])
            ax3d.tick_params(axis='both', which='major', length=4, pad=-2)
            ax3d.tick_params(axis='x', labelsize=10, labelfontfamily='Arial')
            ax3d.tick_params(axis='y', labelsize=10, labelfontfamily='Arial')
            ax3d.set_xlabel(r'R$_{\mathrm{\mathsf{LID-CORE}}}(\mathrm{\AA})$', fontsize=12, labelpad=-5.5)
            ax3d.set_ylabel(r'R$_{\mathrm{\mathsf{NMP-CORE}}}(\mathrm{\AA})$', fontsize=12, labelpad=-5.5)
            
            ax3d.view_init(elev=20, azim=-40)
            ax3d.dist = 8
            
            cbar_ax = fig3d.add_axes([0.73, 0.35, 0.02, 0.35])  # [left, bottom, width, height]
            cbar = fig3d.colorbar(surf, cax=cbar_ax)
            cbar.set_ticks(np.asarray([0, 3, 6, 9]) + space)
            cbar.set_ticklabels([0, 3, 6, 9], size=8, font='Arial')
            cbar.set_label(r"$\mathrm{\mathsf{k_B T}}$", fontsize=10, labelpad=2)
            cbar.ax.yaxis.set_ticks_position('right') 
            cbar.ax.yaxis.set_label_position('right')  
            cbar.ax.tick_params(direction='out', length=3, width=1, colors='k', pad=3)  
            output_filename_3d = fign.replace('.png', '_3d.png')
            plt.savefig(output_filename_3d, dpi=600)
            plt.close(fig3d)
            print(output_filename_3d)
    cb_ax = fig.add_axes([0.865, 0.16, 0.015, 0.4])
    cb_ax.spines['top'].set_linewidth(10)
    cb_ax.spines['top'].set_color('red')
    cbar = fig.colorbar(cntr, cax=cb_ax, ticks=np.arange(0, 8, 2))
    cbar.ax.tick_params(axis='y', direction='out', labelsize=10, length=2)
    fig.text(0.875, 0.595, '$\mathrm{k_{B}T}$', fontsize=12)
    fig.savefig(fign, dpi=my_dpi)
    plt.close("all")