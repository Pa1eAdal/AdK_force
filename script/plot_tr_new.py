# -*- codeing = utf-8 -*-
# Time : 2022/9/18 22:13
# File : plot_tr_new.py
# Software : PyCharm
"""
Input: .npz file from ake_efficient.py and .txt file from ake_efficient.F90
"""
import itertools
import numpy as np
from matplotlib import pyplot as plt
# import proplot as plt
from matplotlib import rc
from matplotlib.lines import Line2D
import os
import method
from scipy.optimize import curve_fit
from time import localtime
import matplotlib.colors as mcolors

control = "force_c"

if "force" in control:
    project_name = "akesi_force"
elif "dna" in control:
    project_name = "akesi_dna"
else:
    project_name = "akesi"

task_name = "turnover_rate" + "_" + control
if control == "force_c":
    dv_list = [14, 16]
    # dv_list = [0]
elif control == "dna_c":
    dv_list = [0, 14]
atp_list = [1000]
amp_list = [9, 22, 44, 66, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000]
adp_list = [30]
cut_list = ["rmsd"]
conc_list = [30, 50, 100, 200, 300, 500, 1000, 2000, 3000]
site = [22, 42, 55, 73, 148, 177]
site_list = [148]
lod_list = [40, 60, 70]
kof_list = [0.1]
cf_list = [4]
cb_list = [5]
cb = 5
ba_list = [[1, 1]]
io_list = [0.05]
dp_list = [0]

# k_list = [0.01]
# rp_list = [3, 12]
rmsd = "0.8"
pale = 0
chi = 0.5
n_seed = 20
special = f"pic_{localtime().tm_year % 100}_{localtime().tm_mon}_{localtime().tm_mday}"
data_dir = "txtdata/"

color_list = list(mcolors.TABLEAU_COLORS.keys())
cf_color_dict = dict(zip(cf_list, color_list))

if pale == 1:
    p1 = 0.11
    p7 = 0.41
elif pale == 0:
    p1 = 0.09
    p7 = 0.36

my_dpi = 350

rect = [0.175, 0.14, 0.813, 0.85]
if "both" in special:
    rect = [0.06, 0.11, 0.9, 0.8]
    y_dict = dict()


def rdata(file, ctrl):
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
        tr = tr * 2e5 / 112 / 2
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
    if ctrl == 3:
        # random.choice 20 data * 20
        avelist = method.sampling_r(tr)
        ave = np.mean(tr)
        aar.append(ave)
        rmsd = np.std(avelist)
        aar.append(rmsd)

    return aar


def funcsi(conc, vm, km, ki):
    return vm * conc / (km + conc + conc * conc / ki)


def funcnosi(conc, vm, km):
    return vm * conc / (km + conc)


# main
plot_setting = method.Plot_setting()

position_list = ["LID(148)-CORE(22)", "NMP(42)-LID(144)", "NMP(55)-CORE(166)", "NMP(73)-LID(142)", "LID(148)-CORE(177)",
                 "NMP(55)-CORE(177)"]
position_dict = dict(zip(site, position_list))

dvlist = [0, 1, 2, 10, 14, 15, 16]
pcloselist = [0.28, 0.58, 0.78, 0.12, 0.99, 0.02, 0.01]
pclosedict = dict(zip(dvlist, pcloselist))

outdir = "./Figures/tr/pdf"
out_data_dir = './value'

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)
method.mk_outdir(out_data_dir)

if control == "force_c":
    out_data_name = f'{out_data_dir}/turnover_rate.dat'
    with open(out_data_name, 'w') as data_file:
        for dv, s in itertools.product(dv_list, site_list):
            if (dv == 0) & (s == 148):
                force_list = [7, 10, 15, 20]
                force_im_list = [10, 15, 20]
                force_color_list = ['tab:purple', "tab:red", "tab:blue", 'tab:green']
            elif (dv == 0) & (s == 177):
                force_list = [2, 5, 10, 15]
                force_im_list = [5, 15]
                force_color_list = ["tab:red", "tab:blue", 'tab:green', "tab:orange"]
            elif (dv == 0) & (s == 73):
                force_list = [2, 5, 10, 15]
                force_im_list = [2, 10]
                force_color_list = ["tab:red", "tab:blue", 'tab:green', "tab:orange"]
            elif dv == 16:
                force_list = [2, 5, 10]
                force_im_list = force_list
                force_color_list = ["tab:red", 'tab:purple', "tab:blue", 'tab:green', 'm', 'yellow', 'brown', "palegreen"]
            elif dv == 14:
                force_list = [10, 15, 20]
                force_im_list = force_list
                force_color_list = ["tab:blue", 'tab:green', 'tab:brown', "tab:red", "palegreen"]
            force_color_dict = dict(zip(force_list, force_color_list))
            picture_name = f"{project_name}_{task_name}_dv{dv}_rm{rmsd}_pi{pale}_s{s}"
            pn = special + "_" + picture_name
            fign = outdir.strip("pdf") + pn + ".png"
            pdfn = outdir + "/" + pn + ".pdf"
            fig = plt.figure(figsize=(2.25, 2), dpi=my_dpi)
            ax = fig.add_axes(rect)
            ax.set_ylabel('Turnover rate(ms$^{-1}$)')
            ax.set_xlabel('[AMP]($\mu$M)')
            ax.set_xlim([-100, 3200])
            ax.set_ylim([0, 0.52])
            ax.set_yticks(np.arange(0, 0.6, 0.1))

            for f in force_list:
                if (f in [10, 15]) & (s == 148) & (dv == 0):
                    alpha = 1
                    zorder = 100 - f
                    lw = 2
                elif (f in [2, 5, 10]) & (s == 148) & (dv == 16):
                    alpha = 1
                    zorder = 100 - f
                    lw = 2
                elif (f in [5, 10, 15, 20]) & (s == 148) & (dv == 14):
                    alpha = 1
                    zorder = 100 - f
                    lw = 2
                elif (f in [2, 10]) & (s == 73):
                    alpha = 1
                    zorder = 100 - f
                    lw = 2
                elif (f in [5, 15]) & (s == 177):
                    alpha = 1
                    zorder = 100 - f
                    lw = 2
                else:
                    alpha = 1
                    lw = 2
                    zorder = f
                y = []

                for co in conc_list:
                    c_list = []
                    c = str(co).zfill(4)

                    # Write header for each condition
                    data_file.write(f"dv{dv}\tforce{f}\tsite{s}\t[AMP]{co}\n")
                    
                    for seed in range(n_seed):
                        n = f"{seed}".zfill(3)
                        filename = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.npz"
                        fn = data_dir + filename
                        data = np.load(fn)
                        turnover_rate = np.squeeze(data["turnover_rate"]) * 2e5 / 112 / 2
                        c_list.append(turnover_rate)
                        # Write data for each seed
                        data_file.write(f"n{n}\t{turnover_rate}\n")
                    
                    y.append(np.mean(c_list))
                y = np.asarray(y)

                if len(conc_list) > 4:
                    if dv == 17:
                        popt, pcov = curve_fit(funcnosi, conc_list, y, maxfev=800000)
                        vm = popt[0]
                        km = popt[1]
                        xlist = np.linspace(0, 3000, 300)
                        ylist = funcnosi(xlist, vm, km)
                        ax.plot(xlist, ylist, '-', color=force_color_dict[f], lw=lw, alpha=alpha, zorder=zorder)
                    else:
                        popt, pcov = curve_fit(funcsi, conc_list, y, maxfev=800000)
                        vm = popt[0]
                        km = popt[1]
                        ki = popt[2]
                        xlist = np.linspace(0, 3000, 300)
                        ylist = funcsi(xlist, vm, km, ki)
                        ax.plot(xlist, ylist, '-', color=force_color_dict[f], lw=lw, alpha=alpha, zorder=zorder)

                else:
                    ax.plot(conc_list, y, lw=lw, color=force_color_dict[f])
                scat = ax.scatter(conc_list, y, s=35, label=f"{f}pN", color=force_color_dict[f], alpha=alpha, zorder=zorder,
                                  edgecolor='none')

            if (s == 73) or (s == 177):
                ax.legend(loc="upper right", ncol=2, columnspacing=0.1, handletextpad=0.05, borderaxespad=0.03)

            # background
            ctrl = 0
            nod = ctrl * 40 + 20
            b_conc_list = [10, 30, 70, 100, 200, 300, 500, 1000, 2000, 3000]
            b_f = 0
            b_alpha = 0.7
            b_zorder = 100
            b_y = []
            
            # Write background header
            data_file.write(f"\n# Background data\n")
            
            for co in b_conc_list:
                b_c = str(co).zfill(4)
                if co == 2000:
                    tr_list = []
                    # Write header for this concentration
                    data_file.write(f"dv{dv}\tforce{b_f}\t[AMP]{co}\tsite{s}\n")
                    for seed in range(n_seed):
                        n = f"{seed}".zfill(3)
                        fn = os.path.expanduser(f"~/soft/cafemol-ake-7_rmsd/data/txtdata/akesi_force_dv{dv}_rm{rmsd}_pi{pale}_c{b_c}_s148_f{b_f}_n{n}.npz")
                        data = np.load(fn)
                        turnover_rate = np.squeeze(data["turnover_rate"])
                        tr_list.append(turnover_rate)
                        # Write background data for 2000
                        data_file.write(f"n{n}\t{turnover_rate * 2e5 / 112 / 2}\n")
                    b_y.append(np.mean(tr_list) * 2e5 / 112 / 2)
                else:
                    b_filename = os.path.expanduser(f"~/soft/cafemol-ake-7_rmsd/data/txtdata/rakesi{dv}_c{b_c}_f{b_f}.txt")
                    with open(b_filename, 'r') as bf:
                        tr_list = []
                        # Write header for this concentration
                        data_file.write(f"dv{dv}\tforce{b_f}\t[AMP]{co}\tsite{s}\n")
                        for line in bf:
                            # Each line has format: n[cycle] [turnover_rate]
                            values = line.strip().split()
                            turnover_rate = float(values[-1])
                            tr_list.append(turnover_rate)
                            # Write background data with cycle number
                            data_file.write(f"n{len(tr_list):03d}\t{turnover_rate}\n")
                    b_y.append(np.mean(tr_list))

            if dv == 17:
                popt, pcov = curve_fit(funcnosi, b_conc_list, b_y, maxfev=800000)
                vm = popt[0]
                km = popt[1]
                xlist = np.linspace(0, 3000, 300)
                ylist = funcnosi(xlist, vm, km)
            else:
                popt, pcov = curve_fit(funcsi, b_conc_list, b_y, maxfev=800000)
                vm = popt[0]
                km = popt[1]
                ki = popt[2]
                xlist = np.linspace(0, 3000, 300)
                ylist = funcsi(xlist, vm, km, ki)
            ax.plot(xlist, ylist, '-', color='k', lw=2, alpha=b_alpha, zorder=b_zorder)
            ax.plot(b_conc_list, b_y, 'o', ms=6.5, mew=0, label=f'0pN', color='k', alpha=b_alpha, zorder=b_zorder)
            if (s == 73) or (s == 177):
                pass
            else:
                ax.legend(loc="best", ncol=2, columnspacing=0.1, handletextpad=0.05, borderaxespad=0.03)
            fig.savefig(fign, dpi=my_dpi)
            print(fign)
if control == "dna_c":
    out_data_name = f'{out_data_dir}/turnover_rate_dna.dat'
    with open(out_data_name, 'w') as data_file:
        lod_color_list = ["tab:green", 'tab:blue', 'tab:red', 'yellow', 'brown']
        lod_color_dict = dict(zip(lod_list, lod_color_list))
        for dv, s, dp, io in itertools.product(dv_list, site_list, dp_list, io_list):
            picture_name = f"{project_name}_{task_name}_dv{dv}_rm{rmsd}_pi{pale}_s{s}_dp{dp}_io{io}"
            pn = special + "_" + picture_name
            fign = outdir.strip("pdf") + pn + ".png"
            pdfn = outdir + "/" + pn + ".pdf"
            fig = plt.figure(figsize=(2.25, 2), dpi=my_dpi)
            ax = fig.add_axes(rect)
            ax.set_ylabel('Turnover rate(ms$^{-1}$)')
            ax.set_xlabel('[AMP]($\mu$M)')
            ax.set_xlim([-100, 3200])
            ax.set_ylim([0, 0.55])
            ax.set_yticks(np.arange(0, 0.6, 0.1))

            for lod in lod_list:
                alpha = 1
                zorder = 100 - lod
                y = []
                if (lod == 60) & (dv == 0):
                    conc_list = [30, 50, 200, 300, 500, 1000, 2000, 3000]
                else:
                    conc_list = [30, 50, 100, 200, 300, 500, 1000, 2000, 3000]

                for co in conc_list:
                    c_list = []
                    c = str(co).zfill(4)

                    # Write header for each condition
                    data_file.write(f"dv{dv}\tlod{lod}\tsite{s}\t[AMP]{co}\n")

                    for seed in range(n_seed):
                        n = f"{seed}".zfill(3)
                        filename = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}_n{n}.npz"
                        fn = data_dir + filename
                        data = np.load(fn)
                        turnover_rate = np.squeeze(data["turnover_rate"]) * 2e5 / 112 / 2
                        c_list.append(turnover_rate)
                        # Write data for each seed
                        data_file.write(f"n{n}\t{turnover_rate}\n")

                    y.append(np.mean(c_list))
                y = np.asarray(y)

                if len(conc_list) > 3:
                    popt, pcov = curve_fit(funcsi, conc_list, y, maxfev=800000)
                    vm = popt[0]
                    km = popt[1]
                    ki = popt[2]
                    xlist = np.linspace(0, 3000, 300)
                    ylist = funcsi(xlist, vm, km, ki)
                    ax.plot(xlist, ylist, '-', color=lod_color_dict[lod], lw=2, alpha=alpha, zorder=zorder)
                else:
                    ax.plot(conc_list, y, lw=2, color=lod_color_dict[lod])
                ax.scatter(conc_list, y, s=25, label=f"{lod}bp", color=lod_color_dict[lod], edgecolor='none', zorder=zorder)

            # background
            ctrl = 0
            nod = ctrl * 40 + 20
            b_conc_list = [10, 30, 70, 100, 200, 300, 500, 1000, 2000, 3000]
            b_f = 0
            b_alpha = 0.6
            b_zorder = 1
            b_y = []
            for co in b_conc_list:
                b_c = str(co).zfill(4)
                if co == 2000:
                    tr_list = []
                    for seed in range(n_seed):
                        n = f"{seed}".zfill(3)
                        fn = os.path.expanduser(f"~/soft/cafemol-ake-7_rmsd/data/txtdata/akesi_force_dv{dv}_rm{rmsd}_pi{pale}_c{b_c}_s148_f{b_f}_n{n}.npz")
                        data = np.load(fn)
                        turnover_rate = np.squeeze(data["turnover_rate"])
                        tr_list.append(turnover_rate)
                    b_y.append(np.mean(tr_list) * 2e5 / 112 / 2)
                else:
                    b_filename = os.path.expanduser(f"~/soft/cafemol-ake-7_rmsd/data/txtdata/rakesi{dv}_c{b_c}_f{b_f}.txt")
                    aar = rdata(b_filename, ctrl)
                    b_y.append(aar[0])
            popt, pcov = curve_fit(funcsi, b_conc_list, b_y, maxfev=800000)
            vm = popt[0]
            km = popt[1]
            ki = popt[2]
            xlist = np.linspace(0, 3000, 300)
            ylist = funcsi(xlist, vm, km, ki)
            ax.plot(xlist, ylist, '-', color='k', lw=2, alpha=b_alpha, zorder=b_zorder)
            ax.plot(b_conc_list, b_y, 'o', ms=6, mew=0, label=f'free', color='k', alpha=b_alpha, zorder=b_zorder)

            ax.legend(loc="upper right", ncol=2, columnspacing=0.1, handletextpad=0.05, borderaxespad=0.03)
            fig.savefig(fign, dpi=my_dpi)
            print(fign)

if control == "atp":
    dv_list = [0, 14, 16]
    atp_list = [10, 30, 100, 300, 1000, 3000]
    amp_list = [300]
    adp_list = [0]
    lw = 2
    alpha = 1
    zorder = 10
    for amp, ba, cf, cb, cut, adp in itertools.product(amp_list, ba_list, cf_list, cb_list, cut_list, adp_list):
        if cut == "rmsd":
            cut_val = 0.8
        m = f"{amp:04d}"
        d = f"{adp:04d}"
        bas = ba[0]
        bap = ba[1]
        dv_color_list = ["black", "tab:blue", "tab:red"]
        dv_color_dict = dict(zip(dv_list, dv_color_list))
        picture_name = f"{control}_{project_name}_{cut}{cut_val}_p1_{p1}_p7_{p7}_bas{bas}_bap{bap}_cf{cf}_cb{cb}_m{m}_d{d}"
        pn = special + "_" + picture_name
        fign = outdir.strip("pdf") + pn + ".png"
        pdfn = outdir + "/" + pn + ".pdf"
        fig = plt.figure(figsize=(2.25, 2), dpi=my_dpi)
        ax = fig.add_axes(rect)
        # ax.set_title(f"$P_{{close}} = {pclosedict[dv]}$, {position_dict[s]}")
        ax.set_ylabel('Turnover rate(ms$^{-1}$)')
        ax.set_xlabel('[ATP]($\mu$M)')
        ax.set_xlim([-100, 3200])
        ax.set_ylim([0, 0.55])
        ax.set_yticks(np.arange(0, 0.6, 0.1))
        # ax.set_xscale('log')
        for dv in dv_list:
            y = []
            for co in atp_list:
                c_list = []
                c = str(co).zfill(4)
                for seed in range(n_seed):
                    n = f"{seed}".zfill(3)
                    filename = f"{project_name}_dv{dv}_{cut}{cut_val}_p1_{p1}_p7_{p7}_bas{bas}_bap{bap}_cf{cf}_cb{cb}_t{c}_m{m}_d{d}_n{n}.npz"
                    fn = data_dir + filename
                    data = np.load(fn)
                    turnover_rate = np.squeeze(data["turnover_rate"])
                    c_list.append(turnover_rate)
                y.append(np.mean(c_list))
            y = np.asarray(y)
            print(y)
            y *= 2e8 / 112 / 2
            print("=====")
            print(data["nframe"][-1])
            print(y)
            if len(atp_list) > 4:
                if dv == 17:
                    popt, pcov = curve_fit(funcnosi, atp_list, y, maxfev=800000)
                    vm = popt[0]
                    km = popt[1]
                    xlist = np.linspace(0, 3000, 300)
                    ylist = funcnosi(xlist, vm, km)
                    ax.plot(xlist, ylist, '-', color=dv_color_dict[dv], lw=lw, alpha=alpha, zorder=zorder)
                else:
                    popt, pcov = curve_fit(funcnosi, atp_list, y, maxfev=800000)
                    vm = popt[0]
                    km = popt[1]
                    xlist = np.linspace(0, 3000, 300)
                    ylist = funcnosi(xlist, vm, km)
                    ax.plot(xlist, ylist, '-', color=dv_color_dict[dv], lw=lw, alpha=alpha, zorder=zorder)

            else:
                ax.plot(atp_list, y, lw=lw, color=dv_color_dict[dv])
            scat = ax.scatter(atp_list, y, s=35, label=f"$P_{{close}}$={pclosedict[dv]}", color=dv_color_dict[dv],
                              alpha=alpha, zorder=zorder, edgecolor='none')

        ax.legend(loc="best", ncol=1, columnspacing=0.1, handletextpad=0.05, borderaxespad=0.03)

        fig.savefig(fign, dpi=my_dpi)
        print(fign)
