# -*- coding: utf-8 -*-
# Author : RuaKo
# Time : 2024/9/18 10:24
# File : plot_tr_force.py
# Software : PyCharm
"""
Input: .npz file from ake_efficient.py and .txt file from ake_efficient.F90
"""
import itertools
import numpy as np
from matplotlib import pyplot as plt
import method
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors
import enzyme_setting

control = "force_c"

if "force" in control:
    project_name = "akesi_force"
elif "dna" in control:
    project_name = "akesi_dna"
else:
    project_name = "akesi"

task_name = "turnover_rate" + "_" + control
adk = enzyme_setting.Enzyme_adk(cut="rmsd")
site_list = [148]
plot_setting = method.Plot_setting()
time = method.time_y_m_d()
adk.dv_list = [0, 14, 16]
adk.amp_list = [300, 1000, 3000]
rect = [0.165, 0.14, 0.815, 0.84]

funcsi = enzyme_setting.func_MM_si
funcnosi = enzyme_setting.func_MM

data_dir = "txtdata/"
outdir = "./Figures/tr_force/pdf"

method.mk_outdir(outdir)


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
        tr = tr * 2e5 / 112 /2
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


if control == "force_c":
    for dv, s in itertools.product(adk.dv_list, site_list):
        force_setting = enzyme_setting.Force_setting(dv)
        amp_color_list = ["tab:red", "tab:blue", 'tab:green', 'm', 'yellow', 'brown', "palegreen"]
        amp_color_dict = dict(zip(adk.amp_list, amp_color_list))
        picture_name = f"{project_name}_{task_name}_dv{dv}_rm{adk.cut_val}_pi{adk.pale}_s{s}"
        pn = time + "_" + picture_name
        fign = outdir.strip("pdf") + pn + ".png"
        pdfn = outdir + "/" + pn + ".pdf"
        fig = plt.figure(figsize=plot_setting.figsize_111, dpi=plot_setting.dpi)
        ax = fig.add_axes(rect)
        # ax.set_title(f"$P_{{close}} = {pclosedict[dv]}$, {position_dict[s]}")
        ax.set_ylabel('Turnover rate(ms$^{-1}$)')
        ax.set_xlabel('Force(pN)')
        if dv == 0:
            ax.set_xlim([-1, 21])
            ax.set_ylim([0, 0.52])
            ax.set_yticks(np.arange(0, 0.6, 0.1))
        elif dv == 16:
            ax.set_xlim([-0.5, 10.5])
            ax.set_ylim([0, 0.25])
            ax.set_yticks(np.arange(0, 0.3, 0.1))
            ax.set_xticks(np.arange(0, 11, 2))
        elif dv == 14:
            ax.set_xlim([-1, 26])
            ax.set_ylim([0, 0.62])
            ax.set_yticks(np.arange(0, 0.62, 0.1))
        # ax.set_xscale('log')

        for amp in adk.amp_list:
            m = f"{amp:04d}"
            alpha = 1
            zorder = 4000 - amp
            lw = 2
            y = []
            force_setting.force_list.remove(0)
            for f in force_setting.force_list:
                temp_list = []
                for seed in range(adk.n_seed):
                    n = f"{seed}".zfill(3)
                    filename = f"{project_name}_dv{dv}_rm{adk.cut_val}_pi{adk.pale}_c{m}_s{s}_f{f}_n{n}.npz"
                    fn = data_dir + filename
                    data = np.load(fn)
                    turnover_rate = np.squeeze(data["turnover_rate"])
                    temp_list.append(turnover_rate)
                y.append(np.mean(temp_list))

            ctrl = 0
            nod = ctrl * 40 + 20
            b_f = 0
            b_alpha = 0.7
            b_zorder = 100
            b_y = []
            if amp == 2000:
                tr_list = []
                for seed in range(adk.n_seed):
                    n = f"{seed}".zfill(3)
                    fn = os.path.expanduser(f"~/soft/cafemol-ake-7_rmsd/data/txtdata/akesi_force_dv{dv}_rm{adk.cut_val}_pi{adk.pale}_c{m}_s148_f{b_f}_n{n}.npz")
                    data = np.load(fn)
                    turnover_rate = np.squeeze(data["turnover_rate"])
                    tr_list.append(turnover_rate)
                b_y.append(np.mean(tr_list) * 2e5 / 112 / 2)
            else:
                b_filename = os.path.expanduser(f"~/soft/cafemol-ake-7_rmsd/data/txtdata/rakesi{dv}_c{m}_f{b_f}.txt")
                aar = rdata(b_filename, ctrl)
                b_y.append(aar[0])

            y.insert(0, b_y[0])
            y = np.asarray(y)
            y[1:] *= 2e5 / 112 / 2
            force_setting.force_list.insert(0, 0)
            ax.plot(force_setting.force_list, y, lw=lw, color=amp_color_dict[amp], zorder=zorder)
            ax.scatter(force_setting.force_list, y, s=35, label=f"{amp}$\mu$ M", color=amp_color_dict[amp], alpha=alpha,
                       zorder=zorder, edgecolor='none')

        # background
        ax.legend(loc="best", ncol=1, columnspacing=0.1, handletextpad=0.05, borderaxespad=0.03)
        fig.savefig(fign, dpi=plot_setting.dpi)
        print(fign)
