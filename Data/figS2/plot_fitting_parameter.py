# -*- codeing = utf-8 -*-
# Time : 2022/11/5 1:46
# File : plot_fit_tr.py
# Software : PyCharm
"""
Input: .npz file from ake_efficient.py and .txt file from ake_efficient.F90
"""
import itertools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import os
from scipy.optimize import curve_fit
import method
import enzyme_setting
from time import localtime

control = "force_c"
cut = 'rmsd'

if "force" in control:
    project_name = "akesi_force"
elif "dna" in control:
    project_name = "akesi_dna"
else:
    project_name = "akesi"

adk = enzyme_setting.Enzyme_adk(cut=cut, control=control)
dv_list = [14, 16]
conc_list = [30, 50, 100, 200, 300, 500, 1000, 2000, 3000]
site = [22, 42, 55, 73, 148, 177]
site_list = [148]
lod_list = [70, 60, 40]
io_list = [0.05]
dp_list = [0]
force_list = [2, 3, 5, 7, 10, 15, 20]
# k_list = [0.01]
# rp_list = np.asarray([-2, 2])
rmsd = "0.8"
pale = 0
n_seed = 20
special = f"pic_{localtime().tm_year % 100}_{localtime().tm_mon}_{localtime().tm_mday}"
data_dir = "txtdata/"
fit_dict = {}
final_d = {}

my_dpi = 350
rect = [0.175, 0.14, 0.815, 0.85]
if "both" in special:
    rect = [0.06, 0.11, 0.9, 0.8]
    y_dict = dict()
width = 0.8


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
        ave = np.mean(avelist)
        aar.append(ave)
        rmsd = np.std(avelist)
        aar.append(rmsd)
    if ctrl == 4:
        # random.choice 20 data
        avelist = method.sampling_r(tr, times=1, r_seed=None)
        ave = np.mean(avelist)
        aar.append(ave)
        rmsd = np.std(avelist)
        aar.append(rmsd)
    return aar


def funcsi(conc, vm, km, ki):
    return vm * conc / (km + conc + conc * conc / ki)


# main
plot_setting = method.Plot_setting()
rect = [0.175, 0.14, 0.815, 0.85]

task_name = control

position_list = ["LID(148)-CORE(22)", "NMP(42)-LID(144)", "NMP(55)-CORE(166)", "NMP(73)-LID(142)", "LID(148)-CORE(177)",
                 "NMP(55)-CORE(177)"]
position_dict = dict(zip(site, position_list))

dvlist = [0, 1, 2, 14, 16]
pcloselist = [0.28, 0.58, 0.78, 0.99, 0.01]
pclosedict = dict(zip(dvlist, pcloselist))

lod_dis_list = [3, 2, 1]
lod_dis_dict = dict(zip(lod_list, lod_dis_list))

cycle_end = [cycle[i + 1] for i in range(len(cycle) - 1)]
cycle_end.append(cycle[0])

kma_dir = './kmatrix/'
outdir = "./Figures/tr/pdf"
method.mk_outdir(outdir)

if control == "force_c":
    fit_n_list = ['vm', 'km', 'ki']
    fit_n_label = ['$k_{cat}(ms^{-1}$)', '$K_{M}(\mu M$)', '$K_{I}(mM$)']
    fit_n_dict = dict(zip(fit_n_list, fit_n_label))
    fitn_color_list = ["tab:blue", "tab:blue", "tab:blue"]
    fitn_color_dict = dict(zip(fit_n_list, fitn_color_list))
    for dv in dv_list:
        force_setting = enzyme_setting.Force_setting()
        fit_dict[f'vm'] = []
        fit_dict[f'km'] = []
        fit_dict[f'ki'] = []
        fit_dict[f'vmr'] = []
        fit_dict[f'kmr'] = []
        fit_dict[f'kir'] = []
        if dv == 0:
            pass
        if dv == 16:
            force_list = [1, 2, 3, 5, 10]
        elif dv == 14:
            force_list = [5, 10, 15, 20, 25]

        for f in force_list:
            y = []
            fit_dict[f'vm{f}'] = []
            fit_dict[f'km{f}'] = []
            fit_dict[f'ki{f}'] = []
            if f in [2, 5, 10, 15, 20, 25]:
                conc_list = [30, 50, 100, 200, 300, 500, 1000, 2000, 3000]
            else:
                conc_list = [30, 100, 300, 1000, 3000]
            for co in conc_list:
                c = str(co).zfill(4)
                y_list = []

                for seed in range(n_seed):
                    n = f"{seed}".zfill(3)
                    filename = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.npz"
                    fn = data_dir + filename
                    data = np.load(fn)
                    turnover_rate = np.squeeze(data["turnover_rate"])
                    y_list.append(turnover_rate)
                final_list = method.sampling_r(y_list)
                y.append(final_list)
            y = np.asarray(y)
            y *= 2e5 / 112 / 2

            if len(conc_list) > 4:
                for i in range(len(final_list)):
                    popt, pcov = curve_fit(funcsi, conc_list, y[:, i], maxfev=800000)
                    vm = popt[0]
                    km = popt[1]
                    ki = popt[2]
                    fit_dict[f'vm{f}'].append(vm)
                    fit_dict[f'km{f}'].append(km)
                    fit_dict[f'ki{f}'].append(ki)

            # final_vm = method.sampling_r(fit_dict[f'vm{f}'])
            # final_km = method.sampling_r(fit_dict[f'km{f}'])
            # final_ki = method.sampling_r(fit_dict[f'ki{f}'])

            fit_dict[f'vm'].append(np.average(fit_dict[f'vm{f}']))
            fit_dict[f'km'].append(np.average(fit_dict[f'km{f}']))
            fit_dict[f'ki'].append(np.average(fit_dict[f'ki{f}']))
            fit_dict[f'vmr'].append(np.std(fit_dict[f'vm{f}']))
            fit_dict[f'kmr'].append(np.std(fit_dict[f'km{f}']))
            fit_dict[f'kir'].append(np.std(fit_dict[f'ki{f}']))

        # background
        fit_dict[f'vm0'] = []
        fit_dict[f'km0'] = []
        fit_dict[f'ki0'] = []
        ctrl = 4
        nod = ctrl * 40 + 20
        b_conc_list = [10, 30, 100, 200, 300, 500, 800, 1000, 3000]
        b_conc_list2 = [200, 500, 800]
        b_f = 0
        b_alpha = 0.5
        b_zorder = 1
        y = []
        for i in range(20):
            b_y = []
            for co in b_conc_list:
                b_c = str(co).zfill(4)

                b_filename = f"/fsa/home/ww_zhangzy/soft/cafemol-ake-7_rmsd/data/txtdata/rakesi8{dv}_c{b_c}_s42_f{b_f}.txt"
                if co in b_conc_list2:
                    b_c = str(co).zfill(4)
                    b_filename = f"/fsa/home/ww_zhangzy/soft/cafemol-ake-7_rmsd/data/txtdata/rakesi8{dv}_p0_c{b_c}_s42_f{b_f}.txt"
                aar = rdata(b_filename, ctrl)
                b_y.append(aar[0])
            y.append(b_y)
        y = np.asarray(y)
        for j in range(len(b_y)):
            popt, pcov = curve_fit(funcsi, b_conc_list, y[j, :], maxfev=800000)
            vm = popt[0]
            km = popt[1]
            ki = popt[2]
            fit_dict[f'vm0'].append(vm)
            fit_dict[f'km0'].append(km)
            fit_dict[f'ki0'].append(ki)

        fit_dict[f'vm0r'] = np.std(fit_dict[f'vm0'])
        fit_dict[f'km0r'] = np.std(fit_dict[f'km0'])
        fit_dict[f'ki0r'] = np.std(fit_dict[f'ki0'])

        force_list.insert(0, 0)
        for fit_n in ['vm', 'km', 'ki']:

            picture_name = f"{project_name}_{task_name}_dv{dv}_rm{rmsd}_pi{pale}_s{s}_{fit_n}"
            pn = special + "_" + picture_name
            fign = outdir.strip("pdf") + pn + ".png"
            pdfn = outdir + "/" + pn + ".pdf"
            fig = plt.figure(figsize=(2.25, 2), dpi=my_dpi)
            #if fit_n == "ki":
            #    rect = [0.143, 0.13, 0.85, 0.85]
            ax = fig.add_axes(rect)
            if fit_n == "vm":
                if dv == 0:
                    ax.set_ylim([0.35, 0.85])
                    ax.set_yticks([0.4, 0.5, 0.6, 0.7, 0.8])
                if dv == 14:
                    ax.set_yticks(np.arange(0.5, 1.2, 0.2))
                # ax.set_title("$V_{m}$")
            if fit_n == "km":
                if dv == 0:
                    ax.set_ylim([0, 280])
                # ax.set_title("$K_{M}$")
            if fit_n == "ki":
                if dv == 0:
                    ax.set_ylim([0, 3.2])
                    ax.set_yticks(np.arange(0, 3.2, 0.5))
                if dv == 14:
                    ax.set_yticks(np.arange(0.5, 0.9, 0.1))
            if dv == 16:
                ax.set_xticks(np.arange(0, 11, 2))
            elif dv == 14:
                ax.set_xticks(np.arange(0, 26, 5))
                # ax.set_title("$K_{i}$")
            ax.set_ylabel(f'{fit_n_dict[fit_n]}')
            ax.set_xlabel(f'Force(pN)')
            # ax.get_xaxis().set_visible(False)

            # ax.plot(0, fit_dict[f'{fit_n}0'], '--', color='k', lw=3.5, alpha=b_alpha, zorder=b_zorder)
            # ax.plot(0, fit_dict[f'{fit_n}0'], 'o', ms=10, mew=0, color='k', alpha=b_alpha,
            #         zorder=b_zorder)
            # ax.axhline(y=fit_dict[f'{fit_n}0'], ls='--', lw=3, alpha=0.5, color='k')
            # plt.show()  # 画出图像
            # ax.legend(loc="best")

            fit_dict[f'{fit_n}'].insert(0, np.average(fit_dict[f'{fit_n}0']))
            fit_dict[f'{fit_n}r'].insert(0, fit_dict[f'{fit_n}0r'])
            if fit_n == "ki":
                fit_dict[f'{fit_n}'] = np.asarray(fit_dict[f'{fit_n}'])
                fit_dict[f'{fit_n}r'] = np.asarray(fit_dict[f'{fit_n}r'])
                fit_dict[f'{fit_n}'] /= 1000
                fit_dict[f'{fit_n}r'] /= 1000
            # print(force_list, '\n', fit_dict[f'{fit_n}'])

            ax.errorbar(force_list, fit_dict[f'{fit_n}'], yerr=fit_dict[f'{fit_n}r'], ms=7, capsize=3.8,
                        color=fitn_color_dict[fit_n])
            ax.scatter(force_list, fit_dict[f"{fit_n}"], s=35, color="none", edgecolor='tab:red', zorder=100)

            fig.savefig(fign, dpi=my_dpi)
            print(fign)

if control == "dna_c":
    dna_color_list = ['tab:green', 'tab:blue', "tab:red"]
    dna_color_dict = dict(zip(lod_list, dna_color_list))
    fit_n_list = ['vm', 'km', 'ki']
    fit_n_label = ['$k_{cat}(ms^{-1}$)', '$K_{M}(\mu M$)', '$K_{I}(\mu M$)']
    fit_n_dict = dict(zip(fit_n_list, fit_n_label))
    # posi_color_list = ["r", "k", 'm', 'yellow', "green", 'brown']
    # posi_color_dict = dict(zip(site, posi_color_list))
    for dv, s, io, dp in itertools.product(dv_list, site_list, io_list, dp_list):
        fit_dict[f'vm'] = []
        fit_dict[f'km'] = []
        fit_dict[f'ki'] = []
        fit_dict[f'vmr'] = []
        fit_dict[f'kmr'] = []
        fit_dict[f'kir'] = []
        for lod in lod_list:

            fit_dict[f'vm{lod}'] = []
            fit_dict[f'km{lod}'] = []
            fit_dict[f'ki{lod}'] = []
            fit_dict[f'vm0'] = []
            fit_dict[f'km0'] = []
            fit_dict[f'ki0'] = []

            y = []

            for co in conc_list:
                c = str(co).zfill(4)
                y_list = []

                for seed in range(n_seed):
                    n = f"{seed}".zfill(3)
                    filename = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_dp{dp}_l{lod}_io{io}_n{n}.npz"
                    fn = data_dir + filename
                    data = np.load(fn)
                    turnover_rate = np.squeeze(data["turnover_rate"])
                    y_list.append(turnover_rate)
                final_list = method.sampling_r(y_list)
                y.append(final_list)
            y = np.asarray(y)
            y *= 2e5 / 112 / 2

            if len(conc_list) > 4:
                for i in range(len(final_list)):
                    popt, pcov = curve_fit(funcsi, conc_list, y[:, i], maxfev=800000)
                    vm = popt[0]
                    km = popt[1]
                    ki = popt[2]
                    fit_dict[f'vm{lod}'].append(vm)
                    fit_dict[f'km{lod}'].append(km)
                    fit_dict[f'ki{lod}'].append(ki)
                    
                    
            fit_dict[f'vm'].append(np.average(fit_dict[f'vm{lod}']))
            fit_dict[f'km'].append(np.average(fit_dict[f'km{lod}']))
            fit_dict[f'ki'].append(np.average(fit_dict[f'ki{lod}']))
            fit_dict[f'vmr'].append(np.std(fit_dict[f'vm{lod}']))
            fit_dict[f'kmr'].append(np.std(fit_dict[f'km{lod}']))
            fit_dict[f'kir'].append(np.std(fit_dict[f'ki{lod}']))

            # background
            fit_dict[f'vm0'] = []
            fit_dict[f'km0'] = []
            fit_dict[f'ki0'] = []
            ctrl = 4
            nod = ctrl * 40 + 20
            b_conc_list = [10, 30, 100, 200, 300, 500, 800, 1000, 3000]
            b_conc_list2 = [200, 500, 800]
            b_f = 0
            b_alpha = 0.5
            b_zorder = 1
            y = []
            for i in range(20):
                b_y = []
                for co in b_conc_list:
                    b_c = str(co).zfill(4)

                    b_filename = f"/fsa/home/ww_zhangzy/soft/cafemol-ake-7_rmsd/data/txtdata/rakesi8{dv}_c{b_c}_s42_f{b_f}.txt"
                    if co in b_conc_list2:
                        b_c = str(co).zfill(4)
                        b_filename = f"/fsa/home/ww_zhangzy/soft/cafemol-ake-7_rmsd/data/txtdata/rakesi8{dv}_p0_c{b_c}_s42_f{b_f}.txt"
                    aar = rdata(b_filename, ctrl)
                    b_y.append(aar[0])
                y.append(b_y)
            y = np.asarray(y)
            for j in range(len(b_y)):
                popt, pcov = curve_fit(funcsi, b_conc_list, y[j, :], maxfev=800000)
                vm = popt[0]
                km = popt[1]
                ki = popt[2]
                fit_dict[f'vm0'].append(vm)
                fit_dict[f'km0'].append(km)
                fit_dict[f'ki0'].append(ki)

            fit_dict[f'vm0r'] = np.std(fit_dict[f'vm0'])
            fit_dict[f'km0r'] = np.std(fit_dict[f'km0'])
            fit_dict[f'ki0r'] = np.std(fit_dict[f'ki0'])

        for fit_n in ['vm', 'km', 'ki']:
            fit_dict[f'{fit_n}'].insert(0, np.average(fit_dict[f'{fit_n}0']))
            fit_dict[f'{fit_n}r'].insert(0, fit_dict[f'{fit_n}0r'])

            picture_name = f"{project_name}_{task_name}_dv{dv}_rm{rmsd}_pi{pale}_s{s}_dp{dp}_io{io}_{fit_n}"
            pn = special + "_" + picture_name
            fign = outdir.strip("pdf") + pn + ".png"
            pdfn = outdir + "/" + pn + ".pdf"
            fig = plt.figure(figsize=(2.25, 2), dpi=my_dpi)
            ax = fig.add_axes(rect)
            ax.set_ylabel(f'{fit_n_dict[fit_n]}')
            # ax.get_xaxis().set_visible(False)
            if fit_n == "vm":
                if dv == 0:
                    ax.set_ylim([0.35, 0.85])
                elif dv == 14:
                    ax.set_ylim([0.25, 0.65])
                    ax.set_yticks(np.arange(0.3, 0.7, 0.1))
                # ax.set_title("$V_{m}$")
            if fit_n == "km":
                if dv == 0:
                    ax.set_ylim([0, 280])
                elif dv == 14:
                    ax.set_ylim([0, 160])
                    ax.set_yticks(np.arange(0, 160, 30))
                # ax.set_title("$K_{M}$")
            if fit_n == "ki":
                if dv == 0:
                    # ax.set_ylim([0.9, 2.85])
                    ax.set_yticks(np.arange(1.2, 2.8, 0.5))
                elif dv == 14:
                    ax.set_ylim([0.48, 0.82])
                    ax.set_yticks(np.arange(0.5, 0.9, 0.1))
                # ax.set_title("$K_{i}$")
                fit_dict[f'{fit_n}'] = np.asarray(fit_dict[f'{fit_n}'])
                fit_dict[f'{fit_n}r'] = np.asarray(fit_dict[f'{fit_n}r'])
                fit_dict[f'{fit_n}'] /= 1000
                fit_dict[f'{fit_n}r'] /= 1000
            ax.set_xlim([-0.25, 3.25])

            # ax.set_xlabel('Separation$\AA$')
            #     ax.set_xticks([1, 2, 3], ["40", '60', '70'])
            # for lod in lod_list:
            # ax.bar(lod_dis_dict[lod], fit_dict[f'{fit_n}{lod}'][0], width=width, color=dna_color_dict[lod])
            # ax.plot([0, 1, 2, 3], fit_dict[f'{fit_n}{lod}'], lw=1, color=dna_color_dict[lod])

            xlist = np.arange(0, 4)
            ax.errorbar(xlist, fit_dict[f"{fit_n}"], yerr=fit_dict[f'{fit_n}r'], ms=7, capsize=3.8,
                        color="tab:blue")
            ax.scatter(xlist, fit_dict[f"{fit_n}"], s=35, color="none", edgecolor='tab:red', zorder=100)
            # ax.plot(xlist, ylist, '-', color='tab:blue', lw=3.5, alpha=b_alpha, zorder=b_zorder)
            # ax.bar(0, fit_dict[f'{fit_n}0'], width=width, color='k')
            # ax.axhline(y=fit_dict[f'{fit_n}0'], ls='--', lw=3, alpha=0.5, color='k')
            ax.set_xticks([3, 2, 1, 0], ["40bp", "60bp", "70bp", "Free"])
            # plt.show()  # 画出图像
            # ax.legend(loc="best")
            fig.savefig(fign, dpi=my_dpi)
            plt.close("all")
            print(fign)
