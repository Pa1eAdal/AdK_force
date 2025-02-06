# -*- codeing = utf-8 -*-
# Time : 2023/2/26 9:33
# File : plot_oo2cc.py
# Software : PyCharm
"""
Input: .npz file from calculate_distribution.py
"""
import itertools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import os
import method

control = "force_c"
special = method.time_y_m_d()
special2 = "no_confor_standard"


dv_list = [0, 14]
amp_list = [300, 3000]
site = [148]
# length_of_dna = [60]
# k_list = [0.01]
# rp_list = np.arange(0, 13, 2)
force_list = [0, 5, 7, 10, 15, 20]
n_seed = 20
rmsd = 0.8
pale = 0

project_name = "akesi_force"
work_dir = "./distribution/"
tr_dir = "./txtdata/"

plot_setting = method.Plot_setting()

all_alpha = 0.5
my_dpi = 350
rect = [0.15, 0.14, 0.82, 0.84]
figsize = (2.25, 2)

position_list = ["LID(148)-CORE(22)", "NMP(42)-LID(144)", "NMP(55)-CORE(166)", "NMP(73)-LID(142)"]
position_dict = dict(zip(site, position_list))

conc_color_list = ["tab:red", "tab:blue", 'b', 'm', 'yellow', 'brown']
conc_color_dict = dict(zip(amp_list, conc_color_list))

outdir = "./Figures/oo2cc/pdf"

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

dvn_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16]
p_close = [0.28, 0.52, 0.77, 0.15, 0.05, 0.02, 9e-3, 3e-3, 0.89, 0.94, 0.12, 0.10, 0.09, 0.986, 0.994, 0.008]
p_dict = dict(zip(dvn_list, p_close))

# Create output directory if it doesn't exist
out_data_dir = "./value"
if not os.path.exists(out_data_dir):
    os.makedirs(out_data_dir, exist_ok=True)

# Open data files
o2c_data_name = f'{out_data_dir}/AT_o2c.dat'
c2o_data_name = f'{out_data_dir}/AT_c2o.dat'
all_data_name = f'{out_data_dir}/AT_all.dat'
with open(o2c_data_name, 'w') as o2c_file, open(c2o_data_name, 'w') as c2o_file, open(all_data_name, 'w') as all_file:
    for dv, s in itertools.product(dv_list, site):
        fig_c2o = plt.figure(figsize=figsize, dpi=my_dpi)
        fig_o2c = plt.figure(figsize=figsize, dpi=my_dpi)
        ax_c2o = fig_c2o.add_axes(rect)
        ax_o2c = fig_o2c.add_axes(rect)
        if dv == 14:
            force_list = np.arange(0, 26, 5)
            ylim = (0, 16)
            ytick = np.arange(0, 16, 3)
            xtick = np.arange(0, 26, 5)
        elif dv == 0:
            ylim = (0, 7)
            ytick = np.arange(0, 8, 2)
            xtick = np.arange(0, 21, 5)

        for conc in amp_list:
            c = f"{conc}".zfill(4)
            c2o_ar_list = []
            c2o_rmsd_list = []
            o2c_ar_list = []
            o2c_rmsd_list = []
            all_ar_list = []
            all_rmsd_list = []
            zorder = conc
            for f in force_list:
                x_list_c2o = []
                x_list_o2c = []
                x_list_all = []
                data_n = f"{project_name}_{special2}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}.npz"

                np.load.__defaults__ = (None, True, True, 'ASCII')
                fn = f"{work_dir}" + data_n
                data = np.load(fn)
                np.load.__defaults__ = (None, False, True, 'ASCII')
                interval_list = []

                # Write header for each force/conc combination
                o2c_file.write(f"dv{dv}\tforce{f}\tsite{s}\t[AMP]{conc}\n")
                c2o_file.write(f"dv{dv}\tforce{f}\tsite{s}\t[AMP]{conc}\n")
                all_file.write(f"dv{dv}\tforce{f}\tsite{s}\t[AMP]{conc}\n")

                for n_s in range(n_seed):
                    n = f"{n_s:03d}"
                    tr_filename = f"{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.npz"
                    tr_fn = tr_dir + tr_filename
                    tr_data = np.load(tr_fn)
                    nframe_list = tr_data['nframe']
                    tem_list = nframe_list[::2]
                    for i in range(len(tem_list) - 1):
                        interval_list.append(tem_list[i + 1] - tem_list[i])

                c2o_list = data["c2o_list"]
                o2c_list = data["o2c_list"]
                all_list = np.asarray(interval_list)

                # Write raw data with conversion
                for n in range(len(o2c_list)):
                    o2c_file.write(f"{o2c_list[n] * 112 / 2e5}\n")
                for n in range(len(c2o_list)):
                    c2o_file.write(f"{c2o_list[n] * 112 / 2e5}\n")
                for n in range(len(all_list)):
                    all_file.write(f"{all_list[n] * 112 / 2e8}\n")

                c2o_list = c2o_list * 112 / 2e5
                o2c_list = o2c_list * 112 / 2e5
                all_list = all_list * 112 / 2e8

                o2c_tot_time = method.sampling_r(o2c_list, times=50, size=50)
                c2o_tot_time = method.sampling_r(c2o_list, times=50, size=50)
                all_tot_time = method.sampling_r(all_list, times=50, size=50)

                c2o_data = method.data_analysis(c2o_tot_time)
                o2c_data = method.data_analysis(o2c_tot_time)
                all_data = method.data_analysis(all_tot_time)

                c2o_ar_list.append(c2o_data[0])
                o2c_ar_list.append(o2c_data[0])
                all_ar_list.append(all_data[0])

                c2o_rmsd_list.append(c2o_data[1])
                o2c_rmsd_list.append(o2c_data[1])
                all_rmsd_list.append(all_data[1])

            ax_c2o.errorbar(force_list, c2o_ar_list, yerr=c2o_rmsd_list, capsize=2, color=conc_color_dict[conc],
                            zorder=zorder)
            ax_o2c.errorbar(force_list, o2c_ar_list, yerr=c2o_rmsd_list, capsize=2, color=conc_color_dict[conc],
                            zorder=zorder)
            ax_c2o.errorbar(force_list, all_ar_list, yerr=all_rmsd_list, capsize=2, color=conc_color_dict[conc],
                            zorder=zorder, alpha=all_alpha, ls='--')
            ax_o2c.errorbar(force_list, all_ar_list, yerr=all_rmsd_list, capsize=2, color=conc_color_dict[conc],
                            zorder=zorder, alpha=all_alpha, ls='--')
            ax_c2o.scatter(force_list, c2o_ar_list, s=35, color=conc_color_dict[conc], label=f'{conc}$\mu $M',
                        zorder=zorder)
            ax_o2c.scatter(force_list, o2c_ar_list, s=35, color=conc_color_dict[conc], label=f'{conc}$\mu $M',
                        zorder=zorder)
            ax_c2o.scatter(force_list, all_ar_list, s=35, color=conc_color_dict[conc],
                        zorder=zorder, alpha=all_alpha)
            ax_o2c.scatter(force_list, all_ar_list, s=35, color=conc_color_dict[conc],
                        zorder=zorder, alpha=all_alpha)

    picture_name = f"{project_name}_{special}_{special2}_dv{dv}_rm{rmsd}_pi{pale}"

    fig_c2on = outdir.strip("pdf") + picture_name + "_c2o" + ".png"
    pdf_c2on = outdir + "/" + picture_name + "_c2o" + ".pdf"
    ax_c2o.set_ylabel("Time(ms)")
    ax_c2o.set_ylim(ylim)
    ax_c2o.set_yticks(ytick)
    ax_c2o.set_xticks(xtick)
    ax_c2o.set_xlabel("Force(pN)")
    if dv == 14:
        ax_o2c.legend(loc='best')
    fig_c2o.savefig(fig_c2on, dpi=my_dpi)
    print(fig_c2on)

    fig_o2cn = outdir.strip("pdf") + picture_name + "_o2c" + ".png"
    pdf_o2cn = outdir + "/" + picture_name + "_o2c" + ".pdf"
    ax_o2c.set_ylabel("Time(ms)")
    ax_o2c.set_ylim(ylim)
    ax_o2c.set_yticks(ytick)
    ax_o2c.set_xticks(xtick)
    ax_o2c.set_xlabel("Force(pN)")
    fig_o2c.savefig(fig_o2cn, dpi=my_dpi)
    print(fig_o2cn)

    plt.close('all')