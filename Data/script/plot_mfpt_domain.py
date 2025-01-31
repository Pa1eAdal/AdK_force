# -*- codeing = utf-8 -*-
# Time : 2023/2/16 12:48
# File : plot_mfpt_domain.py
# Software : PyCharm

"""
the input file from calculate_mfpt_domain.py
"""

import itertools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import method
import os

plot_setting = method.Plot_setting()

project_name = "mfpt_akesi_force_o2c_real"
control = "c_force"
task_name = "domain_mfpt"

special = method.time_y_m_d()

dv = 0
conc_list = [300, 3000]
site = [148]
# lod_list = [70]
# dp = 0
# io = 0.05
# k_list = [0.01]
# rp_list = np.arange(0, 13, 2)
force_list = [0, 5, 7, 10, 15, 20]
n_seed = 400
rmsd = 0.8
pale = 0

frame = 10000

limit = frame + 2

if "o2c" in project_name:
    start_state = ["OO"]
    start_state2 = ["EE"]
    end_state = ["CC"]
    chi_cut = 0.9
elif "c2o" in project_name:
    start_state = ["CC"]
    start_state2 = ["DD"]
    end_state = ["OO"]
    chi_cut = 0.8
domain_list = ["nmp_time_phy", "lid_time_phy"]
# domain_list = ["oo"]
domain_label = [{300: "NMP 300$\mu$M", 3000: "NMP 3000$\mu$M"}, {300: "LID 300$\mu$M", 3000: "LID 3000$\mu$M"}]
label_dict = dict(zip(domain_list, domain_label))

my_dpi = 350
rect = [0.15, 0.135, 0.83, 0.845]
figsize = (2.25, 2)

position_list = ["LID(148)-CORE(22)", "NMP(42)-LID(144)", "NMP(55)-CORE(166)", "NMP(73)-LID(142)"]
position_dict = dict(zip(site, position_list))

outdir = f"./Figures/{task_name}/pdf"

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

dvn_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16]
p_close = ["0.28 (WT)", 0.52, 0.77, 0.15, 0.05, 0.02, 9e-3, 3e-3, 0.89, 0.94, 0.12, 0.10, 0.09, 0.986, 0.994, 0.008]
p_dict = dict(zip(dvn_list, p_close))

if "c_force" in control:
    # Create output directory if it doesn't exist
    out_data_dir = "./value"
    if not os.path.exists(out_data_dir):
        os.makedirs(out_data_dir, exist_ok=True)
        
    out_data_name = f'{out_data_dir}/MFPT_time.dat'
    with open(out_data_name, 'w') as data_file:
        color_list = ["tab:red", "tab:blue", 'tab:brown', 'tab:purple']
        color_dict = dict(zip(domain_list, color_list))
        fig = plt.figure(figsize=figsize, dpi=my_dpi)
        ax = fig.add_axes(rect)
        for conc in conc_list:
            c = f"{conc}".zfill(4)
            if conc == 300:
                alpha = 1
                fmt = "-"
            elif conc == 3000:
                alpha = 0.4
                fmt = "--"
            for ind, s, domain in itertools.product(range(len(start_state)), site, domain_list):
                zorder = domain_list.index(domain) + 5
                data_list = []
                rmsd_list = []

                if "phy" in domain:
                    for f in force_list:
                        fn = f"{task_name}/{task_name}_{project_name}_cut{chi_cut}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_{start_state[ind]}_{end_state[ind]}.npz"
                        data = np.load(fn)
                        oo_tot_time = data[f'{domain}']
                        
                        # Write header with metadata
                        data_file.write(f"{domain}\tdv{dv}\tforce{f}\tsite{s}\t[AMP]{conc}\n")
                        # Write raw data with seed numbers and conversion
                        for n in range(len(oo_tot_time)):
                            converted_time = oo_tot_time[n] * 112 / 2e5
                            data_file.write(f"n{n:03d}\t{converted_time}\n")
                            
                        oo_tot_time_old = method.sampling_r_limited(oo_tot_time, limit=limit, times=20, size=20)
                        oo_data_old = method.data_analysis(oo_tot_time_old)
                        data_list.append(oo_data_old[0])
                        rmsd_list.append(oo_data_old[1])
                else:
                    for f in force_list:
                        fn = f"ake_MFPT/ake_MFPT_mfpt_akesi_force_c2o_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_DDCC_OO.npz"
                        np.load.__defaults__ = (None, True, True, 'ASCII')
                        data = np.load(fn)
                        np.load.__defaults__ = (None, False, True, 'ASCII')
                        oo_tot_time = data[f'{domain}']
                        
                        # Write header with metadata
                        data_file.write(f"{domain}\tdv{dv}\tforce{f}\tsite{s}\t[AMP]{conc}\n")
                        # Write raw data with seed numbers and conversion
                        for n in range(len(oo_tot_time)):
                            converted_time = oo_tot_time[n] * 112 / 2e5
                            data_file.write(f"n{n:03d}\t{converted_time}\n")
                            
                        oo_tot_time_old = method.sampling_r_limited(oo_tot_time, limit=limit, times=20, size=20)
                        oo_data_old = method.data_analysis(oo_tot_time_old)
                        data_list.append(oo_data_old[0])
                        rmsd_list.append(oo_data_old[1])

                rmsd_list = np.asarray(rmsd_list)
                rmsd_list = rmsd_list * 112 / 2e5
                data_list = np.asarray(data_list)
                data_list = data_list * 112 / 2e5
                ax.errorbar(force_list, data_list, fmt=fmt, yerr=rmsd_list, capsize=3, color=color_dict[domain],
                            zorder=zorder, alpha=alpha)
                if "phy" in domain:
                    if conc == 300:
                        ax.plot(force_list, data_list, 's', ms=5, color=color_dict[domain],
                                label=f"{label_dict[domain][conc]}", zorder=zorder, alpha=alpha)
                    elif conc == 3000:
                        ax.plot(force_list, data_list, 's', ms=5, color=color_dict[domain],
                                label=f"{label_dict[domain][conc]}",
                                zorder=zorder, markerfacecolor="none", alpha=alpha)
                else:
                    ax.plot(force_list, data_list, 'o', ms=5, color=color_dict[domain], label=f"{label_dict[domain]}",
                            zorder=zorder, markerfacecolor="none")

        picture_name = f"{task_name}_{project_name}_{special}_{control}_dv{dv}_rm{rmsd}_pi{pale}_{start_state[ind]}_{end_state[ind]}"
        fign = outdir.strip("pdf") + picture_name + ".png"
        pdfn = outdir + "/" + picture_name + ".pdf"

        ax.set_ylabel("Time(ms)")
        ax.set_xlabel("Force(pN)")
        ax.set_ylim([0, 5])
        legend_elements = [Line2D([0], [0], color='tab:red', marker='s', ms=5, label='NMP'),
                           Line2D([0], [0], color='tab:blue', marker='s', ms=5, label='LID')]
        if "c2o" in project_name:
            pass
        else:
            ax.legend(handles=legend_elements, loc='best')
        fig.savefig(fign, dpi=my_dpi)
        print(fign)