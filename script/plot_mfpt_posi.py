# -*- codeing = utf-8 -*-
# Time : 2022/10/6 19:28
# File : plot_mfpt_posi.py
# Software : PyCharm
"""
the input file from calculate_ake_MFPT.py
"""

import itertools
import numpy as np
from matplotlib import pyplot as plt
import method
import os

special = ['o2c', 'c2o']
control = "paper"
task_name = "ake_MFPT"

dv = 0
conc_list = [300, 3000]
site = [148]
# lod_list = [70]
# dp = 0
# io = 0.05
# k_list = [0.01]
# rp_list = np.arange(0, 13, 2)
force_list = [0, 5, 7, 10, 15, 20]
n_seed_s = 20
rmsd = 0.8
pale = 0

# domain_list = ["nmp", "lid", "oo"]
domain_list = ["oo"]

all_alpha = 0.4
my_dpi = 350
rect = [0.15, 0.135, 0.83, 0.845]
figsize = (2.25, 2)

plot_setting = method.Plot_setting()

position_list = ["LID(148)-CORE(22)", "NMP(42)-LID(144)", "NMP(55)-CORE(166)", "NMP(73)-LID(142)"]
position_dict = dict(zip(site, position_list))

outdir = "./Figures/mfpt_mcc/pdf"
tr_dir = "./txtdata/"

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

dvn_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16]
p_close = [0.28, 0.52, 0.77, 0.15, 0.05, 0.02, 9e-3, 3e-3, 0.89, 0.94, 0.12, 0.10, 0.09, 0.986, 0.994, 0.008]
p_dict = dict(zip(dvn_list, p_close))

if "paper" in control:
    # Create output directory if it doesn't exist
    out_data_dir = "./value"
    if not os.path.exists(out_data_dir):
        os.makedirs(out_data_dir, exist_ok=True)
        
    color_list = ["tab:red", "tab:blue", 'tab:blue', 'm', 'yellow']
    color_dict = dict(zip(conc_list, color_list))
    for domain, spec in itertools.product(domain_list, special):
        out_data_name = f'{out_data_dir}/MFPT_all_{spec}.dat'
        with open(out_data_name, 'w') as data_file:

            if spec == "o2c":
                project_name = "mfpt_akesi_force_o2c_long"
            elif spec == "c2o":
                project_name = "mfpt_akesi_force_c2o_real"
            if "o2c" in project_name:
                start_state = ["EEOO"]
                start_state2 = ["EE"]
                end_state = ["CC"]
                frame = 10000
            elif "c2o" in project_name:
                start_state = ["DDCC"]
                start_state2 = ["DD"]
                end_state = ["OO"]
                frame = 40000

            limit = frame + 2

            fig = plt.figure(figsize=figsize, dpi=my_dpi)
            ax = fig.add_axes(rect)
            for ind, s, conc in itertools.product(range(len(start_state)), site, conc_list):
                zorder = conc
                c = f"{conc}".zfill(4)
                data_list = []
                rmsd_list = []
                all_ar_list = []
                all_rmsd_list = []
                x_list_all = []
                for f in force_list:
                    fn = f"{task_name}/{task_name}_{project_name}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_{start_state[ind]}_{end_state[ind]}.npz"
                    data = np.load(fn)
                    interval_list = []
                    print(fn)
                    for n_s in range(n_seed_s):
                        n = f"{n_s:03d}"
                        tr_filename = f"akesi_force_dv{dv}_rm{rmsd}_pi{pale}_c{c}_s{s}_f{f}_n{n}.npz"
                        tr_fn = tr_dir + tr_filename
                        tr_data = np.load(tr_fn)
                        nframe_list = tr_data['nframe']
                        tem_list = nframe_list[::2]
                        for i in range(len(tem_list) - 1):
                            interval_list.append(tem_list[i + 1] - tem_list[i])
                    all_list = np.asarray(interval_list)
                    all_list = all_list * 112 / 2e8
                    all_tot_time = method.sampling_r(all_list, times=50, size=50)
                    all_data = method.data_analysis(all_tot_time)
                    all_ar_list.append(all_data[0])
                    all_rmsd_list.append(all_data[1])

                    if (f == 21) & ("o2c" in project_name):
                        print(spec)
                        print(project_name)
                        oo_data_old = data["oo_data_old"]
                    else:
                        if (f == 20) and "o2c" in project_name and "real" in project_name:
                            print(spec)
                            print(project_name)
                            oo_data_old = data["oo_data_old"]
                            oo_tot_time = data['oo_data']
                        else:
                            oo_tot_time = data['all_time']
                            
                        # Write header with metadata
                        data_file.write(f"dv{dv}\tforce{f}\tsite{s}\t[AMP]{conc}\n")
                        # Write raw data with seed numbers and conversion
                        for n in range(len(oo_tot_time)):
                            converted_time = oo_tot_time[n] * 112 / 2e5
                            data_file.write(f"n{n:03d}\t{converted_time}\n")
                            
                        oo_tot_time_old = method.sampling_r_limited(oo_tot_time, limit=limit)
                        oo_tot_time = method.sampling_r(oo_tot_time)

                        oo_data_old = method.data_analysis(oo_tot_time_old)
                        oo_data = method.data_analysis(oo_tot_time)

                    data_list.append(oo_data[0])
                    rmsd_list.append(oo_data[1])

                rmsd_list = np.asarray(rmsd_list)
                rmsd_list = rmsd_list * 112 / 2e5
                data_list = np.asarray(data_list)
                data_list = data_list * 112 / 2e5

                ax.errorbar(force_list, data_list, yerr=rmsd_list, capsize=2, color=color_dict[conc], zorder=zorder)
                ax.errorbar(force_list, all_ar_list, yerr=all_rmsd_list, capsize=2, color=color_dict[conc],
                            zorder=zorder, alpha=all_alpha, ls='--')
                ax.scatter(force_list, data_list, s=35, color=color_dict[conc], label=f"{conc}$\mu$M", zorder=zorder)
                ax.scatter(force_list, all_ar_list, s=35, facecolor="none", edgecolor=color_dict[conc], zorder=zorder, alpha=all_alpha)

                picture_name = f"{task_name}_{project_name}_{method.time_y_m_d()}_{control}_dv{dv}_rm{rmsd}_pi{pale}_c{c}_{domain}_{start_state[ind]}_{end_state[ind]}"
                fign = outdir.strip("pdf") + picture_name + ".png"
                pdfn = outdir + "/" + picture_name + ".pdf"

            ax.set_ylabel("Time(ms)")
            ax.set_xlabel("Force(pN)")
            ax.set_ylim([0, 9.5])
            ax.set_yticks(np.arange(0, 10, 3))
            ax.legend(loc="upper left")
            fig.savefig(fign, dpi=my_dpi)
            print(fign)
            plt.close('all')