"""
Created on Wed Mar  3 20:43:40 2021

@author: yyzhang

Modified on Fri Aug  5 23:50 2022

@developer: zzy

need no chemical reaction
"""
"""
Input: .xtc file from dcd2xtc (from mdtraj)
"""
import numpy as np
import mdtraj as md
import itertools
import os
import enzyme_setting

# i_bed_lid=122
# i_stop_lid=160

# i_beg_nmp=30
# i_stop_nmp=75

# i_beg_core1=1
# i_stop_core1=29

# i_beg_core2=76
# i_stop_core2=121

# i_beg_core3=161
# i_stop_core3=214


def get_rlc(t):
    core = md.compute_center_of_mass(t.atom_slice(t.top.select('chainid 0 2 4')))
    lid = md.compute_center_of_mass(t.atom_slice(t.top.select('chainid 3')))
    return np.linalg.norm(lid - core, axis=1)


def get_rnc(t):
    core = md.compute_center_of_mass(t.atom_slice(t.top.select('chainid 0 2 4')))
    nmp = md.compute_center_of_mass(t.atom_slice(t.top.select('chainid 1')))
    return np.linalg.norm(nmp - core, axis=1)


task_name = "akesi"
adk = enzyme_setting.Enzyme_adk(cut='chi', dv_list=[0])
adk.amp_list = [300]

topn = os.path.expanduser(f"~/data/pdb/4ake_5unit0_CA.pdb")
data_dir = os.path.expanduser(f"~/mnt/portable_storage")

outdir = 'pclose'

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

# for co in conclist:
if "dna" in task_name:

    lod_list = [70]
    dp_list = [0]
    io_list = [0.05]
    for dv, amp, s, dp, lod, io in itertools.product(dv_list, amp_list, site_list, dp_list, lod_list, io_list):
        topn = (os.path.expanduser(f"~/data/top.pdb"))
        m = f"{amp:04d}"
        rlc = np.empty((0,))
        rnc = np.empty((0,))
        for ns in np.arange(20):
            n = f"{ns}".zfill(3)
            # trjn = dirn + 'dv%d/ake_dv%d_n%03d.xtc' % (i, i, n)
            # trjn = f"{data_dir}/{task_name}/dv{dv}/{cut_val_c}{cut_val}/pi{pale}/c{m}/s{s}/f{f}/{task_name}_dv{dv}_{cut_val_c}{cut_val}_pi{pale}_c{m}_s{s}_f{f}_n{n}.xtc"
            trjn = (f"{task_name}/dv{dv}/{cut}{cut_val}/pi{pale}/c{m}/s{s}/dp{dp}/l{lod}/io{io}/"
                    f"{task_name}_dv{dv}_{cut}{cut_val}_pi{pale}_c{m}_s{s}_dp{dp}_l{lod}_io{io}_n{n}.xtc")

            data = md.load(trjn, top=topn)
            print(trjn)
            rlc = np.hstack((rlc, get_rlc(data)))
            rnc = np.hstack((rnc, get_rnc(data)))
        print(dv, pale, m, dp, lod, io, rlc.shape, rnc.shape)
        # np.savez('pclose/rlc-rnc_dv%d' % (i), rlc=rlc, rnc=rnc)
        np.savez(f"pclose/rlc-rnc_{task_name}_dv{dv}_{cut}{cut_val}_pi{pale}_c{m}_s{s}_dp{dp}_l{lod}_io{io}", rlc=rlc, rnc=rnc)
else:
    for dv, ba, cf, cb, atp, amp, adp in itertools.product(adk.dv_list, adk.ba_list, adk.cf_list, adk.cb_list,
                                                           adk.atp_list, adk.amp_list, adk.adp_list):
        bas = ba[0]
        bap = ba[1]
        m = f"{amp:04d}"
        t = f"{atp:04d}"
        d = f"{adp:04d}"
        rlc = np.empty((0,))
        rnc = np.empty((0,))
        xtc_dir = f"{task_name}/dv{dv}/{adk.cut}{adk.cut_val}/p1_{adk.p1}/p7_{adk.p7}/bas{bas}/bap{bap}/cf{cf}/cb{cb}/t{t}/m{m}/d{d}"
        prefix = f"{task_name}_dv{dv}_{adk.cut}{adk.cut_val}_p1_{adk.p1}_p7_{adk.p7}_bas{bas}_bap{bap}_cf{cf}_cb{cb}_t{t}_m{m}_d{d}"
        for ns in np.arange(20):
            n = f"{ns}".zfill(3)
            trjn = f"{xtc_dir}/{prefix}_n{n}.xtc"
            data = md.load(trjn, top=topn)
            print(trjn)
            rlc = np.hstack((rlc, get_rlc(data)))
            rnc = np.hstack((rnc, get_rnc(data)))
        print(dv, adk.p1, adk.p7, t, m, d, rlc.shape, rnc.shape)
        # np.savez('pclose/rlc-rnc_dv%d' % (i), rlc=rlc, rnc=rnc)
        np.savez(f"{outdir}/rlc-rnc_{prefix}", rlc=rlc, rnc=rnc)
