#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 22:44:46 2024
The code creates network plots for the case of circular boundary.
@author: abhinav
"""

import numpy as np
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.pyplot as plt
import math
#import networkx as nx
from matplotlib.collections import LineCollection
import matplotlib.ticker as ticker

# angle of a vector with horizontal
# -----------------------------------------------------------------------------
# Small vectorized helpers for GitHub-ready analysis code.
# These replace repeated Python loops without changing the simulation outputs.
# -----------------------------------------------------------------------------
def angle(x1, y1, x2, y2):
    """Angle of a vector with respect to the horizontal axis.

    Works with either scalars or NumPy arrays.
    """
    return np.arctan2(y2 - y1, x2 - x1)


def flatten_segments(x_start, y_start, x_end, y_end):
    """Return flattened x/y arrays and Nx2x2 segments for LineCollection."""
    segments = np.stack(
        (
            np.column_stack((x_start, y_start)),
            np.column_stack((x_end, y_end)),
        ),
        axis=1,
    )
    x_flat = np.column_stack((x_start, x_end)).ravel()
    y_flat = np.column_stack((y_start, y_end)).ravel()
    return x_flat, y_flat, segments


def build_bond_plot_arrays(xpos, ypos, conn_node, outer_nodes, L, unconnected_const, pbd_flag):
    """Build bond coordinate arrays used by LineCollection.

    The old code repeatedly used np.hstack inside nested loops. This version first
    finds all valid bonds once, then slices positions in vectorized form for each
    saved force step.
    """
    n_steps, num_pts = xpos.shape
    nodes = np.arange(num_pts)
    rows = nodes // L
    cols = nodes % L
    boundary_nodes = np.unique(nodes[(rows == 0) | (rows == L - 1) | (cols == 0) | (cols == L - 1)])

    valid = conn_node[:, :3] != unconnected_const
    src, bond = np.where(valid)
    dst = conn_node[src, bond]

    outer_mask_src = np.isin(src, outer_nodes)
    outer_mask_dst = np.isin(dst, outer_nodes)
    plot_mask = ~(outer_mask_src | outer_mask_dst)

    x_loc_arr, y_loc_arr, x_loc_plot, y_loc_plot = [], [], [], []
    for step in range(n_steps):
        x0 = xpos[step, src]
        y0 = ypos[step, src]
        x1 = xpos[step, dst].copy()
        y1 = ypos[step, dst]

        if pbd_flag == 1:
            left = cols[src] == 0
            right = cols[src] == L - 1
            x1 = np.where(left & ((x1 - x0) > L / 2), x1 - L, x1)
            x1 = np.where(right & ((x0 - x1) > L / 2), x1 + L, x1)

        x_all, y_all, _ = flatten_segments(x0, y0, x1, y1)
        x_plot, y_plot, _ = flatten_segments(x0[plot_mask], y0[plot_mask], x1[plot_mask], y1[plot_mask])

        x_loc_arr.append(x_all)
        y_loc_arr.append(y_all)
        x_loc_plot.append(x_plot)
        y_loc_plot.append(y_plot)

    node_bond_single = np.column_stack((src[plot_mask], bond[plot_mask])).astype(int)
    node_bond_plot = np.broadcast_to(node_bond_single, (n_steps, *node_bond_single.shape)).copy()
    flat_plot_index = node_bond_single[:, 0] * 3 + node_bond_single[:, 1]

    return x_loc_arr, y_loc_arr, x_loc_plot, y_loc_plot, boundary_nodes, node_bond_plot, flat_plot_index


def smooth_circular_vectors(values):
    """Three-point circular moving average used for boundary-force quivers."""
    return (np.roll(values, 1) + values + np.roll(values, -1)) / 3.0


#==============================================================================
# =============================================================================
# USER SETTINGS
# =============================================================================

only_last_flag = 0
only_last_num = 2

L = 64
num_pts = L * L

srand_flag = 1
srand = 112 if srand_flag else None

num_center = 5
num_dip = num_center * 6

num = 10 + 1  # total number of force steps

pbond = 0.55

mu = 1
mu_c = 1
tol = 1e-7
kappa = 1e-5
force_val = 0.01
strain_thresh = 1e-3 if pbond == 1 else 1e-6

rlen = 0.9
rlen_txt = f"{rlen:.4f}"

base = "/home/abhinav/david/"

# =============================================================================
# FLAGS
# =============================================================================

rlen_flag = 0
rlen_hex_flag = 0
rlen_hex_sep_flag = 0
sep_flag = 0
pbd_flag = 0
buckle_flag = 0
radial_flag = 0
inner_radial_flag = 0
inner_radial_hex_arp_flag = 0
inner_radial_arp_flag_bash = 0
radial_only_flag = 0
cluster_flag = 0
cluster_new_flag = 0
cluster_radial_flag = 1
hex_flag = 0
hex_rand_flag = 0
test_hex = 0
test_hex_rand = 0
anisotropic_flag = 0
force_flag = 0
diff_flag = 1
auto_flag = 0

if anisotropic_flag:
    num_dip = num_center

if force_flag:
    force_val = 0.01115
    force_txt = f"{force_val:.9f}"

if buckle_flag:
    mu_c = 0.1

# =============================================================================
# NETWORK SEEDS / COMMAND-LINE OVERRIDES
# =============================================================================

if diff_flag:
    diff_network1 = 111111
    diff_network2 = 101111
    diff_network = f"{diff_network1},{diff_network2}"

if auto_flag:
    import sys

    if pbond < 1:
        print(sys.argv[1], " ** ", sys.argv[2], " ** ", sys.argv[3])
        diff_network = sys.argv[1]
        num_center = int(sys.argv[2])
        if srand_flag:
            srand = int(sys.argv[3])
        print("Diff network:", diff_network)

    else:
        print(sys.argv[1], " ** ", sys.argv[2])
        num_center = int(sys.argv[1])
        if srand_flag:
            srand = int(sys.argv[2])

    num_dip = num_center * 6

# =============================================================================
# STRINGS USED IN FILE NAMES
# =============================================================================

pbond_string = f"{pbond:.2f}"
tol_str = f"{tol:.2e}"
kappa_str = "0.00e+00" if kappa == 0 else f"{kappa:.2e}"
force_str = f"{force_val:.4f}"
mu_str = f"{mu:.4f}"
mu_c_str = f"{mu_c:.4f}"

# =============================================================================
# FOLDER SELECTION
# =============================================================================

def kappa_folder(kappa_value):
    """Return the folder name used for a given bending stiffness."""
    mapping = {
        0: "kappa2_0/",
        2e-7: "kappa2_2e-7/",
        5e-7: "kappa2_5e-7/",
        1e-6: "kappa2_e-6/",
        2e-6: "kappa2_2e-6/",
        5e-6: "kappa2_5e-6/",
        1e-5: "kappa2_e-5/",
        2e-5: "kappa2_2e-5/",
        5e-5: "kappa2_5e-5/",
        1e-4: "kappa2_e-4/",
        2e-4: "kappa2_2e-4/",
        5e-4: "kappa2_5e-4/",
        1e-3: "kappa2_e-3/",
        2e-3: "kappa2_2e-3/",
        5e-3: "kappa2_5e-3/",
        1e-2: "kappa2_e-2/",
    }
    return mapping[kappa_value]


def active_folder_key():
    """Return the active simulation-folder key based on the enabled flags."""
    flag_to_folder = [
        (radial_flag, "lattice_nopbd_bndry_all_clamp_restlength_circle_radial/"),
        (inner_radial_flag, "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial/"),
        (inner_radial_hex_arp_flag, "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/"),
        (radial_only_flag, "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/"),
        (inner_radial_arp_flag_bash, "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/"),
        (cluster_flag, "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/"),
        (cluster_new_flag, "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/"),
        (cluster_radial_flag, "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/"),
        (hex_flag, "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/"),
        (hex_rand_flag, "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/"),
        (force_flag, "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_force/"),
        (test_hex, "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/"),
        (test_hex_rand, "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/"),
        (anisotropic_flag, "lattice_nopbd_bndry_all_clamp_restlength_circle_arp_anisotropic/"),
    ]

    for is_active, folder_name in flag_to_folder:
        if is_active:
            return folder_name

    if diff_flag and srand_flag:
        return "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed/"

    return "lattice_nopbd_bndry_all_clamp_restlength_circle/"


folder = "lattice_nopbd_bndry_all_clamp_restlength_circle/"

if pbond == 1:
    folder = active_folder_key()

    if L == 64 and srand == 300:
        folder = "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_large/"

else:
    ranseed = "667720,601210"
    kappa_fname = kappa_folder(kappa)

    folder_root = active_folder_key()

    if srand_flag and diff_flag:
        folder = f"{folder_root}{kappa_fname}{ranseed}/{srand}/{diff_network}/"
    elif srand_flag:
        folder = f"{folder_root}{kappa_fname}{ranseed}/{srand}/"
    else:
        folder = f"{folder_root}{kappa_fname}{ranseed}/"
        
pos_folder = base + folder + "txt/strain/Lattice_" 
strain_pos_folder = base + folder + "txt/strain/Strain_Lattice_" 
rlen_folder = base + folder + "txt/rlen/rlen_" 
node_en_folder = base + folder + "energy/node/Lattice_node_" 

# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1:
    pos_folder = base + folder + "txt/strain/Lattice_" + "srand_" + str(srand) + "_"
    strain_pos_folder = base + folder + "txt/strain/Strain_Lattice_"+ "srand_"  + str(srand) + "_"
    rlen_folder = base + folder + "txt/rlen/rlen_"+ "srand_"  + str(srand) + "_" 
    node_en_folder = base + folder + "energy/node/Lattice_node_srand_" + str(srand) + "_" 

if L != 64 and srand_flag == 1:
    pos_folder = base + folder + "txt/strain/Lattice_"+str(L)+'_' 
    strain_pos_folder = base + folder + "txt/strain/Strain_Lattice_"+str(L)+'_' 
    rlen_folder = base + folder + "txt/rlen/rlen_"+str(L)+'_' 
    
if force_flag == 1:
    pos_folder = base + folder + "txt/strain/Lattice_" + "srand_" + str(srand) + "_"
    strain_pos_folder = base + folder + "txt/strain/Strain_Lattice_"+ "srand_"  + str(srand) + "_"
    rlen_folder = base + folder + "txt/rlen/rlen_"+ "srand_"  + str(srand) + "_" 
    
if L == 64 and pbond == 1 and srand == 300:
    pos_folder = base + folder + "txt/strain/Lattice_" + "srand_" + str(srand) + "_"
    strain_pos_folder = base + folder + "txt/strain/Strain_Lattice_"+ "srand_"  + str(srand) + "_"
    rlen_folder = base + folder + "txt/rlen/rlen_"+ "srand_"  + str(srand) + "_" 
    node_en_folder = base + folder + "energy/node/Lattice_node_srand_" + str(srand) + "_" 

    
all_data = []
strain_all_data = []
rlen_all_data = []

if only_last_flag == 0:
    for i in range(0,num):     # for all positions
        # file for restlength changes
        fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"        # input file name
        #fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+"_force.txt"        # input file name
        strain_fname = strain_pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"  # strain file
        #strain_fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+"_force.txt"  # strain file
                    
        rlen_fname = rlen_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"  # rlen file
        
    
        if force_flag == 1:
            fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"        # input file name
            #fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+"_force.txt"        # input file name
            strain_fname = strain_pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"  # strain file
            #strain_fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+"_force.txt"  # strain file
                        
            rlen_fname = rlen_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"  # rlen file
            
   
        rlen_data = np.loadtxt(rlen_fname)
        rlen_all_data.append(rlen_data)
    
    
        print("Node position file name: ",fname)
        print("strain file name: ", strain_fname)
        print("rlen file name: ", rlen_fname)
    
    
        data = np.loadtxt(fname)
    #    xpos = data[:,1]
    #    ypos = data[:,2]
        all_data.append(data)
        
        strain_data = np.loadtxt(strain_fname)
        strain_all_data.append(strain_data)

else:
    for i in range(0,only_last_num):     # for all positions
        # file for restlength changes
        fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"        # input file name
        strain_fname = strain_pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"  # strain file                    
        rlen_fname = rlen_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"  # rlen file

        if i == only_last_num-1 and (srand != 123 or srand != 124):
            fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+"_force.txt"        # input file name
            rlen_fname = rlen_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+".txt"  # rlen file
             
        rlen_data = np.loadtxt(rlen_fname)
        rlen_all_data.append(rlen_data)
        print("* Node position file name: ",fname)
        print("* strain file name: ", strain_fname)
        print("* rlen file name: ", rlen_fname)
    
        data = np.loadtxt(fname)
    #    xpos = data[:,1]
    #    ypos = data[:,2]
        all_data.append(data)
        
        strain_data = np.loadtxt(strain_fname)
        strain_all_data.append(strain_data)
    

#all_data = np.array(all_data)

all_data = np.array(all_data)
xpos = all_data[:,:,1]
ypos = all_data[:,:,2]

rlen_all_data = np.array(rlen_all_data)
rlen_all_data = rlen_all_data[:,:,:3]

if only_last_flag == 0:
    rlen_data_val = np.reshape(rlen_all_data,(num,3*L*L))
else:
    rlen_data_val = np.reshape(rlen_all_data,(only_last_num,3*L*L))

###############################################################################
# reading the connection array
fname_connect = base+folder+'txt/strain/Lattice_connect_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
#fname_connect = base+folder+'txt/strain/34Lattice_connect_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond != 1:
    fname_connect = base+folder+'txt/strain/Lattice_connect_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond == 1:
    fname_connect = base+folder+'txt/strain/Lattice_connect_'+ "srand_" + str(srand) + '_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
if srand_flag == 1 and pbond == 1 and L ==64:
    fname_connect = base+folder+'txt/strain/Lattice_connect_'+ "srand_" + str(srand) + '_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
if force_flag == 1:
    fname_connect = base+folder+'txt/strain/Lattice_connect_'+ "srand_" + str(srand) + '_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
        
        
print("connection filename: ",fname_connect)
conn_data = np.loadtxt(fname_connect)
conn_node = conn_data[:,3:9]
conn_node = (np.rint(conn_node)).astype(int)


###############################################################################
# reading the bending energy at each node
node_en_fname = node_en_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+"_force.txt"  # rlen file
if only_last_flag == 1:
    node_en_fname = node_en_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(only_last_num-1)+"_"+str(num-1)+"_force.txt"  # rlen file

print("energy per node filename: ",node_en_fname)
en_data = np.loadtxt(node_en_fname)
bend_en_node_all = en_data[:,7:]
bend_en_node = np.sum(bend_en_node_all, axis = 1)

###############################################################################    
# setting the circular boundary
###############################################################################
# making a circular region to identify the boundary nodes so that I can fix them in the simulation
radius = 25
dr = 0.5
center = 2080

if L == 128:
    center = 8383
    radius = 50

# reading the nodes on the circular boundary
fname_boundary = base+folder+'txt/area/Boundary_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond != 1:
    fname_boundary = base+folder+'txt/area/Boundary_nodes_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond == 1:
    fname_boundary = base+folder+'txt/area/Boundary_nodes_'+ "srand_" + str(srand) + '_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64:
    fname_boundary = base+folder+'txt/area/Boundary_nodes_'+ "srand_" + str(srand) + '_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if force_flag == 1:
    fname_boundary = base+folder+'txt/area/Boundary_nodes_'+ "srand_" + str(srand) + '_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
    
        
print("circular boundary filename: ",fname_boundary)
boundary_data = np.loadtxt(fname_boundary)
boundary_nodes_circle = boundary_data[:,1]
boundary_nodes_circle = boundary_nodes_circle.astype(int)

# manually checking the nodes that should not be send to the energy minimizer
# outer_nodes_circle = []
# for i in range(0,num_pts):
#     dx = xpos[0][center] - xpos[0][i]
#     dy = ypos[0][center] - ypos[0][i]
#     dist_center = np.sqrt(dx*dx + dy*dy)
#     if (dist_center > radius+dr):
#         outer_nodes_circle.append(i)

# reading the nodes outside the inner region :: this includes the circular boundary nodes btw
fname_outer = base+folder+'txt/area/Outer_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond != 1:
    fname_outer = base+folder+'txt/area/Outer_nodes_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond == 1:
    fname_outer = base+folder+'txt/area/Outer_nodes_'+ "srand_" + str(srand) +'_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64:
    fname_outer = base+folder+'txt/area/Outer_nodes_'+ "srand_" + str(srand) +'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if force_flag == 1:    
    fname_outer = base+folder+'txt/area/Outer_nodes_'+ "srand_" + str(srand) +'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
        
print("outer nodes filename: ",fname_outer)
outer_data = np.loadtxt(fname_outer)
#outer_nodes = outer_data
outer_nodes = outer_data.astype(int)

# reading the nodes in the inner region :: 
fname_inner = base+folder+'txt/area/Inner_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond != 1:
    fname_inner = base+folder+'txt/area/Inner_nodes_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond == 1:
    fname_inner = base+folder+'txt/area/Inner_nodes_'+ "srand_" + str(srand) +'_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64:
    fname_inner = base+folder+'txt/area/Inner_nodes_'+ "srand_" + str(srand) +'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if force_flag == 1:
    fname_inner = base+folder+'txt/area/Inner_nodes_'+ "srand_" + str(srand) +'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
    
        
if srand_flag == 1 and diff_flag == 1:
    print("inner nodes filename: ",fname_inner)
    inner_data = np.loadtxt(fname_inner)
    #outer_nodes = outer_data
    inner_nodes = inner_data.astype(int)
    inner_nodes = inner_nodes[:,1]

# Remove circular-boundary nodes from the outer-node list using vectorized set logic.
outer_nodes = np.setdiff1d(outer_nodes, boundary_nodes_circle, assume_unique=False)

#==============================================================================
# setting the node and line locations for plotting
#==============================================================================

xpos = np.asarray(xpos)
ypos = np.asarray(ypos)

unconnected_const = -1 if L == 128 else 9999

(
    x_loc_arr,
    y_loc_arr,
    x_loc_plot,
    y_loc_plot,
    boundary_nodes,
    node_bond_plot,
    flat_plot_index,
) = build_bond_plot_arrays(xpos, ypos, conn_node, outer_nodes, L, unconnected_const, pbd_flag)

# finding the list of inner and ring nodes including the boundary nodes
all_nodes = np.arange(num_pts)
circle_nodes = np.setdiff1d(all_nodes, outer_nodes, assume_unique=False)   # nodes inside/on the outer boundary
# ###############################################################################    
# # setting up the strain array
# ###############################################################################
# reading the strain data

strain_all_data = np.array(strain_all_data)
bond_len = strain_all_data[:,:,1:4]

if only_last_flag == 0:
    bond_len = np.reshape(bond_len,(num,3*L*L))
else:
    bond_len = np.reshape(bond_len,(only_last_num,3*L*L))
#bond_len = np.reshape(bond_len,(num-6,3*L*L))

rlen_circle_data = rlen_all_data[:,circle_nodes,:3]

if only_last_flag == 0:
    rlen_data_val_circle = np.reshape(rlen_circle_data,(num,3*len(circle_nodes)))
else:
    rlen_data_val_circle = np.reshape(rlen_circle_data,(only_last_num,3*len(circle_nodes)))

bond_len_circle = strain_all_data[:,circle_nodes,1:4]     # bond length of bonds associated with nodes in the outer circle, including the boundary
# bond_len_circle = np.reshape(bond_len_circle,(num,3*len(circle_nodes)))

# Vectorized strain calculation.  flat_plot_index is aligned with node_bond_plot,
# x_loc_plot and y_loc_plot, so later plotting/stress code keeps the same ordering.
strain_circle_node = node_bond_plot[0].copy()
strain_circle = bond_len[:, flat_plot_index] - rlen_data_val[:, flat_plot_index]

# Full strain array, keeping zero-length/nonexistent bonds at zero.
strain_all = np.zeros_like(bond_len)
valid_bond_lengths = bond_len != 0
strain_all[valid_bond_lengths] = (bond_len - rlen_data_val)[valid_bond_lengths]
strain = [row[mask] for row, mask in zip(strain_all, valid_bond_lengths)]
            

strain = np.asarray(strain)
strain_circle = np.asarray(strain_circle)

# just for now a hack:: REMOVE LATER !!!!!
strain[0][:] = 0.0
strain_circle[0][:] = 0.0

###############################################################################
# plotting strain histograms
hist_base = base + folder + "png/strain_hist/"

###############################################################################    
###############################################################################    
# finding the stress on the inner boundary
# this is done for bonds that cross a circular region defined as the boundary
# of the inner reigon
###############################################################################    
###############################################################################

# setting the location of the inner radius
inner_radius = 12      # plus one because 
if L == 128:
    inner_radius = 24
    inner_radius = 12
theta_list = np.arange(0,2*np.pi,0.01)
# x_inner_circle = np.zeros((len(xpos)))
# y_inner_circle = np.zeros((len(ypos)))
x_inner_circle = []
y_inner_circle = []
inner_radius_for_plot = 12.5
x_inner_circle_13 = []
y_inner_circle_13 = []
inner_radius_for_plot_13 = 13
for i in range(0,len(x_loc_arr)):
    # x_inner_circle.append(xpos[i][center] + inner_radius*np.cos(theta_list))
    # y_inner_circle.append(ypos[i][center] + inner_radius*np.sin(theta_list))
    x_inner_circle.append(xpos[i][center] + inner_radius_for_plot*np.cos(theta_list))
    y_inner_circle.append(ypos[i][center] + inner_radius_for_plot*np.sin(theta_list))
    x_inner_circle_13.append(xpos[i][center] + inner_radius_for_plot_13*np.cos(theta_list))
    y_inner_circle_13.append(ypos[i][center] + inner_radius_for_plot_13*np.sin(theta_list))

# locating the bonds that cross the inner radius
strain_bond_inner_circle = []
node_cross = []
bond_cross = []
bond_cross_index = []
angle_cross_list = []
strain_bond_cross_x_list = []
strain_bond_cross_y_list = []
bond_mid_x_list = []
bond_mid_y_list = []
radial_force_cross_bond_list = []

for k in range(node_bond_plot.shape[0]):      # keep one readable loop over saved steps
    nodes = node_bond_plot[k, :, 0]
    bonds = node_bond_plot[k, :, 1]
    connected_nodes = conn_node[nodes, bonds]

    dx_node = xpos[k, nodes] - xpos[k, center]
    dy_node = ypos[k, nodes] - ypos[k, center]
    dr_node = np.hypot(dx_node, dy_node)

    dx_connected = xpos[k, connected_nodes] - xpos[k, center]
    dy_connected = ypos[k, connected_nodes] - ypos[k, center]
    dr_connected = np.hypot(dx_connected, dy_connected)

    cross_mask = ((dr_node < inner_radius) & (dr_connected > inner_radius)) | ((dr_node > inner_radius) & (dr_connected < inner_radius))
    cross_idx = np.flatnonzero(cross_mask)

    # multiplying by -1 because stretched bonds give negative stress and vice versa as per the original convention
    strain_circle[k, cross_idx] *= -1
    strain_cross = strain_circle[k, cross_idx]
    nodes_cross = nodes[cross_idx]
    bonds_cross = bonds[cross_idx]
    connected_cross = connected_nodes[cross_idx]

    angle_cross = angle(xpos[k, nodes_cross], ypos[k, nodes_cross], xpos[k, connected_cross], ypos[k, connected_cross])
    strain_x = strain_cross * np.cos(angle_cross)
    strain_y = strain_cross * np.sin(angle_cross)

    bond_mid_x = 0.5 * (xpos[k, nodes_cross] + xpos[k, connected_cross])
    bond_mid_y = 0.5 * (ypos[k, nodes_cross] + ypos[k, connected_cross])
    radial_norm = np.hypot(bond_mid_x, bond_mid_y)
    radial_force = np.divide(strain_x * bond_mid_x + strain_y * bond_mid_y, radial_norm, out=np.zeros_like(strain_cross), where=radial_norm != 0)

    strain_bond_inner_circle.append(strain_cross)
    node_cross.append(nodes_cross)
    bond_cross.append(bonds_cross)
    bond_cross_index.append(cross_idx)
    angle_cross_list.append(angle_cross)
    strain_bond_cross_x_list.append(strain_x)
    strain_bond_cross_y_list.append(strain_y)
    bond_mid_x_list.append(bond_mid_x)
    bond_mid_y_list.append(bond_mid_y)
    radial_force_cross_bond_list.append(radial_force)

total_radial_stress = [np.sum(values) / (2 * np.pi * inner_radius) for values in radial_force_cross_bond_list]
# writing stress results to file for later analysis

# writing all the bond data : for bonds that cross the inner boundary
stress_bond_outfname = base+folder+"txt/force/"+"inner_bndry_bond_radial_force_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if L != 64:
    stress_bond_outfname = base+folder+"txt/force/"+"inner_bndry_bond_radial_force_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    
print('here all individual radial forces on all bonds that cross the inner boundary are written to file : ',stress_bond_outfname)    
heading = '    node          bond       bond angle        pos_x         pos_y       strain_x      strain_y        strain        radial_strain'
fmt = '%12d', '%12d', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e'
np.savetxt(stress_bond_outfname, np.column_stack((node_cross[-1], bond_cross[-1], angle_cross_list[-1], bond_mid_x_list[-1], bond_mid_y_list[-1], strain_bond_cross_x_list[-1], strain_bond_cross_y_list[-1], strain_bond_inner_circle[-1], radial_force_cross_bond_list[-1])), header = heading, fmt = fmt)

# writing total stress data : sum of all bonds that cross the inner boundary
stress_bond_outfname = base+folder+"txt/force/"+"inner_bndry_bond_radial_force_total_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if L != 64:
    stress_bond_outfname = base+folder+"txt/force/"+"inner_bndry_bond_radial_force_total_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
print('here all individual radial forces on all bonds that cross the inner boundary are written to file : ',stress_bond_outfname)    
heading = 'total nodes      radial stress'
fmt = '%12d', '%15.7e'
np.savetxt(stress_bond_outfname, np.column_stack((len(node_cross[-1]), total_radial_stress[-1])), header = heading, fmt = fmt)

###############################################################################    
#finding dipole nodes and those outside the inner circle
###############################################################################
print("Plotting begins below")

#reading the location of nodes in outer ring
if force_flag == 1:
    fname_ring = base+folder+'txt/area/ring_nodes_'+ "srand_" + str(srand) + "_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
    print("outer ring nodes filename: ",fname_ring)        
    ring_nodes = np.loadtxt(fname_ring)
    ring_nodes = ring_nodes[:,1]   # the 0th column is just a count
    ring_nodes = ring_nodes.astype(int)
    
# reading location of dipole nodes
#dip_node_arr = [119,120,134,135,136,151,152]
fname_dip = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond != 1:
    fname_dip = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64 and pbond == 1:
    fname_dip = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+ "srand_" + str(srand) + "_"+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64:
    fname_dip = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+ "srand_" + str(srand) + "_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if force_flag == 1:
    fname_dip = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+ "srand_" + str(srand) + "_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
    

print("dipole nodes filename: ",fname_dip)
dip_loc = np.loadtxt(fname_dip)
if num_center > 0:
    dip_node_arr = dip_loc[:,1]
    dip_node_arr = (np.rint(dip_loc[:,1])).astype(int)
else:
    dip_node_arr = np.array([0])


# to a first approximation, finding the coordination number of non dipole nodes: includes boundary nodes (should be removed later)
non_dip_nodes = np.setdiff1d(all_nodes, dip_node_arr)
non_dip_nodes = np.setdiff1d(non_dip_nodes, outer_nodes)
coord_num = conn_data[:,-1]
mean_coord = np.mean(coord_num[non_dip_nodes])
std_coord = np.std(coord_num[non_dip_nodes])

coord_num_outer_dip_nodes = []
for i in range(0,len(dip_node_arr)):
    if i%7 != 0:
        coord_num_outer_dip_nodes.append(coord_num[dip_node_arr[i]])
coord_num_outer_dip_nodes = np.array(coord_num_outer_dip_nodes)        

# writing coordination numbers in non dipole nodes
coord_outfname = base+folder+"txt/area/"+"mean_coord_num_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if L != 64:
    coord_outfname = base+folder+"txt/area/"+"mean_coord_num_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
print('here the mean coordination number is written to file : ', coord_outfname)    
heading = 'total nodes      mean z          std z             mean inner z           std inner z'
fmt = '%12d', '%15.5e', '%15.5e'
np.savetxt(coord_outfname, np.column_stack((len(non_dip_nodes), mean_coord, std_coord)), header = heading, fmt = fmt)


###############################################################################    
#plotting
###############################################################################

colors=np.full((num_pts), 'grey')
sc_size=np.full((num_pts), 0.2)


def absmax(x):
    maximum = abs(max(x))
    minimum = abs(min(x))
    
    if maximum > minimum:
        return maximum
    else:
        return minimum
    
def scaledown(x):
    for k in range(0,len(x)):
        if x[k] > 0.005:
            x[k] = np.exp(x[k])
        else:
            x[k] = 0.2
#        if x[k] <= 0.1:
#            x[k] = np.exp(x[k])/5
            
    return x
            
def binary(x):
    maxx = max(x)
    minn = np.min(x)
    for k in range(0,len(x)):
        if x[k] > 0.5*maxx:
            x[k] = 2
        elif x[k] < 0.1*minn:
            x[k] = 2
        else:
            x[k] = 0.2
            
    return x
    
    
    

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)



from mpl_toolkits.axes_grid1 import make_axes_locatable

log_flag = 0    # to apply log scale to strains
lw = 0.5 # line width of bonds in plot

#custom colormap:
#top = cm.get_cmap('Oranges_r', 128)
#bottom = cm.get_cmap('Blues', 128)

autumn_r_big = cm.get_cmap('Reds_r', 512)
autumn_small = ListedColormap(autumn_r_big(np.linspace(0.0, 0.7, 128)))
#autumn_small = ListedColormap(autumn_r_big(np.linspace(0.0, 1.0, 128)))

winter_big = cm.get_cmap('Blues', 512)
winter_small = ListedColormap(winter_big(np.linspace(0.3, 1.0, 128)))
#winter_small = ListedColormap(winter_big(np.linspace(0.0, 1.0, 128)))

top = cm.get_cmap('winter', 128)
bottom = cm.get_cmap('autumn_r', 128)

#newcolors = np.vstack((top(np.linspace(0, 1, 128)), bottom(np.linspace(0, 1, 128))))
newcolors = np.vstack((autumn_small(np.linspace(0, 1, 128)), winter_small(np.linspace(0, 1, 128))))
#newcmp = ListedColormap(newcolors, name='OrangeBlue')
newcmp = ListedColormap(newcolors, name='winter-cool')

# From Joe Kington: This one gives two different linear ramps:
import matplotlib.colors as mcolors
class MidpointNormalize(mcolors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
    
    
# =============================================================================
# Strain, bending-energy, and boundary-force visualizations
# =============================================================================

from scipy.interpolate import griddata
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib import ticker

plot_folder = base + folder + "png/strain/"

movie_flag = 0
color_flag = 0
test_flag = 0
bend_heat_flag = 0
test_force_flag = 1

if movie_flag:
    lw = 2.5


# -----------------------------------------------------------------------------
# Small plotting utilities
# -----------------------------------------------------------------------------

def sci_notation(value, pos=None):
    """Format colorbar labels as mantissa × 10^exponent."""
    if value == 0:
        return "0"

    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10**exponent
    return rf"${mantissa:.1f}\times10^{{{exponent}}}$"


def circular_line_segments(x, y):
    """Convert flattened x/y bond endpoint arrays into LineCollection segments."""
    return np.stack(
        (
            np.column_stack((x[0::2], y[0::2])),
            np.column_stack((x[1::2], y[1::2])),
        ),
        axis=1)


def transform_strain_values(values):
    """Apply optional log transform and compression scaling."""
    values = values.copy()

    if log_flag:
        positive = values > 0
        negative = values < 0

        values[positive] = -np.log(values[positive])
        values[negative] = np.log(np.abs(values[negative]))

    return np.where(values < 0, values * mu_c, values)


def apply_buckling_force_scaling():
    """Convert negative strain to force-like values for buckling runs."""
    global strain

    if buckle_flag:
        strain = np.where(strain < 0, strain * mu_c, strain)


def strain_line_collection(step):
    """Create strain-colored or binary-colored bond LineCollection."""
    threshold = abs(strain_thresh)
    raw_strain = strain_circle[step]
    plot_strain = transform_strain_values(raw_strain)

    x = x_loc_plot[step]
    y = y_loc_plot[step]
    segments = circular_line_segments(x, y)

    stretch = raw_strain > threshold
    compress = raw_strain < -threshold
    neutral = np.abs(raw_strain) < threshold

    line_widths = np.full(raw_strain.shape, lw)
    line_widths[stretch | compress] = 4 * lw

    if color_flag:
        vlim = np.max(np.abs(plot_strain))
        line_collection = LineCollection(segments, lw=line_widths, array=plot_strain, cmap=newcmp, norm=plt.Normalize(vmin=-vlim, vmax=vlim))
    else:
        line_colors = np.full(raw_strain.shape, "k", dtype=object)
        line_colors[stretch] = "b"
        line_colors[compress] = "r"
        line_colors[neutral] = "k"

        line_collection = LineCollection(segments, lw=line_widths, colors=line_colors)

    return line_collection, threshold


def node_plot_style(highlight_ring=False):
    """Return marker sizes and colors for network nodes."""
    node_colors = np.full(num_pts, "grey", dtype=object)
    node_sizes = np.full(num_pts, 0.2)

    if srand_flag and diff_flag:
        node_sizes[inner_nodes] = 3.0 if L == 128 else 10.0
        node_colors[inner_nodes] = "k"

    if highlight_ring:
        node_sizes[nodes_in_ring] = 20.0
        node_colors[nodes_in_ring] = "m"

    node_sizes[dip_node_arr] = 200.0 if movie_flag else (10.0 if L == 128 else 50.0)
    node_colors[dip_node_arr] = "g"

    node_sizes[outer_nodes] = 0.0
    node_colors[outer_nodes] = "w"

    node_sizes[boundary_nodes_circle] = 20.0
    node_colors[boundary_nodes_circle] = "m"

    if force_flag:
        node_sizes[ring_nodes] = 20.0
        node_colors[ring_nodes] = "b"

    return node_sizes, node_colors


def add_strain_colorbar(fig, ax, line_collection):
    """Attach strain colorbar when color_flag is enabled."""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = fig.colorbar(
        line_collection,
        ax=ax,
        format=ticker.FuncFormatter(fmt),
        cax=cax,
    )
    cbar.ax.tick_params(labelsize=30)
    cbar.locator = ticker.MaxNLocator(nbins=3)
    cbar.update_ticks()


def format_network_axes(ax):
    """Apply common axis formatting for network plots."""
    if movie_flag:
        ax.set_xlim([30, 50])
        ax.set_ylim([35, 50])

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(axis="both", which="both", labelsize=18, left=False, bottom=False)


def plot_filename(prefix, step, threshold_txt, include_quiver_scale=False, quiver_scale_txt=None):
    """Build output filename while preserving the old naming convention."""
    value_txt = force_txt if force_flag else rlen_txt

    if movie_flag:
        return f"{plot_folder}img{step:03d}.jpg"

    if force_flag:
        stem = (
            f"srand_{srand}_{L}_{num_center}_{num_dip}_{pbond_string}_{mu_str}_"
            f"{mu_c_str}_{kappa_str}_{tol_str}_{value_txt}_{threshold_txt}_{L}_{step}_{num - 1}")

    elif L == 128 and srand_flag:
        stem = (
            f"srand_{srand}_{L}_{num_center}_{num_dip}_{pbond_string}_{mu_str}_"
            f"{mu_c_str}_{kappa_str}_{tol_str}_{value_txt}_{threshold_txt}_{L}_{step}_{num - 1}")

    elif L == 128:
        stem = (
            f"{L}_{num_center}_{num_dip}_{pbond_string}_{mu_str}_{mu_c_str}_"
            f"{kappa_str}_{tol_str}_{value_txt}_{threshold_txt}_{L}_{step}_{num - 1}")

    elif srand_flag and pbond == 1:
        stem = (
            f"srand_{srand}_{num_center}_{num_dip}_{pbond_string}_{mu_str}_"
            f"{mu_c_str}_{kappa_str}_{tol_str}_{value_txt}_{threshold_txt}_{L}_{step}_{num - 1}")

    else:
        stem = (
            f"{num_center}_{num_dip}_{pbond_string}_{mu_str}_{mu_c_str}_"
            f"{kappa_str}_{tol_str}_{value_txt}_{threshold_txt}")

        if include_quiver_scale and quiver_scale_txt is not None:
            stem += f"_{quiver_scale_txt}"

        stem += f"_{L}_{step}_{num - 1}"

    return f"{plot_folder}{prefix}_{stem}.png"


# -----------------------------------------------------------------------------
# Circular strain plot
# -----------------------------------------------------------------------------

def plot_circle_strain(step):
    """Plot binary or continuous strain on the circular network."""
    line_collection, threshold = strain_line_collection(step)
    node_sizes, node_colors = node_plot_style(highlight_ring=False)

    fig, ax = plt.subplots(figsize=(15, 9), dpi=100)

    ax.add_collection(line_collection)
    ax.scatter(xpos[step], ypos[step], s=node_sizes, marker="o", color=node_colors, alpha=0.5)

    if color_flag:
        add_strain_colorbar(fig, ax, line_collection)

    format_network_axes(ax)
    plt.subplots_adjust(right=0.85, left=0.05, top=0.95, bottom=0.05)

    prefix = "circle_colored_strain" if color_flag else "circle_binary_strain"
    threshold_txt = f"{threshold:.2e}"
    output_file = plot_filename(prefix, step, threshold_txt)

    plt.savefig(output_file, dpi=100 if movie_flag else 1200)
    plt.close(fig)

    print("plotted to:", output_file)


# -----------------------------------------------------------------------------
# Bending-energy heatmap
# -----------------------------------------------------------------------------

def bending_energy_grid(grid_size=300):
    """Interpolate nodal bending energy to a regular grid for heatmap plotting."""
    xi = np.linspace(np.min(xpos[-1]), np.max(xpos[-1]), grid_size)
    yi = np.linspace(np.min(ypos[-1]), np.max(ypos[-1]), grid_size)
    xi, yi = np.meshgrid(xi, yi)

    zi = griddata((xpos[-1], ypos[-1]), bend_en_node, (xi, yi), method="cubic")

    return xi, yi, zi


def bending_energy_ticks():
    """Return manually tuned bending-energy colorbar ticks."""
    tick_map = {
        (5, 1e-6, 113): [0, 1e-8, 2e-8, 3e-8],
        (5, 1e-5, 113): [0, 1e-7, 2e-7, 3e-7],
        (5, 1e-4, 113): [0, 1e-6, 2e-6, 3e-6],
        (35, 1e-6, 113): [0, 1e-7, 2e-7, 3e-7, 4e-7, 5e-7],
        (35, 1e-5, 113): [0, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6],
        (35, 1e-4, 113): [0, 1e-5, 2e-5, 3e-5, 4e-5],
    }

    return tick_map.get((num_center, kappa, srand), [0, 1e-8, 2e-8, 3e-8, 4e-8])


def plot_bending_heatmap(step):
    """Plot bending-energy heatmap underneath the network strain plot."""
    line_collection, threshold = strain_line_collection(step)
    node_sizes, node_colors = node_plot_style(highlight_ring=True)

    xi, yi, zi = bending_energy_grid()

    fig, ax = plt.subplots(figsize=(15, 9), dpi=300)

    cmap_orange = LinearSegmentedColormap.from_list(
        "fade_white_orange",
        [(1, 1, 1), (1, 0.5, 0)])

    heatmap = ax.contourf(xi, yi, zi, levels=100, cmap=cmap_orange, vmin=0)

    cbar = fig.colorbar(heatmap, ax=ax)
    cbar.set_label("Bending Energy", fontsize=25, labelpad=10)
    cbar.ax.set_ylim([0, np.nanmax(zi)])
    cbar.ax.tick_params(labelsize=20)
    cbar.set_ticks(bending_energy_ticks())
    cbar.ax.yaxis.set_major_formatter(ticker.FuncFormatter(sci_notation))

    ax.add_collection(line_collection)
    ax.scatter(xpos[step], ypos[step], s=node_sizes, marker="o", color=node_colors, alpha=0.5 )

    if color_flag:
        add_strain_colorbar(fig, ax, line_collection)

    format_network_axes(ax)
    plt.subplots_adjust(right=0.85, left=0.05, top=0.95, bottom=0.05)

    threshold_txt = f"{threshold:.2e}"
    output_file = plot_filename("circle_heatmap", step, threshold_txt)

    plt.savefig(output_file, dpi=300)
    plt.close(fig)

    print("plotted to:", output_file)


# -----------------------------------------------------------------------------
# Boundary-force plot
# -----------------------------------------------------------------------------

def boundary_force_filename():
    """Return input file containing boundary-node forces."""
    if only_last_flag:
        stem = (
            f"{num_center}_{num_dip}_{pbond_string}_{tol_str}_{kappa_str}_"
            f"{rlen_txt}_{mu_str}_{mu_c_str}_{num - 1}")

    elif force_flag or (srand_flag and pbond == 1 and L == 64):
        stem = (
            f"srand_{srand}_{num_center}_{num_dip}_{pbond_string}_{tol_str}_"
            f"{kappa_str}_{rlen_txt}_{mu_str}_{mu_c_str}_{num - 1}")

    elif srand_flag and L != 64:
        stem = (
            f"srand_{srand}_{L}_{num_center}_{num_dip}_{pbond_string}_{tol_str}_"
            f"{kappa_str}_{rlen_txt}_{mu_str}_{mu_c_str}_{num - 1}")

    else:
        stem = (
            f"{num_center}_{num_dip}_{pbond_string}_{tol_str}_{kappa_str}_"
            f"{rlen_txt}_{mu_str}_{mu_c_str}_{num - 1}")

    return f"{base}{folder}txt/area/bndry_node_force_{stem}.txt"


def boundary_quiver_data(step):
    """Sort boundary nodes by angle and return smoothed, clipped force vectors."""
    x_boundary = xpos[step, boundary_nodes_circle]
    y_boundary = ypos[step, boundary_nodes_circle]

    x0 = np.mean(x_boundary)
    y0 = np.mean(y_boundary)

    theta = np.mod(np.arctan2(y_boundary - y0, x_boundary - x0), 2 * np.pi)
    sort_idx = np.argsort(theta)
    sample_idx = np.arange(0, len(sort_idx), 2)

    x_plot = x_boundary[sort_idx][sample_idx]
    y_plot = y_boundary[sort_idx][sample_idx]

    fx_plot = smooth_circular_vectors(force_x_bndry[sort_idx])[sample_idx]
    fy_plot = smooth_circular_vectors(force_y_bndry[sort_idx])[sample_idx]

    if num_center in {5, 15} and pbond == 0.55:
        quiver_scale = 1e-5
    elif num_center == 20 and pbond == 0.55:
        quiver_scale = 1e-6
    else:
        quiver_scale = 100 * np.mean(np.hypot(fx_plot, fy_plot))

    if num_center == 5 and kappa == 1e-6:
        max_mag = 1e-6
    elif num_center == 5 and kappa == 1e-5:
        max_mag = 2e-6
    elif num_center == 20:
        max_mag = 1.5e-7
    else:
        max_mag = 2e-6

    magnitude = np.hypot(fx_plot, fy_plot)
    scale = np.minimum(1.0, max_mag / np.maximum(magnitude, 1e-30))

    return x_plot, y_plot, fx_plot * scale, fy_plot * scale, quiver_scale


def plot_boundary_forces(step):
    """Plot network strain with boundary force vectors."""
    line_collection, threshold = strain_line_collection(step)
    node_sizes, node_colors = node_plot_style(highlight_ring=False)

    x_q, y_q, fx_q, fy_q, quiver_scale = boundary_quiver_data(step)

    fig, ax = plt.subplots(figsize=(15, 9), dpi=100)

    ax.add_collection(line_collection)
    ax.scatter(xpos[step], ypos[step], s=node_sizes, marker="o", color=node_colors, alpha=0.5)
    
    ax.scatter(xpos[step, dip_node_arr], ypos[step, dip_node_arr], s=20, marker="o", color="g", zorder=50)

    ax.quiver(x_q, y_q, fx_q, fy_q, scale=quiver_scale, zorder=10)

    if color_flag:
        add_strain_colorbar(fig, ax, line_collection)

    ax.set_axis_off()
    format_network_axes(ax)
    plt.subplots_adjust(right=0.85, left=0.05, top=0.95, bottom=0.05)

    threshold_txt = f"{threshold:.2e}"
    quiver_scale_txt = f"{quiver_scale:.4f}"

    output_file = plot_filename(
        "force_bndry",
        step,
        threshold_txt,
        include_quiver_scale=True,
        quiver_scale_txt=quiver_scale_txt,
    )

    plt.savefig(output_file, dpi=100 if movie_flag else 1200)
    plt.close(fig)

    print("plotted to:", output_file)


# -----------------------------------------------------------------------------
# Region setup and plot execution
# -----------------------------------------------------------------------------

apply_buckling_force_scaling()

ring_rad = 25 if L == 128 else 12
dist_all = np.hypot(xpos[0] - xpos[0, center], ypos[0] - ypos[0, center])
nodes_in_ring = np.flatnonzero((dist_all > ring_rad - 0.5) & (dist_all < ring_rad + 0.5))

if test_flag:
    for step in range(len(strain_circle)):
        plot_circle_strain(step)

if bend_heat_flag:
    for step in range(len(strain_circle) - 1, len(strain_circle)):
        plot_bending_heatmap(step)

if test_force_flag:
    force_fname = boundary_force_filename()
    print("Boundary-node force file:", force_fname)

    force_data = np.loadtxt(force_fname)
    force_x_bndry = force_data[:, 1]
    force_y_bndry = force_data[:, 2]
    perp_force_bndry = force_data[:, -1]

    for step in range(len(strain) - 1, len(strain)):
        plot_boundary_forces(step)


###############################################################################    
# making a ring to find the decay of radial stress vs radius
# should only do it when the dipole is in the center!!
###############################################################################
outer_radius = 25
radius_bins = np.arange(1,outer_radius+1)
x_center = xpos[0, center]
y_center = ypos[0, center]
drr = 0.5

# Vectorized ring membership for all radii.
dist_from_center = np.hypot(x_center - xpos[-1], y_center - ypos[-1])
radius_ring_nodes = [
    np.flatnonzero((dist_from_center > radius - drr) & (dist_from_center < radius + drr))
    for radius in radius_bins
]


#==============================================================================
# testing if we can produce the rings :: do only for when 1 dipole is at center
#==============================================================================
# plotting a circle as well

ring_flag = 0
colors=np.full((num_pts), 'grey')
sc_size=np.full((num_pts), 40)

# setting different color for each ring
cmap = plt.get_cmap('jet')
num_colors = len(radius_ring_nodes)
colors_list = [cmap(i / (num_colors - 1)) for i in range(num_colors)]
colors_list = ['b','c','g','y','m','r','b','c','g','y','m','r','b','c','g','y','m','r','b','c','g','y','m','r','b']

for i in range(0,len(radius_ring_nodes)):
    colors[radius_ring_nodes[i]] = colors_list[i]

if ring_flag == 1:
    # for i in np.arange(0,len(strain_circle)):
#    for i in np.arange(0,2):
    for i in np.arange(0,1):
        plt_thrshld = abs(strain_thresh)       # threshold to plot red adn blue colors for expanding and contracting bonds
        plt_thrshld_str = "%.2e" % plt_thrshld 
                
        x = x_loc_plot[i]
        y = y_loc_plot[i]
        
        # initial locations
        x_init = x_loc_plot[0]
        y_init = y_loc_plot[0]

        xy = []

        for j in np.arange(0,len(x)-1,2):
            xy.append([(x[j], y[j]),(x[j+1],y[j+1])])
            
# better method to remove certain bonds so that we can get a nice colorbar (also must remove the colorbar for these lines)
        high_strain = np.array(xy) # copy of xy to be potted separately to have sanity of colorbar
        remove_loc = np.squeeze(np.where(strain_circle[i] < 0.5*np.min(strain_circle[i])))
#        remove_loc = np.squeeze(np.where(strain < 0.5*np.min(strain)))

        lc_arr = np.zeros_like(strain_circle[i])
#        lc_arr = np.zeros_like(strain)
        lc_arr[:] = strain_circle[i][:]
#        lc_arr[:] = strain[:]
        if log_flag == 1:
            for k in range(len(lc_arr)):
                if lc_arr[k] > 0: lc_arr[k] = -1*np.log(lc_arr[k])
                if lc_arr[k] < 0: lc_arr[k] = np.log(abs(lc_arr[k]))
#                lc_arr[lc_arr < 0.0] = -1*np.log(abs(lc_arr))
        
        col = np.array(['k']*len(strain_circle[0]))
        coltest = np.array([(0.0,0.0,0.0,1.0)]*len(strain_circle[0]))

        lwidth = np.array([lw]*len(strain_circle[i]))

        lc_arr = list(lc_arr)
        lwidth = list(lwidth)
        remove_loc = np.atleast_1d(np.array(remove_loc))

        lc_arr = np.array(lc_arr)
        for kk in range(0,len(lc_arr)):
            if lc_arr[kk] < 0:
                lc_arr[kk] = lc_arr[kk]*mu_c
        lwidth = np.array(lwidth)

        stretch_index = np.where(lc_arr>plt_thrshld)
        compress_index = np.where(lc_arr<(-1.0*plt_thrshld))

        stretch_index = np.where(strain_circle[i]>plt_thrshld)
        compress_index = np.where(strain_circle[i]<(-1.0*plt_thrshld))

        zero_index = np.where(abs(lc_arr) < plt_thrshld)
        zero_index = np.where(abs(strain_circle[i]) < plt_thrshld)

        strain_pos_max = np.max(lc_arr)
        strain_neg_max = np.min(lc_arr)
        
        if abs(strain_neg_max) > strain_pos_max:
            vlim = -1.*strain_neg_max
        else:
            vlim = strain_pos_max
        
        if (i==0):
            coltest[stretch_index,2] = 1.0
            coltest[compress_index,0] = 1.0
    
        if srand_flag == 1 and diff_flag == 1:
            sc_size[inner_nodes] = 40.0
            # colors[inner_nodes] = 'k'
            if L == 128:
                sc_size[inner_nodes] = 3.0
                # colors[inner_nodes] = 'k'

        sc_size[dip_node_arr] = 50.0
        if movie_flag == 1:        
            sc_size[dip_node_arr] = 200.0              # for movie
        if L ==128:
            sc_size[dip_node_arr] = 10.0              # for large L            
        # colors[dip_node_arr] = 'g'

        # sc_size[boundary_nodes] = 20.0
        # colors[boundary_nodes] = 'm'
        
        sc_size[outer_nodes] = 0.0
        colors[outer_nodes] = 'w'

        # sc_size[boundary_nodes_circle] = 20.0
        # colors[boundary_nodes_circle] = 'm'

        col[stretch_index] = 'b'
        col[compress_index] = 'r'
        col[zero_index] = 'k'
        lwidth[stretch_index] = 4*lw
        lwidth[compress_index] = 4*lw
                        
#        lc = LineCollection(xy, lw = lwidth, cmap = newcmp, norm=plt.Normalize(vmin=-vlim,vmax=vlim))
        if color_flag == 1:
            lc = LineCollection(xy, lw = lwidth, array = lc_arr, cmap = newcmp, norm=plt.Normalize(vmin=-vlim,vmax=vlim))#, alpha = 0.5)#,vmin=-0.05))
        else:
            lc = LineCollection(xy, lw = lwidth, colors = col)#, alpha = 0.5)#,vmin=-0.05))
        
        fig = plt.figure(figsize=(15, 9), dpi=100)
        ax1 = fig.add_subplot(1, 1, 1)

        if color_flag == 1:
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes("right", size="5%", pad=0.05)
    
            cbar = fig.colorbar(lc, ax = ax1, format=ticker.FuncFormatter(fmt), cax=cax)
            for t in cbar.ax.get_yticklabels():
                    t.set_fontsize(30)
            tick_locator = ticker.MaxNLocator(nbins=3)
            cbar.locator = tick_locator
            cbar.update_ticks()        
        ax1.add_collection(lc)
        
        ax1.scatter(xpos[i], ypos[i], s = sc_size, marker = "o", color = colors, alpha = 0.8)
        
        ax1.set_aspect('equal')
#        ax1.set_title("Strain", fontsize=16)

        ax1.tick_params(axis='both', which='major', labelsize = 18, left = False, bottom = False)
        ax1.tick_params(axis='both', which='minor', labelsize = 18,  left = False, bottom = False)
        
        ax1.set_xticks([])
        ax1.set_yticks([])

        plt.subplots_adjust(right=0.995, left = 0.0, top = 0.9, bottom = 0.05)
        plt.subplots_adjust(right=0.85,left=0.05,top=0.95,bottom=0.05)

        plot_fname = plot_folder + 'ring_binary_strain_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
        if srand_flag == 1 and pbond == 1:
            plot_fname = plot_folder + 'ring_binary_strain_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if L ==128:
            plot_fname = plot_folder + 'ring_binary_strain_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            # if srand_flag == 1 and num_center > 1 and pbond == 1:
            if srand_flag == 1:
                plot_fname = plot_folder + 'ring_binary_strain_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if color_flag == 1:
            plot_fname = plot_folder + 'ring_colored_strain_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if L != 64:
                plot_fname = plot_folder + 'ring_colored_strain_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if srand_flag == 1:
                plot_fname = plot_folder + 'ring_colored_strain_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
                
        if force_flag == 1:
            plot_fname = plot_folder + 'ring_binary_strain_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+force_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if movie_flag != 1:
            plt.savefig(plot_fname, dpi = 300)
                
        print("plotted to: ", plot_fname)
        plt.clf()
        plt.close('all')


#==============================================================================


#==============================================================================
# writing the number of bonds out
#==============================================================================
# bond_count_inner_p1 = 0 # stores the number of bonds in p=1 simulations
inner_conn = conn_node[inner_nodes, :3]
valid_inner = inner_conn != unconnected_const
inner_node_idx, inner_bond_idx = np.where(valid_inner)
connected_inner = inner_conn[inner_node_idx, inner_bond_idx]

dist_center = np.hypot(xpos[0, center] - xpos[0, connected_inner], ypos[0, center] - ypos[0, connected_inner])
inside_inner_radius = dist_center <= inner_radius

temp_bond = np.column_stack((inner_nodes[inner_node_idx[inside_inner_radius]], inner_bond_idx[inside_inner_radius])).astype(int)
count_inner_bonds = len(temp_bond)


#==============================================================================
# plotting the displacements
#==============================================================================
x_disp = xpos[-1,:] - xpos[0,:]
y_disp = ypos[-1,:] - ypos[0,:]
disp = np.sqrt(x_disp**2 + y_disp**2)

init_radial_dx = xpos[0,:] - xpos[0,center] 
init_radial_dy = ypos[0,:] - ypos[0,center] 
init_radial_dist = np.sqrt(init_radial_dx**2 + init_radial_dy**2)

final_radial_dx = xpos[-1,:] - xpos[-1,center] 
final_radial_dy = ypos[-1,:] - ypos[-1,center] 
final_radial_dist = np.sqrt(final_radial_dx**2 + final_radial_dy**2)

radial_disp = final_radial_dist - init_radial_dist

# xlim = radius-1+dr
xlim = radius+dr

# scatter plot of all nodes' radial displacement
fig = plt.figure(figsize=(8, 4), dpi=300)
plt.scatter(init_radial_dist, radial_disp)
plt.xlabel('Distance from center', fontsize = 15)
plt.ylabel('Radial Displacement', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.ylim([-5e-5,5e-5])
plt.xlim([1,xlim])
plt.subplots_adjust(top=0.92, left = 0.13, bottom = 0.15, right = 0.98)
if L == 64:
    plt.savefig(base+folder+"png/displacement_scatter_"+'srand_'+str(srand)+"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
elif L != 64:
    plt.savefig(base+folder+"png/displacement_scatter_"+ 'srand_'+ str(srand) + "_" + str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
plt.close()
plt.clf()



fig = plt.figure(figsize=(8, 4), dpi=300)
plt.scatter(init_radial_dist, radial_disp, alpha = 0.3)
plt.xlabel('Distance from center', fontsize = 15)
plt.ylabel('Radial Displacement', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
if pbond == 1 and L == 64:
    plt.ylim([-0.004,0])   # for p = 1; N = 1
if pbond == 0.9:
    plt.ylim([-0.02,0.01])   # for p = 0.9; N = 1
if pbond == 0.8:
    plt.ylim([-0.025,0.015])   # for p = 0.8; N = 1

plt.xlim([inner_radius,xlim])
plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
if L == 64:
    plt.savefig(base+folder+"png/displacement_scatter_outer_region_"+'srand_'+str(srand)+"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
elif L != 64:
    plt.savefig(base+folder+"png/displacement_scatter_outer_region_"+ 'srand_'+ str(srand) + "_" + str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
plt.close()
plt.clf()


# binning :: this is for p<1 networks
bins = np.arange(1, xlim, 1)
dd = 0.5

# Broadcasted bin mask keeps the original [bin-dd, bin+dd) definition.
bin_mask = (init_radial_dist[None, :] >= (bins[:, None] - dd)) & (init_radial_dist[None, :] < (bins[:, None] + dd))
counts = bin_mask.sum(axis=1)
radial_disp_masked = np.where(bin_mask, radial_disp[None, :], np.nan)

index_list = [np.flatnonzero(mask) for mask in bin_mask]
binned_disp = [radial_disp[idx] for idx in index_list]
binned_mean_disp = np.divide(np.nansum(radial_disp_masked, axis=1), counts, out=np.full(len(bins), np.nan), where=counts > 0)
binned_median_disp = np.nanmedian(radial_disp_masked, axis=1)
binned_std_disp = np.nanstd(radial_disp_masked, axis=1)


# writing results to file for later analysis
# writing all the node data
disp_outfname = base+folder+"txt/displacement/"+"bndry_node_radial_disp_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if L != 64:
    disp_outfname = base+folder+"txt/displacement/"+"bndry_node_radial_disp_"+ 'srand_'+str(srand)+"_" + str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
print('here all individual radial displacements on all nodes are written to file : ',disp_outfname)    
heading = 'node       Initial Radial dist        final Radial dist        Radial disp'
fmt = '%12d', '%15.7e', '%15.7e', '%15.7e'
np.savetxt(disp_outfname, np.column_stack((np.arange(0,num_pts), init_radial_dist, final_radial_dist, radial_disp)), header = heading, fmt = fmt)

# writing binned data
disp_mean_outfname = base+folder+"txt/displacement/"+"bndry_node_radial_disp_mean_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if L != 64:
    disp_mean_outfname = base+folder+"txt/displacement/"+"bndry_node_radial_disp_mean_"+ 'srand_'+str(srand)+"_" + str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('here all mean (after binning) radial displacements on all nodes are written to file : ',disp_mean_outfname)    
heading = 'bin mid       Mean Radial dist        Meadian Radial dist        Std Radial dist'
fmt = '%12d', '%15.7e', '%15.7e', '%15.7e'
np.savetxt(disp_mean_outfname, np.column_stack((bins, binned_mean_disp, binned_median_disp, binned_std_disp)), header = heading, fmt = fmt)

# scatter plot - outer region
r1 = inner_radius + 1
r2 = 25
if L != 64:
    # r1 = 24+1
    # r1 = 12+1
    r2 = 48
    # r2 = 50
    
bulk_m = 0.5*np.sqrt(3)*mu
shear_m = 0.25*np.sqrt(3)*mu

if pbond == 1:
    sigma1 = -5.8e-4
    sigma1_outer = -5.8e-4

elif pbond == 0.9:
    sigma1 = -6.2e-4
    sigma1_outer = -4.4e-4

    shear_m = 3.04e-1
    bulk_m = 6.08e-01

elif pbond == 0.8:
    sigma1 = -6.6e-4
    sigma1_outer = -9.2e-4

    shear_m = 1.87e-1
    bulk_m = 3.74e-01
    
elif pbond == 0.55:
    sigma1 = -9.0e-4
    sigma1_outer = -3.2e-3
    
fit_x = np.linspace(bins[0],bins[-1],100)
fit_x_outer = np.linspace(r1, r2, 100)     # outer region x range
u_r = ( (sigma1 * r1**2) / (2*(bulk_m * r1**2+ shear_m * r2**2)) ) * ( (r2**2 - fit_x**2)/fit_x )       #  theory ar + b/r fit for whole range 
u_r_outer = ( (sigma1_outer * r1**2) / (2*(bulk_m * r1**2+ shear_m * r2**2)) ) * ( (r2**2 - fit_x_outer**2)/fit_x_outer )       #  theory ar + b/r fit for whole range 


from scipy.optimize import curve_fit
def func(x,a):
    # return (a) * (r1**2) / (2*(bulk_m * r1**2+ shear_m * r2**2))  * ( (r2**2 - fit_x_outer**2)/fit_x_outer )
    # return (a) * (r1**2) / (2*(bulk_m * r1**2+ shear_m * r2**2))  * ( (r2**2 - x**2)/x )    
    return (a) * (r1**2) / (2*((np.sqrt(3)/2) * r1**2 + ((np.sqrt(3)/4) * r2**2)))  * ( (r2**2 - x**2)/x )

bin_slice_val = r1-1
# if L == 128:
#     bin_slice_val = 24

xdata = bins[bin_slice_val:]
ydata = binned_mean_disp[bin_slice_val:]

# f2str = '{:.2f}'
popt, pcov = curve_fit(func, xdata, ydata, maxfev=10000)
# label = '$%s x^{-%s}$' % tuple(f2str.format(t) for t in [popt[0], popt[1]])
print('Fit parameters: ', popt)

u_r_outer_fit = func(fit_x_outer, popt[0])

# scatter plot of binned data (wrt initial radial distance) and also fits
fig = plt.figure(figsize=(8, 4), dpi=300)
plt.scatter(bins, binned_mean_disp)
if L == 64:
    plt.plot(fit_x, u_r, c = 'r')
plt.xlabel('Distance from center', fontsize = 15)
plt.ylabel('Radial Displacement', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.ylim([-0.004,0])   # for p = 1; N = 1
# plt.ylim([-0.02,0.01])   # for p = 0.9; N = 1
# plt.ylim([-0.025,0.015])   # for p = 0.8; N = 1

# plt.xlim([12,radius+dr])
plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
if L == 64:
    plt.savefig(base+folder+"png/displacement_scatter_mean_"+'srand_'+str(srand)+"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
if L != 64:
    plt.savefig(base+folder+"png/displacement_scatter_mean_"+'srand_'+str(srand)+"_" + str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
    
plt.close()
plt.clf()

binned_mean_disp_outer = binned_mean_disp[bin_slice_val:]
binned_std_disp_outer = binned_std_disp[bin_slice_val:]
bins_outer = bins[bin_slice_val:]

fig = plt.figure(figsize=(8, 4), dpi=300)
# plt.scatter(bins, binned_mean_disp)
plt.scatter(bins_outer, binned_mean_disp_outer)
# plt.errorbar(bins_outer, binned_mean_disp_outer, yerr= binned_std_disp_outer, ls = "None", marker = 'o')
# plt.plot(fit_x_outer, u_r_outer, c = 'r', label = '$ \\Sigma_{1} = %0.2e $' % sigma1_outer)
plt.plot(fit_x_outer, u_r_outer_fit, c = 'b', label = '$ \\Sigma_{1} = %0.2e $' % popt[0])
plt.xlabel('Distance from center', fontsize = 15)
plt.ylabel('Radial Displacement', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.ylim([-0.004,0])   # for p = 1; N = 1
if pbond == 1 and num_center == 1 and L == 64:
    plt.ylim([-0.004,0.0])   # for p = 0.9; N = 1
# if pbond == 1 and num_center == 2:
#     plt.ylim([-0.008,0.0])   # for p = 0.9; N = 1
if pbond == 0.9:
    plt.ylim([-0.005,0.0])   #   for p = 0.9; N = 1
if pbond == 0.8:
    plt.ylim([-0.007,0.0])   # for p = 0.8; N = 1
if pbond == 0.55:
    plt.ylim([-0.025,0.0])   # for p = 0.8; N = 1

if L == 64:
    plt.ylim([np.min(binned_mean_disp_outer)+0.1*np.min(binned_mean_disp_outer),0.0+1e-4])   # for p = 0.8; N = 1

plt.xlim([inner_radius+0.5,xlim])
# plt.xlim([12+0.5,xlim])
plt.subplots_adjust(top=0.92, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend()
if L == 64:
    plt.savefig(base+folder+"png/displacement_scatter_outer_region_mean_"+'srand_'+str(srand)+"_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data
if L != 64:
    plt.savefig(base+folder+"png/displacement_scatter_outer_region_mean_"+'srand_'+str(srand)+"_" + str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data    

plt.close()
plt.clf()

# writing fit to sigma1 or sigma1/alpha_m to file
disp_outfname = base+folder+"txt/displacement/"+"fit_radial_disp_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if L != 64:
    disp_outfname = base+folder+"txt/displacement/"+"fit_radial_disp_"+'srand_'+str(srand)+"_" + str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
print('here the fits are written to file : ',disp_outfname)    
heading = 'pbond         kappa             N        Sigma1/alpha_m'
fmt = '%12d', '%10.2e', '%12d', '%15.7e'
np.savetxt(disp_outfname, np.column_stack((pbond, kappa, num_center, popt[0])), header = heading, fmt = fmt)


center_list = dip_node_arr[::7]  # take every 7th element starting at index 0
print(center_list)

#==============================================================================
# distance of dipole from center for 1 dipole case

if pbond == 1 and num_center == 1:
    dist_x = xpos[-1,dip_node_arr[0]] - xpos[0,center] 
    dist_y = ypos[-1,dip_node_arr[0]] - ypos[0,center] 
    
    dist_val = np.sqrt(dist_x**2 + dist_y**2)
    
    plt.savetxt()
    
