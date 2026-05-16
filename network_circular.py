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
def angle(x1,y1,x2,y2):
    dx = x2-x1
    dy = y2-y1
    return(math.atan2(dy,dx))


L = 64
#L = 65
#L = 16
# L = 128

srand_flag = 1
if srand_flag == 1:
    srand = 113


num_pts = L*L
num_center = 1
# num_center = 2
# num_center = 3
# num_center = 4
num_center = 5
#num_center = 6
#num_center = 7
# num_center = 10
# num_center = 15 
# #num_center = 19
num_center = 20
# num_center = 25 
# num_center = 30
#num_center = 32
#num_center = 33
#num_center = 34
# num_center = 35
# num_center = 40
# num_center = 45
# num_center = 50
#num_center = 60
#num_center = 70
#num_center = 100
#num_center = 200

# num_center = 160

num_dip = num_center*6
# num = 1 + 1
#num = 5 + 1                                                                           # the total number of force steps
num = 10 + 1                                                                           # the total number of force steps
#num = 15 + 1                                                                          # the total number of force steps for p = 0.7, kappa = 5e-5, rlen = 0.8
#num = 100 + 1                                                                         # the total number of force steps

#num = 1000 + 1

pbond = 1
# pbond = 0.9
# pbond = 0.8
# pbond = 0.7
# pbond = 0.75
#pbond = 0.5
pbond = 0.55
# pbond = 0.6

#pbond = 0.61

#pbond = 0.58
#pbond = 0.59
#pbond = 0.63
#pbond = 0.69
#pbond = 0.7
#pbond = 0.81
#pbond = 0.9


mu = 1
mu_c = 1
# mu = 1e-3
# mu_c = 1e-3


tol = 2e-6
tol = 1e-6
tol = 1e-7

kappa = 0
kappa = 1e-6
#kappa = 2e-6
#kappa = 5e-6
# kappa = 1e-5
#kappa = 2e-5
#kappa = 5e-5

# kappa = 1e-4
#kappa = 1e-2

force_val = 0.01

strain_thresh = 1e-2
strain_thresh = 1e-3
# strain_thresh = 1e-4
#strain_thresh = 5e-4
# strain_thresh = 1e-5
strain_thresh = 1e-6
strain_thresh = 1e-7
# strain_thresh = 1e-8
#

if pbond == 1:
    strain_thresh = 1e-3

rlen_flag = 0
rlen_hex_flag = 0
rlen_hex_sep_flag = 0
sep_flag = 0
pbd_flag = 0
buckle_flag = 0
radial_flag = 0
inner_radial_flag = 0

inner_radial_hex_arp_flag = 0         # for all hex bonds being present
inner_radial_arp_flag_bash = 0
radial_only_flag = 0

cluster_flag = 0
cluster_new_flag = 0         # should be 0 for radial and hex cases from the cluster

cluster_radial_flag = 1    # only radial bonds around central node are present
hex_flag = 0   # these are for cluster 100 simulations where the dipoles have all 12 bonds present - including the outer hex bonds
hex_rand_flag = 0    # azimuthal bonds randomly removed

test_hex = 0     # for testing if the hex code is working right locally
test_hex_rand = 0

anisotropic_flag = 0
if anisotropic_flag == 1:
    num_dip = num_center

force_flag = 0
if force_flag == 1:
    force_val = 0.01115
    force_txt = "%.9f" % force_val

if buckle_flag == 1:
    mu_c = 0.5
    mu_c = 0.1
    
rlen = 0.9
rlen_txt = "%.4f" % rlen


diff_flag = 1
if diff_flag == 1:
    diff_network1 = 111111
    # diff_network1 = 222222
#    diff_network1 = 666666
    # diff_network1 = 444444
    # diff_network1 = 900000

    # diff_network2 = 111111
#    diff_network2 = 222222
#    diff_network2 = 333333
    # diff_network2 = 444444
#    diff_network2 = 555555
#    diff_network2 = 666666

    diff_network2 = 101111
    # diff_network2 = 111111
    # diff_network2 = 121111
    # diff_network2 = 131111
    # diff_network2 = 141111
    # diff_network2 = 151111
    # diff_network2 = 161111
    # diff_network2 = 171111
    # diff_network2 = 181111
    # diff_network2 = 191111

    # diff_network2 = 202222
    # diff_network2 = 212222
#    diff_network2 = 222222
    # diff_network2 = 232222
#    diff_network2 = 242222
    # diff_network2 = 252222
    # diff_network2 = 262222
    # diff_network2 = 272222
    # diff_network2 = 282222
#    diff_network2 = 292222

#    diff_network2 = 616666
    # diff_network2 = 454444

    # diff_network2 = 900003

    diff_network = str(diff_network1)+","+str(diff_network2)
    
auto_flag = 0
if auto_flag == 1 and pbond < 1:
    import sys
    print(sys.argv[1], " ** ", sys.argv[2], " ** ", sys.argv[3])#, "**", sys.argv[1])
    diff_network = sys.argv[1]
    num_center = int(sys.argv[2])    
    num_dip = num_center*6    
    if srand_flag == 1:
        srand = int(sys.argv[3])
    print("Diff network: ", diff_network)

elif auto_flag == 1 and pbond == 1:
    import sys
    print(sys.argv[1], " ** ", sys.argv[2])#, " ** ", sys.argv[3])#, "**", sys.argv[1])
    num_center = int(sys.argv[1])    
    num_dip = num_center*6    
    if srand_flag == 1:
        srand = int(sys.argv[2])


pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
kappa_str = "%.2e" % kappa # to write the filename
if kappa == 0:
    kappa_str = "0.00e+00" # to write the filename
    
force_str = "%.4f" % force_val # to write the filename
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"


folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'
if pbond == 1 and radial_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_radial/'
if pbond == 1 and inner_radial_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial/'
if pbond == 1 and inner_radial_hex_arp_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'     # this is for hex bonds enforeced
if pbond == 1 and radial_only_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/'
if pbond == 1 and cluster_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/'
if pbond == 1 and cluster_new_flag == 1:
    folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'
if pbond == 1 and cluster_radial_flag == 1:
    folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'
if pbond == 1 and test_hex == 1:
    folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'
if force_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_force/'
if anisotropic_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_arp_anisotropic/'

if (pbond != 1):
    ranseed = '667720,601210'
#    ranseed = '152240,193915'
#    ranseed = '167720,101210'
#    ranseed = '712240,203915'
#    ranseed = '537206,701210'

#    ranseed = '399020,800099'
#    ranseed = '367725,927820'
    if kappa == 2e-7:
        kappa_fname = 'kappa2_2e-7/'
    elif kappa == 5e-7:
        kappa_fname = 'kappa2_5e-7/'
    elif kappa == 1e-6:
        kappa_fname = 'kappa2_e-6/'
    elif kappa == 2e-6:
        kappa_fname = 'kappa2_2e-6/'
    elif kappa == 5e-6:
        kappa_fname = 'kappa2_5e-6/'
    elif kappa == 1e-5:
        kappa_fname = 'kappa2_e-5/'
    elif kappa == 2e-5:
        kappa_fname = 'kappa2_2e-5/'
    elif kappa == 5e-5:
        kappa_fname = 'kappa2_5e-5/'
    elif kappa == 1e-4:
        kappa_fname = 'kappa2_e-4/'
    elif kappa == 2e-4:
        kappa_fname = 'kappa2_2e-4/'
    elif kappa == 5e-4:
        kappa_fname = 'kappa2_5e-4/'
    elif kappa == 1e-3:
        kappa_fname = 'kappa2_e-3/'
    elif kappa == 2e-3:
        kappa_fname = 'kappa2_2e-3/'
    elif kappa == 5e-3:
        kappa_fname = 'kappa2_5e-3/'
    elif kappa == 1e-2:
        kappa_fname = 'kappa2_e-2/'
    elif kappa == 0:
        kappa_fname = 'kappa2_0/'
    
    if srand_flag == 1 and diff_flag == 0:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'+ kappa_fname +ranseed+"/"+str(srand)+"/"   # different dipole positions in the same network
    elif diff_flag == 1 and srand_flag == 1 and inner_radial_flag == 0 and inner_radial_hex_arp_flag == 0 and cluster_flag == 0 and \
         inner_radial_arp_flag_bash == 0 and cluster_new_flag == 0 and cluster_radial_flag == 0 and hex_flag == 0 and hex_rand_flag == 0 and \
         test_hex == 0 and test_hex_rand == 0 and anisotropic_flag == 0:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # different dipole positions in the same network
    elif radial_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_radial/'+ kappa_fname +ranseed+"/"   # only radial bonds are used for force dipole
    elif inner_radial_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # only radial bonds are used for force dipole and inner network is difference from outer
    elif inner_radial_hex_arp_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # hex AND radial bonds are used for force dipole and inner network is difference from outer
    elif radial_only_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # only radial bonds are used for force dipole and inner network is difference from outer
    elif inner_radial_arp_flag_bash == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # only radial bonds are used for force dipole and inner network is difference from outer
    elif cluster_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # only radial bonds are used for force dipole and inner network is difference from outer
    elif cluster_new_flag == 1:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # only radial bonds are used for force dipole and inner network is difference from outer
    elif cluster_radial_flag == 1:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # only radial bonds are used for force dipole and inner network is difference from outer
    elif hex_flag == 1:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
    elif hex_rand_flag == 1:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
    elif force_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_force/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"   # only radial bonds are used for force dipole and inner network is difference from outer
    elif test_hex == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
    elif test_hex_rand == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
    elif anisotropic_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_arp_anisotropic/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
    else:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'+ kappa_fname +ranseed+"/"
        
        
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
    
    
all_data = []
strain_all_data = []
rlen_all_data = []

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


#all_data = np.array(all_data)

all_data = np.array(all_data)
xpos = all_data[:,:,1]
ypos = all_data[:,:,2]

rlen_all_data = np.array(rlen_all_data)
rlen_all_data = rlen_all_data[:,:,:3]

rlen_data_val = np.reshape(rlen_all_data,(num,3*L*L))

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
# if L != 64 and pbond != 1:
#     fname_connect = base+folder+'txt/strain/Lattice_connect_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
# if L != 64 and pbond == 1:
#     fname_connect = base+folder+'txt/strain/Lattice_connect_'+ "srand_" + str(srand) + '_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
# if srand_flag == 1 and pbond == 1 and L ==64:
#     fname_connect = base+folder+'txt/strain/Lattice_connect_'+ "srand_" + str(srand) + '_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
# if force_flag == 1:
#     fname_connect = base+folder+'txt/strain/Lattice_connect_'+ "srand_" + str(srand) + '_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
                
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

# boundary_nodes_circle = []
# for i in range(0,num_pts):
#     dx = xpos[0][center] - xpos[0][i]
#     dy = ypos[0][center] - ypos[0][i]
#     dist_center = np.sqrt(dx*dx + dy*dy)
#     if (dist_center >= (radius-dr)) and (dist_center < (radius + dr)):
#         boundary_nodes_circle.append(i)
# #    if (dist_center > radius):
# #        boundary_nodes_circle.append(i)



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

remove_loc_bndry = []
for i in range(0,len(outer_nodes)):
    if (outer_nodes[i] in boundary_nodes_circle):
        remove_loc_bndry.append(i)
        
outer_nodes = np.delete(outer_nodes, remove_loc_bndry)    # remving the boundary nodes from outer boundary array
#==============================================================================
# setting the node and line locations for plotting
#==============================================================================

# define the x-location and y-location array below to start using linecollection
xpos = np.array(xpos)
ypos = np.array(ypos)
x_loc_arr = []
y_loc_arr = []
x_loc_plot = []
y_loc_plot = []
boundary_nodes = []
node_bond_plot = []

unconnected_const = 9999
if L == 128: unconnected_const = -1

for k in range(0,num):
    boundary_nodes = []

    #for k in range(0,num-6):
    x_loc_arr.append([])
    y_loc_arr.append([])
    x_loc_plot.append([])
    y_loc_plot.append([])
    node_bond_plot.append([])
    for i in range(0,num_pts):
        # finding boundary nodes
        row = int(i/L)
        colmn = i%L            
        # populating the list of boundary nodes of the box. This index will be later used for contraction quantification            
        if row == 0 or row == (L-1):
            boundary_nodes.append(i)
        if colmn == 0 or colmn == (L-1):
            boundary_nodes.append(i)                    
        
#        for j in range(0,len(conn_node[0])):            
        for j in range(0,3):   
            if conn_node[i,j] != unconnected_const:
                if pbd_flag == 1:
                    if colmn != 0 and colmn != L-1:
                        xval = [xpos[k,i],xpos[k,conn_node[i,j]]]
                        x_loc_arr[k] = np.hstack((x_loc_arr[k],xval))
                    # pbd along x-axis
                    if colmn == 0:
                        if xpos[k,conn_node[i,j]] - xpos[k,i] > L/2:
                            x_pbd = xpos[k,conn_node[i,j]] - L
                            xval = [xpos[k,i],x_pbd]
                        else:
                            xval = [xpos[k,i],xpos[k,conn_node[i,j]]]
                        x_loc_arr[k] = np.hstack((x_loc_arr[k],xval))
                    if colmn == L-1:
                        if (xpos[k,i] - xpos[k,conn_node[i,j]]) > L/2:
                            x_pbd = xpos[k,conn_node[i,j]] + L
                            xval = [xpos[k,i],x_pbd]
                        else:
                            xval = [xpos[k,i],xpos[k,conn_node[i,j]]]
                        
                        x_loc_arr[k] = np.hstack((x_loc_arr[k],xval))
    
                    yval = [ypos[k,i],ypos[k,conn_node[i,j]]]
                    y_loc_arr[k] = np.hstack((y_loc_arr[k],yval))
                else:    
                    xval = [xpos[k,i],xpos[k,conn_node[i,j]]]
                    x_loc_arr[k] = np.hstack((x_loc_arr[k],xval))
                    yval = [ypos[k,i],ypos[k,conn_node[i,j]]]
                    y_loc_arr[k] = np.hstack((y_loc_arr[k],yval))
                    
                    if (i not in outer_nodes) and (conn_node[i,j] not in outer_nodes):
                        x_loc_plot[k] = np.hstack((x_loc_plot[k],xval))
                        y_loc_plot[k] = np.hstack((y_loc_plot[k],yval))
                        
                        node_bond_plot[k].append([i,j])
                    
node_bond_plot = np.array(node_bond_plot)
# finding the list of inner and ring nodes including the boundary nodes
all_nodes = np.arange(0,num_pts)
circle_nodes = np.delete(all_nodes, outer_nodes)   # these nodes are inside the outer boundary and on the outer boundary
# ###############################################################################    
# # setting up the strain array
# ###############################################################################
# reading the strain data
#if move_flag == 1:
#    strain_fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+move_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+"_force.txt"
#elif rlen_flag == 1:
#    strain_fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_"+str(num-1)+"_force.txt"
#print("strain file name: ", strain_fname)
#strain_data = np.loadtxt(strain_fname)
#strain = strain_data[:,1:4]
#strain = np.reshape(strain,(3*L*L))
strain_all_data = np.array(strain_all_data)
bond_len = strain_all_data[:,:,1:4]

bond_len = np.reshape(bond_len,(num,3*L*L))
#bond_len = np.reshape(bond_len,(num-6,3*L*L))

rlen_circle_data = rlen_all_data[:,circle_nodes,:3]
rlen_data_val_circle = np.reshape(rlen_circle_data,(num,3*len(circle_nodes)))
bond_len_circle = strain_all_data[:,circle_nodes,1:4]     # bond length of bonds associated with nodes in the outer circle, including the boundary
# bond_len_circle = np.reshape(bond_len_circle,(num,3*len(circle_nodes)))

strain_circle = []
bond_len_circle_flattened = []
rlen_circle_flattened = []

for i in range(0,len(bond_len_circle)):
    strain_circle.append([])
    bond_len_circle_flattened.append([])
    rlen_circle_flattened.append([])
    for j in range(0,len(circle_nodes)):
        for k in range(0,np.shape(bond_len_circle)[2]):
            if (conn_node[circle_nodes[j],k] not in outer_nodes):
                bond_len_circle_flattened[i].append(bond_len_circle[i,j,k])
                rlen_circle_flattened[i].append(rlen_circle_data[i,j,k])
                
                if bond_len_circle[i,j,k] != 0:
                    strain_circle[i].append(bond_len_circle[i,j,k] - rlen_circle_data[i,j,k])
                    
        

#strain = strain - 1.0    # subtracting restlength from bond length to get strain
strain = []
for i in range(0,len(bond_len)):
#    del strain[i][rm_loc]
    strain.append([])
    for j in range(0,len(bond_len[0])):
        if bond_len[i,j] != 0:
#            strain[i].append(bond_len[i,j] - 1.0)
            strain[i].append(bond_len[i,j] - rlen_data_val[i,j])


    
# for i in range(0,len(bond_len_circle)):
#     strain_circle.append([])
#     for j in range(0,len(bond_len_circle_flattened[0])):
#         if bond_len_circle_flattened[i,j] != 0:
# #            strain[i].append(bond_len[i,j] - 1.0)
#             strain_circle[i].append(bond_len_circle[i,j] - rlen_data_val_circle[i,j])

#strain = temp_strain
#strain[(strain != 0)] -= 1.0

strain = np.array(strain)
strain_circle = np.array(strain_circle)

# just for now a hack:: REMOVE LATER !!!!!
strain[0][:] = 0.0
strain_circle[0][:] = 0.0

###############################################################################
# plotting strain histograms
hist_base = base + folder + "png/strain_hist/"

# for i in range(0,len(strain)):
#     fig = plt.figure(figsize=(8, 4), dpi=300)
#     #mean_pos_abhi_str = 'Mean = ' + "%.2e" % np.mean(mean_pos_abhi)
#     #max_lim = 0.02      # works for p = 0.55
#     max_lim = np.max(strain)
#     min_lim = np.min(strain)
#     #binwidth = (max_lim-0)/100.
#     binwidth = 0.005
#     binwidth_str = "%.2e" % binwidth
#     bins = np.arange(min_lim-binwidth,max_lim+binwidth,binwidth)
    
#     n,bins,patches = plt.hist(strain[i], bins = bins, histtype = 'step', alpha = 0.5)#, color='blue')#, label = "CG " + mean_pos_abhi_str)
#     plt.yscale('symlog')
#     #plt.title('Positive Strains')
#     #plt.legend()
#     plt.xlabel('Strains')
#     plt.ylabel('Number of Bonds')
    
#     plt.subplots_adjust(top=0.92, left = 0.08, bottom = 0.12, right = 0.96)
    
# #    plt.savefig(hist_base+"strain_hist_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+binwidth_str+"_"+str(i)+"_"+str(num-1)+".png")  # just plotting the last step data
#     plt.close()
#     plt.clf()

# ###############################################################################    
# # setting the circular boundary
# ###############################################################################
# # making a circular region to identify the boundary nodes so that I can fix them in the simulation
# radius = 25
# dr = 0.5
# center = 2080

# if L == 128:
#     center = 8319
#     radius = 50

# # boundary_nodes_circle = []
# # for i in range(0,num_pts):
# #     dx = xpos[0][center] - xpos[0][i]
# #     dy = ypos[0][center] - ypos[0][i]
# #     dist_center = np.sqrt(dx*dx + dy*dy)
# #     if (dist_center >= (radius-dr)) and (dist_center < (radius + dr)):
# #         boundary_nodes_circle.append(i)
# # #    if (dist_center > radius):
# # #        boundary_nodes_circle.append(i)



# # reading the nodes on the circular boundary
# fname_boundary = base+folder+'txt/area/Boundary_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if L != 64:
#     fname_boundary = base+folder+'txt/area/Boundary_nodes_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# # if srand_flag == 1 and num_center > 1 and pbond == 1:
# if srand_flag == 1 and pbond == 1:
#     fname_boundary = base+folder+'txt/area/Boundary_nodes_'+ "srand_" + str(srand) + '_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if force_flag == 1:
#     fname_boundary = base+folder+'txt/area/Boundary_nodes_'+ "srand_" + str(srand) + '_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
    
        
# print("circular boundary filename: ",fname_boundary)
# boundary_data = np.loadtxt(fname_boundary)
# boundary_nodes_circle = boundary_data[:,1]
# boundary_nodes_circle = boundary_nodes_circle.astype(int)

# # manually checking the nodes that should not be send to the energy minimizer
# # outer_nodes_circle = []
# # for i in range(0,num_pts):
# #     dx = xpos[0][center] - xpos[0][i]
# #     dy = ypos[0][center] - ypos[0][i]
# #     dist_center = np.sqrt(dx*dx + dy*dy)
# #     if (dist_center > radius+dr):
# #         outer_nodes_circle.append(i)

# # reading the nodes outside the inner region :: this includes the circular boundary nodes btw
# fname_outer = base+folder+'txt/area/Outer_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if L != 64:
#     fname_outer = base+folder+'txt/area/Outer_nodes_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# # if srand_flag == 1 and num_center > 1 and pbond == 1:
# if srand_flag == 1 and pbond == 1:
#     fname_outer = base+folder+'txt/area/Outer_nodes_'+ "srand_" + str(srand) +'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if force_flag == 1:    
#     fname_outer = base+folder+'txt/area/Outer_nodes_'+ "srand_" + str(srand) +'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
        
# print("outer nodes filename: ",fname_outer)
# outer_data = np.loadtxt(fname_outer)
# #outer_nodes = outer_data
# outer_nodes = outer_data.astype(int)

# # reading the nodes in the inner region :: 
# fname_inner = base+folder+'txt/area/Inner_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if L != 64:
#     fname_inner = base+folder+'txt/area/Inner_nodes_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# # if srand_flag == 1 and num_center > 1 and pbond == 1:
# if srand_flag == 1 and pbond == 1:
#     fname_inner = base+folder+'txt/area/Inner_nodes_'+ "srand_" + str(srand) +'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if force_flag == 1:
#     fname_inner = base+folder+'txt/area/Inner_nodes_'+ "srand_" + str(srand) +'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+force_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
    
        
# if srand_flag == 1 and diff_flag == 1:
#     print("inner nodes filename: ",fname_inner)
#     inner_data = np.loadtxt(fname_inner)
#     #outer_nodes = outer_data
#     inner_nodes = inner_data.astype(int)
#     inner_nodes = inner_nodes[:,1]


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
strain_bond_cross_y_list =[]
bond_mid_x_list = []
bond_mid_y_list = []
radial_force_cross_bond_list = []
for k in range(0,np.shape(node_bond_plot)[0]):      # number of steps
    test = []
    strain_bond_inner_circle.append([])
    node_cross.append([])
    bond_cross.append([])
    bond_cross_index.append([])
    angle_cross_list.append([])
    strain_bond_cross_x_list.append([])
    strain_bond_cross_y_list.append([])
    bond_mid_x_list.append([])
    bond_mid_y_list.append([])
    radial_force_cross_bond_list.append([])
    for l in range(0,np.shape(node_bond_plot)[1]):     # number of bonds
        node = node_bond_plot[k,l,0]
        bond = node_bond_plot[k,l,1]
        connected_node = conn_node[node,bond]      # node connected to this node
                
        # distance from center
        dx_node = xpos[k,node] - xpos[k,center]
        dy_node = ypos[k,node] - ypos[k,center]
        dr_node = np.sqrt(dx_node**2 + dy_node**2)
        
        dx_connected_node = xpos[k,connected_node] - xpos[k,center]
        dy_connected_node = ypos[k,connected_node] - ypos[k,center]
        dr_connected_node = np.sqrt(dx_connected_node**2 + dy_connected_node**2)
        

        if ((dr_node < inner_radius) and (dr_connected_node > inner_radius) or (dr_node > inner_radius) and (dr_connected_node < inner_radius)):
            # multiplying by -1 because stretched bonds give negative stress and vice cersa as per our conventions
            strain_circle[k,l] = -1*strain_circle[k,l]
            strain_bond_inner_circle[k].append(strain_circle[k,l])        
            node_cross[k].append(node)
            bond_cross[k].append(bond)
            bond_cross_index[k].append(l)
                            
            #finding the angle of bond
            angle_bond = angle(xpos[k,node], ypos[k,node], xpos[k,connected_node], ypos[k,connected_node])
            angle_cross_list[k].append(angle_bond)
            # dividing the strain (force) in x and y components
            strain_bond_cross_x = strain_circle[k,l]*np.cos(angle_bond)
            strain_bond_cross_x_list[k].append(strain_bond_cross_x)
            strain_bond_cross_y = strain_circle[k,l]*np.sin(angle_bond)
            strain_bond_cross_y_list[k].append(strain_bond_cross_y)
            strain_bond_cross_vec = [strain_bond_cross_x, strain_bond_cross_y]
            # finding the midpoint of bond
            bond_mid_x = 0.5*(xpos[k,node]+xpos[k,connected_node])
            bond_mid_x_list[k].append(bond_mid_x)
            bond_mid_y = 0.5*(ypos[k,node]+ypos[k,connected_node])
            bond_mid_y_list[k].append(bond_mid_y)
            radial_bond_cross_vec = [bond_mid_x, bond_mid_y]

            radial_force_cross_bond = np.dot(strain_bond_cross_vec, np.transpose(radial_bond_cross_vec)) / np.linalg.norm(radial_bond_cross_vec)
            radial_force_cross_bond_list[k].append(radial_force_cross_bond)
            
            if node == 1307 and k == 10:
                print(node, bond, angle_bond, strain_circle[k,l], strain_bond_cross_vec, radial_bond_cross_vec, radial_force_cross_bond, np.linalg.norm(radial_bond_cross_vec))

total_radial_stress = []
for i in range(0,len(radial_force_cross_bond_list)):
    radial_stress = np.sum(radial_force_cross_bond_list[i])/(2*np.pi*inner_radius)
    total_radial_stress.append(radial_stress)
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
    
    
plot_folder = base + folder + "png/strain/"

movie_flag = 0
color_flag = 0

# for zoomed in movie:
if movie_flag == 1:
    lw = 2.5 # line width of bonds in plot

row_diff = 7
col_diff = -6
test_dip = [1108+L*row_diff+col_diff, 1107+L*row_diff+col_diff, 1043+L*row_diff+col_diff, 1044+L*row_diff+col_diff, 1109+L*row_diff+col_diff, 1172+L*row_diff+col_diff, 1171+L*row_diff+col_diff ]
test_dip1 = [1550,1549,1486,1487,1551,1615,1614]
test_dip2 = [1806,1805,1742,1743,1807,1871,1870]


# for buckling lets convert strain to force and plot that ::
if buckle_flag == 1:
    for i in range(0,len(strain)):
        for j in range(0,len(strain[0])):
            if strain[i,j] < 0:
                test = strain[i,j]
                strain[i,j] = strain[i,j]*mu_c   # converting strain to force
#                print(test, " ** ",strain[i,j])


#==============================================================================
# testing if we can produce just the circular region as a plot
#==============================================================================

dx = xpos[0] - xpos[0,center]
dy = ypos[0] - ypos[0,center]
dist_all = np.sqrt(dx*dx + dy*dy)
ring_rad = 12
if L == 128:
    ring_rad = 25
nodes_in_ring = np.squeeze( np.where( (dist_all > (ring_rad-0.5)) & (dist_all < (ring_rad+0.5)) ) )

test_flag = 1
if test_flag == 1:
    # for i in np.arange(0,len(strain_circle)):
    for i in np.arange(len(strain_circle)-1,len(strain_circle)):
#    for i in np.arange(0,2):
    # for i in np.arange(0,1):
        #plt_thrshld = 5.7e-5       # threshold to plot red adn blue colors for expanding and contracting bonds
#        plt_thrshld = np.sqrt(tol/(3*L*L))       # threshold to plot red adn blue colors for expanding and contracting bonds
#        plt_thrshld = np.max(strain[len(strain)-1])/100.0       # classic case: threshold to plot red adn blue colors for expanding and contracting bonds
#        plt_thrshld = np.max(strain[len(strain)-1])/50.0       # classic case: threshold to plot red adn blue colors for expanding and contracting bonds

        plt_thrshld = abs(strain_thresh)       # threshold to plot red adn blue colors for expanding and contracting bonds
#        plt_thrshld = (5)*np.mean(np.abs(strain[i]))       # threshold to plot red adn blue colors for expanding and contracting bonds
#        plt_thrshld = (10)*np.median(np.abs(strain[i]))       # threshold to plot red adn blue colors for expanding and contracting bonds

#        plt_thrshld = 0.0       # threshold to plot red adn blue colors for expanding and contracting bonds
        plt_thrshld_str = "%.2e" % plt_thrshld 
                
        # x = x_loc_arr[i]
        # y = y_loc_arr[i]
        x = x_loc_plot[i]
        y = y_loc_plot[i]
        
        # initial locations
        # x_init = x_loc_arr[0]
        # y_init = y_loc_arr[0]
        x_init = x_loc_plot[0]
        y_init = y_loc_plot[0]

#        x_ma = np.ma.masked_array(x_arr[i], mask=m_arr)
#        y_ma = np.ma.masked_array(y_arr[i], mask=m_arr)
        
        xy = []

        for j in np.arange(0,len(x)-1,2):
            xy.append([(x[j], y[j]),(x[j+1],y[j+1])])
            

#        for j in range(0,len(xy)):
#            if xy[j][0][0] == x_arr[10][2017*6] and xy[j][1][0] == x_arr[10][2016*6]:
#                print(xy[j][0][0], xy[j][1][0])
#                remove_index = j
#        xy.remove(xy[j])
#        xy.pop(remove_index)

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
#        col = np.array(['k']*len(strain))
#        coltest = np.array([(0.0,0.0,0.0,1.0)]*len(strain))

        lwidth = np.array([lw]*len(strain_circle[i]))
#        lwidth = np.array([lw]*len(strain))

        lc_arr = list(lc_arr)
        lwidth = list(lwidth)
        remove_loc = np.atleast_1d(np.array(remove_loc))

#        for item in sorted(remove_loc, reverse=True):
#            del xy[item]            
#            del lc_arr[item]
#            del lwidth[item]

#            xy.pop(item)            
#            lc_arr.pop(item)
        lc_arr = np.array(lc_arr)
        for kk in range(0,len(lc_arr)):
            if lc_arr[kk] < 0:
                lc_arr[kk] = lc_arr[kk]*mu_c
        lwidth = np.array(lwidth)

#        stretch_index = np.where((strain[i]<plt_thrshld) & (strain[i]>0.0))
#        compress_index = np.where((strain[i]>(-1.0*plt_thrshld)) & (strain[i]<0.0))
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

        colors=np.full((num_pts), 'grey')
        sc_size=np.full((num_pts), 0.2)
        

        if srand_flag == 1 and diff_flag == 1:
            sc_size[inner_nodes] = 10.0
            colors[inner_nodes] = 'k'
            if L == 128:
                sc_size[inner_nodes] = 3.0
                colors[inner_nodes] = 'k'

        sc_size[dip_node_arr] = 50.0
        if movie_flag == 1:        
            sc_size[dip_node_arr] = 200.0              # for movie
        if L ==128:
            sc_size[dip_node_arr] = 10.0              # for large L            
        colors[dip_node_arr] = 'g'

        # sc_size[boundary_nodes] = 20.0
        # colors[boundary_nodes] = 'm'
        
        sc_size[nodes_in_ring] = 20.0
        colors[nodes_in_ring] = 'm'

        # sc_size[node_cross[i]] = 20.0
        # colors[node_cross[i]] = 'm'
        
        
#        sc_size[boundary_nodes_circle] = 20.0
#        colors[boundary_nodes_circle] = 'm'

        sc_size[outer_nodes] = 0.0
        colors[outer_nodes] = 'w'

        sc_size[boundary_nodes_circle] = 20.0
        colors[boundary_nodes_circle] = 'm'
        
        # sc_size[1206] = 20.0
        # colors[1206] = 'r'

        # sc_size[[908, 909,1226,1289, 1353]] = 20.0
        # colors[[908, 909,1226,1289, 1353]] = 'b'

        # for finding the difference in N = 1 and N = 2 networks; unmatch is actually found at the end of this code
        # sc_size[np.unique(unmatch[0])] = 20.0
        # colors[np.unique(unmatch[0])] = 'r'

        # sc_size[1307] = 20.0
        # colors[1307] = 'b'

        if force_flag == 1:
            sc_size[ring_nodes] = 20.0
            colors[ring_nodes] = 'b'
            
        
        # sc_size[10000] = 50
        # colors[10000] = 'r'
        

#        special_nodes = [1963]
#        sc_size[special_nodes] = 40.0
#        colors[special_nodes] = 'r'

        # special_nodes1 = [2488]
        # sc_size[special_nodes1] = 40.0
        # colors[special_nodes1] = 'b'

        # testing by hand the missing nodes on the boundary
        # sc_size[2080+26*L-13] = 20.0
        #colors[2080+26*L-13] = 'g'

        # plotting special nodes separately where we moved the 20th dipole by hand to test transition in stretching energy
        # sc_size[test_dip1] = 20.0
        # colors[test_dip1] = 'r'
        # sc_size[test_dip2] = 20.0
        # colors[test_dip2] = 'b'

        col[stretch_index] = 'b'
        col[compress_index] = 'r'
        col[zero_index] = 'k'
        lwidth[stretch_index] = 4*lw
        lwidth[compress_index] = 4*lw
        
        # col[bond_cross_index[i]] = 'g'
        # lwidth[bond_cross_index[i]] = 4*lw
        # lwidth[bond_cross_index[i]] = 4*lw
        
        # sc_size[223] = 20.0
        # colors[223] = 'b'
        # sc_size[226] = 20.0
        # colors[226] = 'b'
                
#        lc = LineCollection(xy, lw = lwidth, cmap = newcmp, norm=plt.Normalize(vmin=-vlim,vmax=vlim))
        if color_flag == 1:
#            lc = LineCollection(xy, lw = lwidth, cmap = newcmp)#, norm=MidpointNormalize(midpoint=0.))#, alpha = 0.5)#,vmin=-0.05))
#            lc = LineCollection(xy, lw = lwidth, cmap = newcmp, norm=plt.Normalize(vmin=-0.05,vmax=0.05))#, alpha = 0.5)#,vmin=-0.05))
#            lc = LineCollection(xy, lw = lwidth, array = lc_arr, cmap = newcmp, norm=plt.Normalize(vmin=strain_neg_max,vmax=strain_pos_max))#, alpha = 0.5)#,vmin=-0.05))
            lc = LineCollection(xy, lw = lwidth, array = lc_arr, cmap = newcmp, norm=plt.Normalize(vmin=-vlim,vmax=vlim))#, alpha = 0.5)#,vmin=-0.05))
#            lc = LineCollection(xy, lw = lwidth, array = lc_arr, cmap = newcmp, norm=MidpointNormalize(midpoint=0.))
#            lc = LineCollection(xy, lw = lwidth, colors = 'k', alpha = alpha)#,vmin=-0.05))
#            lc = LineCollection(xy, lw = lwidth, cmap = 'jet', norm=plt.Normalize(vmin=-vlim,vmax=vlim))
        else:
            lc = LineCollection(xy, lw = lwidth, colors = col)#, alpha = 0.5)#,vmin=-0.05))

#        lc = LineCollection(xy, lw = lwidth, color=coltest)
        
        fig = plt.figure(figsize=(15, 9), dpi=100)
        ax1 = fig.add_subplot(1, 1, 1)
#        line = ax1.add_collection(lc)

#        lc.set_array(lc_arr)

#        lc.set_alpha(abs(strain[i]))
# plotting the really high negative strains separately to keep the color scheme sane
#        lc_high_strain = LineCollection(high_strain[remove_loc], lw = 2.0, colors = 'k', linestyles = 'dotted')#,vmin=-0.05))

#        lc_test = LineCollection(high_strain[[remove_loc]], lw = 2.0, colors = 'k', linestyles = 'dotted')#, alpha = 0.5)#,vmin=-0.05))

        if color_flag == 1:
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes("right", size="5%", pad=0.05)
#            cbar = fig.colorbar(lc, ax = ax1, format=ticker.FuncFormatter(fmt), cax=cax, ticks=[strain_neg_max, 0, strain_pos_max])
    
            cbar = fig.colorbar(lc, ax = ax1, format=ticker.FuncFormatter(fmt), cax=cax)
    #        cbar.tick_pa
    #        cbar.tick_params(labelsize=10)
            for t in cbar.ax.get_yticklabels():
                    t.set_fontsize(30)
            tick_locator = ticker.MaxNLocator(nbins=3)
            cbar.locator = tick_locator
            cbar.update_ticks()        

#        cbar.ax.locator_params(nbins=6)
        ax1.add_collection(lc)
#        ax1.add_collection(lc_high_strain)

#        ax1.add_collection(lc_test)
        
        # initial positions
        
#        ax1.autoscale()
#        ax1.scatter(all_data[i][:,1], all_data[i][:,2], s = sc_size, marker = "o", color = colors, alpha = 0.5)
#        ax1.scatter(xpos, ypos, s = sc_size, marker = "o", color = colors, alpha = 0.5)

#        sc_size[[3548,3547,3483,3484,3549,3612,3611]] = 20
#        colors[[3548,3547,3483,3484,3549,3612,3611]] = 'r'
        
        ax1.scatter(xpos[i], ypos[i], s = sc_size, marker = "o", color = colors, alpha = 0.5)

        # plotting boundary nodes of undeformed lattice for comparison
#         if (i > 0):
# #            ax1.scatter(all_data[0][boundary_nodes,1], all_data[0][boundary_nodes,2], s = 20, marker = "o", color = 'k', alpha = 0.3)
# #            ax1.scatter(xpos[boundary_nodes], ypos[boundary_nodes], s = 20, marker = "o", color = 'k', alpha = 0.3)
#             ax1.scatter(xpos[0,boundary_nodes], ypos[0,boundary_nodes], s = 20, marker = "o", color = 'k', alpha = 0.3)
        
        # plotting an inner circle:
        # ax1.plot(x_inner_circle[i], y_inner_circle[i])
        # ax1.plot(x_inner_circle_13[i], y_inner_circle_13[i], color = 'r')
        
        # zooming in for movie generation:
        if movie_flag == 1:
            ax1.set_xlim([30,50])
            ax1.set_ylim([35,50])

        ax1.set_aspect('equal')
#        ax1.set_title("Strain", fontsize=16)

        ax1.tick_params(axis='both', which='major', labelsize = 18, left = False, bottom = False)
        ax1.tick_params(axis='both', which='minor', labelsize = 18,  left = False, bottom = False)
        
        ax1.set_xticks([])
        ax1.set_yticks([])

#        ax1.tick_params(left = False, right = False , labelleft = False, labelbottom = False, bottom = False)        
        plt.subplots_adjust(right=0.995, left = 0.0, top = 0.9, bottom = 0.05)
        plt.subplots_adjust(right=0.85,left=0.05,top=0.95,bottom=0.05)

        plot_fname = plot_folder + 'circle_binary_strain_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
        if srand_flag == 1 and pbond == 1:
            plot_fname = plot_folder + 'circle_binary_strain_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if L ==128:
            plot_fname = plot_folder + 'circle_binary_strain_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            # if srand_flag == 1 and num_center > 1 and pbond == 1:
            if srand_flag == 1:
                plot_fname = plot_folder + 'circle_binary_strain_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if color_flag == 1:
            plot_fname = plot_folder + 'circle_colored_strain_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if L != 64:
                plot_fname = plot_folder + 'circle_colored_strain_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if srand_flag == 1:
                plot_fname = plot_folder + 'circle_colored_strain_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
                
        if force_flag == 1:
            plot_fname = plot_folder + 'circle_binary_strain_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+force_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            

        if movie_flag != 1:
            plt.savefig(plot_fname, dpi = 300)
                
#            for movie generation
        if movie_flag == 1:
            fnum = '%03d' % i
            plot_fname = plot_folder + 'img'+fnum+'.jpg'            
            plt.savefig(plot_fname, dpi = 100)

        print("plotted to: ", plot_fname)
        plt.clf()
        plt.close('all')


#==============================================================================
# plotting the bending energy heat map on top of the scatter plot
# it superimposes the bending enregy from last step to all the steps
#==============================================================================
from scipy.interpolate import griddata
from matplotlib.colors import Normalize

bend_heat_flag = 1  

if bend_heat_flag == 1:
    # for i in np.arange(0,len(strain_circle)):
    for i in np.arange(len(strain_circle)-1,len(strain_circle)):

        plt_thrshld = abs(strain_thresh)       # threshold to plot red adn blue colors for expanding and contracting bonds
        plt_thrshld_str = "%.2e" % plt_thrshld 
                
        x = x_loc_plot[i]
        y = y_loc_plot[i]

        x_init = x_loc_plot[0]
        y_init = y_loc_plot[0]
        
        xy = []

        for j in np.arange(0,len(x)-1,2):
            xy.append([(x[j], y[j]),(x[j+1],y[j+1])])
            
        # Define grid for heatmap
        xi = np.linspace(min(xpos[-1]), max(xpos[-1]), 300)
        yi = np.linspace(min(ypos[-1]), max(ypos[-1]), 300)
        xi, yi = np.meshgrid(xi, yi)
        # Interpolate scattered data to grid
        zi = griddata((xpos[-1], ypos[-1]), bend_en_node, (xi, yi), method='cubic')


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
#        col = np.array(['k']*len(strain))
#        coltest = np.array([(0.0,0.0,0.0,1.0)]*len(strain))

        lwidth = np.array([lw]*len(strain_circle[i]))
#        lwidth = np.array([lw]*len(strain))

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

        colors=np.full((num_pts), 'grey')
        sc_size=np.full((num_pts), 0.2)
        
        if srand_flag == 1 and diff_flag == 1:
            sc_size[inner_nodes] = 10.0
            colors[inner_nodes] = 'k'
            if L == 128:
                sc_size[inner_nodes] = 3.0
                colors[inner_nodes] = 'k'

        sc_size[dip_node_arr] = 50.0
        if movie_flag == 1:        
            sc_size[dip_node_arr] = 200.0              # for movie
        if L ==128:
            sc_size[dip_node_arr] = 10.0              # for large L            
        colors[dip_node_arr] = 'g'
        
        sc_size[nodes_in_ring] = 20.0
        colors[nodes_in_ring] = 'm'

        sc_size[outer_nodes] = 0.0
        colors[outer_nodes] = 'w'

        sc_size[boundary_nodes_circle] = 20.0
        colors[boundary_nodes_circle] = 'm'
        
        if force_flag == 1:
            sc_size[ring_nodes] = 20.0
            colors[ring_nodes] = 'b'
            
        col[stretch_index] = 'b'
        col[compress_index] = 'r'
        col[zero_index] = 'k'
        lwidth[stretch_index] = 4*lw
        lwidth[compress_index] = 4*lw
        
                
#        lc = LineCollection(xy, lw = lwidth, cmap = newcmp, norm=plt.Normalize(vmin=-vlim,vmax=vlim))
        if color_flag == 1:
#            lc = LineCollection(xy, lw = lwidth, cmap = newcmp)#, norm=MidpointNormalize(midpoint=0.))#, alpha = 0.5)#,vmin=-0.05))
#            lc = LineCollection(xy, lw = lwidth, cmap = newcmp, norm=plt.Normalize(vmin=-0.05,vmax=0.05))#, alpha = 0.5)#,vmin=-0.05))
#            lc = LineCollection(xy, lw = lwidth, array = lc_arr, cmap = newcmp, norm=plt.Normalize(vmin=strain_neg_max,vmax=strain_pos_max))#, alpha = 0.5)#,vmin=-0.05))
            lc = LineCollection(xy, lw = lwidth, array = lc_arr, cmap = newcmp, norm=plt.Normalize(vmin=-vlim,vmax=vlim))#, alpha = 0.5)#,vmin=-0.05))
#            lc = LineCollection(xy, lw = lwidth, array = lc_arr, cmap = newcmp, norm=MidpointNormalize(midpoint=0.))
#            lc = LineCollection(xy, lw = lwidth, colors = 'k', alpha = alpha)#,vmin=-0.05))
#            lc = LineCollection(xy, lw = lwidth, cmap = 'jet', norm=plt.Normalize(vmin=-vlim,vmax=vlim))
        else:
            lc = LineCollection(xy, lw = lwidth, colors = col)#, alpha = 0.5)#,vmin=-0.05))

#        lc = LineCollection(xy, lw = lwidth, color=coltest)
        
        fig = plt.figure(figsize=(15, 9), dpi=100)
        ax1 = fig.add_subplot(1, 1, 1)

        colors_cmap = [(1, 1, 1), (1, 0.5, 0)]  # white to orange
        cmap_orange = LinearSegmentedColormap.from_list("fade_white_orange", colors_cmap)

        # max_val = 5e-9  # threshold for maximum orange        
        # # Clip zi so all values > max_val are set to max_val
        # zi_clipped = np.clip(zi, a_min=None, a_max=max_val)
        # # Normalize between min and max_val
        # norm = Normalize(vmin=0, vmax=max_val)

        plt.contourf(xi, yi, zi, levels=100, cmap=cmap_orange, vmin=0)  # heatmap

        cbar = plt.colorbar()
        cbar.set_label('Bending Energy', fontsize=25, labelpad=10)  # set font size and padding
        cbar.ax.tick_params(labelsize=20)  # Set tick label font size to 12


        if color_flag == 1:
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes("right", size="5%", pad=0.05)
#            cbar = fig.colorbar(lc, ax = ax1, format=ticker.FuncFormatter(fmt), cax=cax, ticks=[strain_neg_max, 0, strain_pos_max])
    
            cbar = fig.colorbar(lc, ax = ax1, format=ticker.FuncFormatter(fmt), cax=cax)
    #        cbar.tick_pa
    #        cbar.tick_params(labelsize=10)
            for t in cbar.ax.get_yticklabels():
                    t.set_fontsize(30)
            tick_locator = ticker.MaxNLocator(nbins=3)
            cbar.locator = tick_locator
            cbar.update_ticks()        

        ax1.add_collection(lc)
        
        ax1.scatter(xpos[i], ypos[i], s = sc_size, marker = "o", color = colors, alpha = 0.5)
        
        ax1.set_aspect('equal')

        ax1.tick_params(axis='both', which='major', labelsize = 18, left = False, bottom = False)
        ax1.tick_params(axis='both', which='minor', labelsize = 18,  left = False, bottom = False)
        
        ax1.set_xticks([])
        ax1.set_yticks([])

#        ax1.tick_params(left = False, right = False , labelleft = False, labelbottom = False, bottom = False)        
        plt.subplots_adjust(right=0.995, left = 0.0, top = 0.9, bottom = 0.05)
        plt.subplots_adjust(right=0.85,left=0.05,top=0.95,bottom=0.05)

        plot_fname = plot_folder + 'circle_heatmap_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if srand_flag == 1 and pbond == 1:
            plot_fname = plot_folder + 'circle_heatmap_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if L ==128:
            plot_fname = plot_folder + 'circle_heatmap_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            # if srand_flag == 1 and num_center > 1 and pbond == 1:
            if srand_flag == 1:
                plot_fname = plot_folder + 'circle_heatmap_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if color_flag == 1:
            plot_fname = plot_folder + 'circle_heatmap_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if L != 64:
                plot_fname = plot_folder + 'circle_heatmap_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if srand_flag == 1:
                plot_fname = plot_folder + 'circle_heatmap_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
                
        if force_flag == 1:
            plot_fname = plot_folder + 'circle_heatmap_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+force_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            
        plt.savefig(plot_fname, dpi = 300)
                
        print("plotted to: ", plot_fname)
        plt.clf()
        plt.close('all')


#==============================================================================
#--------- Forces on the boundary nodes ---------------------------------------
#==============================================================================

# reading force on each boundary node
force_fname = base+folder+"txt/area/"+"bndry_node_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    force_fname = base+folder+"txt/area/"+"bndry_node_force_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    force_fname = base+folder+"txt/area/"+"bndry_node_force_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    force_fname = base+folder+"txt/area/"+"bndry_node_force_" +  'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

#if diff_flag == 1:
#    force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('here all individual forces on BOUNDARY nodes are read from file : ',force_fname)    
force_data = np.loadtxt(force_fname)
perp_force_bndry = force_data[:,-1]
force_x_bndry = force_data[:,1]
force_y_bndry = force_data[:,2]

test_force_flag = 1
if test_force_flag == 1:
    # for i in np.arange(0,len(strain)):
    for i in np.arange(len(strain)-1,len(strain)):
    # for i in np.arange(0,1):

        plt_thrshld = abs(strain_thresh)       # threshold to plot red adn blue colors for expanding and contracting bonds

#        plt_thrshld = 0.0       # threshold to plot red adn blue colors for expanding and contracting bonds
        plt_thrshld_str = "%.2e" % plt_thrshld 
                
        # x = x_loc_arr[i]
        # y = y_loc_arr[i]
        
        # # initial locations
        # x_init = x_loc_arr[0]
        # y_init = y_loc_arr[0]

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
        # remove_loc = np.squeeze(np.where(strain[i] < 0.5*np.min(strain[i])))
        remove_loc = np.squeeze(np.where(strain_circle[i] < 0.5*np.min(strain_circle[i])))
#        remove_loc = np.squeeze(np.where(strain < 0.5*np.min(strain)))

        # lc_arr = np.zeros_like(strain[i])
        # lc_arr[:] = strain[i][:]
        lc_arr = np.zeros_like(strain_circle[i])
#        lc_arr = np.zeros_like(strain)
        lc_arr[:] = strain_circle[i][:]
#        lc_arr[:] = strain[:]

        if log_flag == 1:
            for k in range(len(lc_arr)):
                if lc_arr[k] > 0: lc_arr[k] = -1*np.log(lc_arr[k])
                if lc_arr[k] < 0: lc_arr[k] = np.log(abs(lc_arr[k]))
        
        col = np.array(['k']*len(strain_circle[0]))
        coltest = np.array([(0.0,0.0,0.0,1.0)]*len(strain_circle[0]))

        # col = np.array(['k']*len(strain[0]))
        # coltest = np.array([(0.0,0.0,0.0,1.0)]*len(strain[0]))

        lwidth = np.array([lw]*len(strain_circle[i]))
        # lwidth = np.array([lw]*len(strain[i]))

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
        # stretch_index = np.where(strain[i]>plt_thrshld)
        # compress_index = np.where(strain[i]<(-1.0*plt_thrshld))

        zero_index = np.where(abs(lc_arr) < plt_thrshld)
        zero_index = np.where(abs(strain[i]) < plt_thrshld)
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

        colors=np.full((num_pts), 'grey')
        sc_size=np.full((num_pts), 0.2)
        

        if srand_flag == 1 and diff_flag == 1:
            sc_size[inner_nodes] = 10.0
            colors[inner_nodes] = 'k'
            if L == 128:
                sc_size[inner_nodes] = 3.0
                colors[inner_nodes] = 'k'

        sc_size[dip_node_arr] = 50.0
        if movie_flag == 1:        
            sc_size[dip_node_arr] = 200.0              # for movie
        if L ==128:
            sc_size[dip_node_arr] = 10.0              # for large L            
        colors[dip_node_arr] = 'g'

        # sc_size[boundary_nodes] = 20.0
        # colors[boundary_nodes] = 'm'
        
        # sc_size[nodes_in_ring] = 20.0
        # colors[nodes_in_ring] = 'm'
        
        sc_size[outer_nodes] = 0.0
        colors[outer_nodes] = 'w'

        sc_size[boundary_nodes_circle] = 20.0
        colors[boundary_nodes_circle] = 'm'
        
        if force_flag == 1:
            sc_size[ring_nodes] = 20.0
            colors[ring_nodes] = 'b'
            
        col[stretch_index] = 'b'
        col[compress_index] = 'r'
        col[zero_index] = 'k'
        lwidth[stretch_index] = 4*lw
        lwidth[compress_index] = 4*lw
        
                
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
        
        ax1.scatter(xpos[i], ypos[i], s = sc_size, marker = "o", color = colors, alpha = 0.5)#, zorder = 5)
        # just the dipole nodes separately
        ax1.scatter(xpos[i, dip_node_arr], ypos[i, dip_node_arr], s = 20, marker = "o", color = 'g', zorder = 50)

        # # plotting boundary nodes of undeformed lattice for comparison
        # if (i > 0):
        #     ax1.scatter(xpos[0,boundary_nodes], ypos[0,boundary_nodes], s = 20, marker = "o", color = 'k', alpha = 0.3)
        
        # zooming in for movie generation:
        if movie_flag == 1:
            ax1.set_xlim([30,50])
            ax1.set_ylim([35,50])

        ax1.set_aspect('equal')
#        ax1.set_title("Strain", fontsize=16)

        ax1.tick_params(axis='both', which='major', labelsize = 18, left = False, bottom = False)
        ax1.tick_params(axis='both', which='minor', labelsize = 18,  left = False, bottom = False)
        
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        
        # Step 1: Extract x and y positions from master list
        x_pos_quiver = xpos[i,boundary_nodes_circle]
        y_pos_quiver = ypos[i,boundary_nodes_circle]
        
        # Step 2: Compute centroid (for polar angle reference)
        x0 = np.mean(x_pos_quiver)
        y0 = np.mean(y_pos_quiver)
        
        # Step 3: Compute theta
        theta = np.arctan2(y_pos_quiver - x0, x_pos_quiver - y0)  # Note: correct order is (y - y0, x - x0)
        theta = np.mod(theta, 2*np.pi)              # optional: wrap to [0, 2π)
        
        # Step 4: Sort by theta
        theta_sort_idx = np.argsort(theta)
        
        # Apply sort to all relevant arrays
        index_sorted = boundary_nodes_circle[theta_sort_idx]
        x_sorted = x_pos_quiver[theta_sort_idx]
        y_sorted = y_pos_quiver[theta_sort_idx]
        fx_sorted = force_x_bndry[theta_sort_idx]
        fy_sorted = force_y_bndry[theta_sort_idx]
        
        # Step 5: Choose every 2nd node for plotting
        n = len(x_sorted)
        idxs_to_plot = np.arange(0, n, 2)
        
        # Step 6: For each index, get average of force with neighbors (with wraparound)
        fx_smoothed = []
        fy_smoothed = []
        
        for ii in idxs_to_plot:
            i_minus = (ii - 1) % n
            i_plus = (ii + 1) % n
            fx_avg = (fx_sorted[i_minus] + fx_sorted[ii] + fx_sorted[i_plus]) / 3
            fy_avg = (fy_sorted[i_minus] + fy_sorted[ii] + fy_sorted[i_plus]) / 3
            fx_smoothed.append(fx_avg)
            fy_smoothed.append(fy_avg)
            
        # Convert to arrays
        x_plot_quiver = x_sorted[idxs_to_plot]
        y_plot_quiver = y_sorted[idxs_to_plot]
        fx_smoothed_quiver = np.array(fx_smoothed)
        fy_smoothed_quiver = np.array(fy_smoothed)
        
        # quiver plot to set the force arrows on the plots
        if num_center == 5 and pbond == 0.55:
            quiver_scale = 2e-6
        elif num_center == 15 and pbond == 0.55:
            quiver_scale = 1e-5
        else:
            # quiver_scale = 5*1e-3
            quiver_scale = 10*np.mean(np.sqrt(fx_smoothed_quiver**2 + fy_smoothed_quiver**2))

        quiver_scale_txt = "%.4f" % quiver_scale

        # plt.quiver(xpos[i,boundary_nodes_circle], ypos[i,boundary_nodes_circle], force_x_bndry, force_y_bndry, scale = quiver_scale, zorder = 10)
        plt.quiver(x_plot_quiver, y_plot_quiver, fx_smoothed_quiver, fy_smoothed_quiver, scale = quiver_scale, zorder = 10)

#        ax1.tick_params(left = False, right = False , labelleft = False, labelbottom = False, bottom = False)        
        plt.subplots_adjust(right=0.995, left = 0.0, top = 0.9, bottom = 0.05)
        plt.subplots_adjust(right=0.85,left=0.05,top=0.95,bottom=0.05)

        plot_fname = plot_folder + 'force_bndry_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+quiver_scale_txt+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
        if srand_flag == 1 and pbond == 1:
            plot_fname = plot_folder + 'force_bndry_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if L ==128:
            plot_fname = plot_folder + 'force_bndry_'+str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if srand_flag == 1:
                plot_fname = plot_folder + 'force_bndry_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            

        if color_flag == 1:
            plot_fname = plot_folder + 'force_bndry_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if L != 64:
                plot_fname = plot_folder + 'force_bndry_'+ 'srand_'+ str(srand) + "_" + str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            if srand_flag == 1 and num_center > 1 and pbond == 1:
                plot_fname = plot_folder + 'force_bndry_'+ "srand_" + str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
                
        if force_flag == 1:
            plot_fname = plot_folder + 'force_bndry_'+ "srand_" + str(srand) + "_" +str(L)+'_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+force_txt+"_"+plt_thrshld_str+"_"+str(L)+"_"+str(i)+"_"+str(num-1)+'.png'            
            
        if movie_flag != 1:
            plt.savefig(plot_fname, dpi = 300)
                
#            for movie generation
        if movie_flag == 1:
            fnum = '%03d' % i
            plot_fname = plot_folder + 'img'+fnum+'.jpg'            
            plt.savefig(plot_fname, dpi = 100)

        print("plotted to: ", plot_fname)
        plt.clf()
        plt.close('all')



###############################################################################    
# making a ring to find the decay of radial stress vs radius
# should only do it when the dipole is in the center!!
###############################################################################
outer_radius = 25
radius_bins = np.arange(1,outer_radius+1)
radius_ring_nodes = []

x_center = xpos[0,center]
y_center = ypos[0,center]
drr = 0.5

for k in range(0,len(radius_bins)):
    radius_ring_nodes.append([])
    for i in range(0,num_pts):
        x1 = xpos[-1,i]
        y1 = ypos[-1,i]
    
        dx = x_center-x1
        dy = y_center-y1
    
        dist = np.sqrt(dx*dx + dy*dy)            
        
        if ((dist>radius_bins[k]-drr) and (dist<radius_bins[k]+drr)):
            radius_ring_nodes[k].append(i)


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
count_inner_bonds = 0   # stores the number of bonds in p<1 simulations
temp_bond = []
for i in range(0,len(inner_nodes)):
    for j in range(0,3):
        if (conn_node[inner_nodes[i],j] != unconnected_const):          # only works for L =< 64 case
            dx_center = xpos[0, center] - xpos[0, conn_node[inner_nodes[i],j]]
            dy_center = ypos[0, center] - ypos[0, conn_node[inner_nodes[i],j]]
            dist_center = np.sqrt(dx_center*dx_center + dy_center*dy_center)
    # 		// finding the number of bonds in inner network for a p = 1 network
            if (dist_center <= inner_radius):
                    count_inner_bonds = count_inner_bonds+1;
                    temp_bond.append([inner_nodes[i],j])
				
temp_bond = np.array(temp_bond)

#============================================================================================================================================
## following hsa been commented out on Apr 9 2025 because the data being written to file is not being used anywhere else

# if L == 64 and pbond < 1:
#     # read how many bonds are present in a p = 1 network for the corresponding srand
#     # bond_count_fname = base+"lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/txt/area/"+"inner_bond_count_srand_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_1.00_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
#     # if L != 64:
#     #     bond_count_fname = base+folder+"txt/area/"+"inner_bond_count_"+ 'srand_'+str(srand)+"_" + str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_1.00_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
#     # print('here the inner bond count for p=1 network is read from file : ',bond_count_fname)    
#     # bond_count_inner_p1 = np.loadtxt(bond_count_fname)
#     # bond_count_inner_p1 = bond_count_inner_p1[1]
#     # bond_ratio = count_inner_bonds/bond_count_inner_p1     # this is total number of bonds in p=1 network
#     bond_ratio = count_inner_bonds/1470     # this is total number of bonds in p=1 network for L = 64
# elif L == 64 and pbond == 1:
#     bond_ratio = 1.0    
    

# bond_count_outfname = base+folder+"txt/area/"+"inner_bond_count_srand_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if L != 64:
#     bond_count_outfname = base+folder+"txt/area/"+"inner_bond_count_"+ 'srand_'+str(srand)+"_" + str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# print('here the inner bond count and bond percentage is written to file : ',bond_count_outfname)    
# heading = 'Inner Nodes      Bond count      Bond ratio'
# fmt = '%12d', '%12d', '%20.5e'
# np.savetxt(bond_count_outfname, np.column_stack((len(inner_nodes), count_inner_bonds, bond_ratio)), header = heading, fmt = fmt)
#============================================================================================================================================

# # writing the inner nodes and bonds - both # this is an important check step because there is some discrepancy
# bond_count_outfname = base+folder+"txt/area/"+"test_inner_bond_count_srand_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if L != 64:
#     bond_count_outfname = base+folder+"txt/area/"+"test_inner_bond_count_"+ 'srand_'+str(srand)+"_" + str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# print('here the inner bond count and bond are written to file : ',bond_count_outfname)    
# heading = 'Inner Node      Bond'
# fmt = '%12d', '%12d'
# np.savetxt(bond_count_outfname, np.column_stack((temp_bond[:,0], temp_bond[:,1])), header = heading, fmt = fmt)

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
bins = np.arange(1,xlim,1)
dd = 0.5
index_list = []
binned_disp = []
binned_mean_disp = []
binned_median_disp = []
binned_std_disp = []
for i in range(0,len(bins)):
    lower_lim = bins[i] - dd
    upper_lim = bins[i] + dd
    disp_in_bin = init_radial_dist[(init_radial_dist >= lower_lim) & (init_radial_dist < upper_lim)]    # finding displacements in the range if lower to upper limit
    index_radial = np.squeeze( np.where((init_radial_dist >= lower_lim) & (init_radial_dist < upper_lim)) )
    # finding displacements in the range if lower to upper limit
    # index_radial = index_radial.astype(int)
    # print(i,lower_lim,upper_lim, len(index_radial))
    index_list.append(index_radial)
    binned_disp.append(radial_disp[index_radial])   # binned displacements
    # binned_disp.append(disp_in_bin)   # binned displacements
    mean_disp_bin = np.mean(radial_disp[index_radial])
    # mean_disp_bin = np.mean(disp_in_bin)
    median_disp_bin = np.median(radial_disp[index_radial])
    # median_disp_bin = np.median(disp_in_bin)
    std_disp_bin = np.std(radial_disp[index_radial])
    binned_mean_disp.append(mean_disp_bin)
    binned_median_disp.append(median_disp_bin)
    binned_std_disp.append(std_disp_bin)


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


###############################################################################
# comparing N = 1 and N= 2 networsk - where exactly are they different?
###############################################################################

# # reading the connection array for N = 1 case
# # fname_connect_N1 = base+folder+'txt/strain/Lattice_connect_'+str(1)+"_"+str(6)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name                
# folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/111111,111111/"   # only radial bonds are used for force dipole and inner network is difference from outer
# fname_connect_N1 = base+folder+'txt/strain/Lattice_connect_'+str(10)+"_"+str(60)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name                
# print("connection filename for N = 1: ",fname_connect_N1)
# conn_data_N1 = np.loadtxt(fname_connect_N1)
# conn_node_N1 = conn_data_N1[:,3:9]
# conn_node_N1 = (np.rint(conn_node_N1)).astype(int)

# # reading the connection array for N = 2 case
# # fname_connect_N2 = base+folder+'txt/strain/Lattice_connect_'+str(2)+"_"+str(12)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name                
# folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/111111,121111/"   # only radial bonds are used for force dipole and inner network is difference from outer
# fname_connect_N2 = base+folder+'txt/strain/Lattice_connect_'+str(10)+"_"+str(60)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name                
# print("connection filename for N = 2: ",fname_connect_N2)
# conn_data_N2 = np.loadtxt(fname_connect_N2)
# conn_node_N2 = conn_data_N2[:,3:9]
# conn_node_N2 = (np.rint(conn_node_N2)).astype(int)

# # getting the index of where they differ
# unmatch = np.where((conn_data_N1 == conn_data_N2) == False)

'''
# ###############################################################################    
# reading different cases of radial displacement to plot on the same plot for a figure panel for the paper
# ###############################################################################    
temp_num_center1 = 5
temp_num_center2 = 10
temp_num_center3 = 15

bin_slice_val = 13
xdata = bins[bin_slice_val:]
fit_x_outer = np.linspace(xdata[0], r2, 100)     # outer region x range
bins_outer = np.arange(xdata[0],r2+1)

# reading binned data
disp_mean_outfname1 = "/home/abhinav/david/cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/kappa2_e-6/667720,601210/112/means/bndry_node_radial_disp_mean_112_"+str(temp_num_center1)+"_"+str(temp_num_center1*6)+"_0.55_1.00e-07_1.00e-06_0.9000_1.0000_1.0000_10.txt"
print('File1: here all mean (after binning) radial displacements on all nodes is read from file : ',disp_mean_outfname1)    
radial_data_all1 = np.loadtxt(disp_mean_outfname1)
radial_data1 = radial_data_all1[bin_slice_val:,1]
radial_data1_std = radial_data_all1[bin_slice_val:,2]
popt1, pcov1 = curve_fit(func, xdata, radial_data1, maxfev=10000)
# label = '$%s x^{-%s}$' % tuple(f2str.format(t) for t in [popt[0], popt[1]])
u_r_outer_fit1 = func(fit_x_outer, popt1[0])

disp_mean_outfname2 = "/home/abhinav/david/cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/kappa2_e-6/667720,601210/112/means/bndry_node_radial_disp_mean_112_"+str(temp_num_center2)+"_"+str(temp_num_center2*6)+"_0.55_1.00e-07_1.00e-06_0.9000_1.0000_1.0000_10.txt"
print('File2: here all mean (after binning) radial displacements on all nodes is read from file : ',disp_mean_outfname2)    
radial_data_all2 = np.loadtxt(disp_mean_outfname2)
radial_data2 = radial_data_all2[bin_slice_val:,1]
radial_data2_std = radial_data_all2[bin_slice_val:,2]
popt2, pcov2 = curve_fit(func, xdata, radial_data2, maxfev=10000)
# label = '$%s x^{-%s}$' % tuple(f2str.format(t) for t in [popt[0], popt[1]])
u_r_outer_fit2 = func(fit_x_outer, popt2[0])
disp_mean_outfname3 = "/home/abhinav/david/cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/kappa2_e-6/667720,601210/112/means/bndry_node_radial_disp_mean_112_"+str(temp_num_center3)+"_"+str(temp_num_center3*6)+"_0.55_1.00e-07_1.00e-06_0.9000_1.0000_1.0000_10.txt"
print('File2: here all mean (after binning) radial displacements on all nodes is read from file : ',disp_mean_outfname3)    
radial_data_all3 = np.loadtxt(disp_mean_outfname3)
radial_data3 = radial_data_all3[bin_slice_val:,1]
radial_data3_std = radial_data_all3[bin_slice_val:,2]
popt3, pcov3 = curve_fit(func, xdata, radial_data3, maxfev=10000)
# label = '$%s x^{-%s}$' % tuple(f2str.format(t) for t in [popt[0], popt[1]])
u_r_outer_fit3 = func(fit_x_outer, popt3[0])
# normalizing both data by their respective first value
radial_data1_norm = radial_data1/radial_data1[0]
radial_data2_norm = radial_data2/radial_data2[0]


from math import log10, floor
def find_exp_base(number):
    exp = floor(log10(abs(number)))
    return round(number/10**exp, 2), exp 

base_num, exp = find_exp_base(popt1[0])
label1_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$'
base_num, exp = find_exp_base(popt2[0])
label2_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$'
base_num, exp = find_exp_base(popt3[0])
label3_txt = '$\\Sigma_{1}/\\mu_{m} =' + str(base_num) + ' \\times 10^{'+ str(exp) + '}$'


# plotting them together
fig = plt.figure(figsize=(8, 6), dpi=300)
# plt.scatter(bins_outer, binned_mean_disp_outer)
# plt.scatter(bins_outer, radial_data1_norm)
# plt.scatter(bins_outer, radial_data2_norm)
# plt.scatter(bins_outer, radial_data1)
# plt.scatter(bins_outer, radial_data2)
# plt.scatter(bins_outer, radial_data3)
plt.errorbar(bins_outer, radial_data1, yerr= radial_data1_std, ls = "None", marker = 'o', ms = 10, c = 'b', label = '$N_d = 5$')
plt.errorbar(bins_outer, radial_data2, yerr= radial_data2_std, ls = "None", marker = 'd', ms = 10, c = 'orange', label = '$N_d = 10$')
plt.errorbar(bins_outer, radial_data3, yerr= radial_data3_std, ls = "None", marker = 's', ms = 10, c = 'g', label = '$N_d = 15$')
# plt.plot(fit_x_outer, u_r_outer_fit1, c = 'b', label = '$ \\Sigma_{1} / \\mu_m = %0.2e $' % popt1[0])
# plt.plot(fit_x_outer, u_r_outer_fit2, c = 'orange', label = '$ \\Sigma_{1} / \\mu_m = %0.2e $' % popt2[0])
# plt.plot(fit_x_outer, u_r_outer_fit3, c = 'g', label = '$ \\Sigma_{1} / \\mu_m = %0.2e $' % popt3[0])
plt.plot(fit_x_outer, u_r_outer_fit1, c = 'b', label = label1_txt)
plt.plot(fit_x_outer, u_r_outer_fit2, c = 'orange', label = label2_txt)
plt.plot(fit_x_outer, u_r_outer_fit3, c = 'g', label = label3_txt)
plt.xlabel('Radial Distance', fontsize = 20)
plt.ylabel('Radial Displacement', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)

plt.xlim([bins_outer[0]-0.5,xlim])
# plt.xlim([12+0.5,xlim])
plt.subplots_adjust(top=0.96, left = 0.17, bottom = 0.15, right = 0.98)
plt.legend(loc='best', fontsize = 15)
plt.savefig('/home/abhinav/david/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'+"png/mean_displacement_scatter_outer_region_mean_"+'srand_'+str(srand)+"_" +str(temp_num_center1)+"_"+str(temp_num_center1*6)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+str(num-1)+".png")  # just plotting the last step data

plt.close()
plt.clf()

# ###############################################################################    
'''