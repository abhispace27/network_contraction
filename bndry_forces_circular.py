#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:25:15 2024
This code calculates the boundary forces for the circular boundary simulations.
It also calculates dipole moments - local and far field.

** also calculates paiwise disntaces and writes to a file under /means/ folder
** also reads the pairwise distance files and makes a plot vs N

@author: abhinav
"""
import numpy as np
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.pyplot as plt
import math

def angle(x1,y1,x2,y2):
    dx = x2-x1
    dy = y2-y1
    return(math.atan2(dy,dx))

############################################################################################################################
# finding the bending forces on all nodes
############################################################################################################################
def calc_bending_forces():
    FILAx = np.zeros((num_pts, 3))  # Forces in x-direction
    FILAy = np.zeros((num_pts, 3))  # Forces in y-direction
    
    for i in range(num_pts):
        for j in range(3):
            
            if conn_node[i,j] != unconnected_const and conn_node[i,(j+3) % 6] != unconnected_const:            
                tempcon1 = conn_node[i,j]
                tempcon2 = conn_node[i,(j+3) % 6]
                                
                delxij = -(xpos[-1,i] - xpos[-1,tempcon1])
                delyij = -(ypos[-1,i] - ypos[-1,tempcon1])
                rij = np.sqrt(delxij**2 + delyij**2)
                
                delxik = -(xpos[-1,i] - xpos[-1,tempcon2])
                delyik = -(ypos[-1,i] - ypos[-1,tempcon2])
                rik = np.sqrt(delxik**2 + delyik**2)
                
                # if rij == 0 or rik == 0:
                #     continue  # Avoid division by zero
                
                cross_jik = -(delxij * delyik - delyij * delxik)
                fjik = np.abs(cross_jik) / (rij * rik)
                
                if np.abs(cross_jik) > 0:
                    DfdX = ((delyik - delyij) * cross_jik / (np.abs(cross_jik) * rij * rik) -
                            fjik * delxij / rij**2 - fjik * delxik / rik**2)
                    DfdY = ((delxij - delxik) * cross_jik / (np.abs(cross_jik) * rij * rik) -
                            fjik * delyij / rij**2 - fjik * delyik / rik**2)
                    
                    FILAx[i, j] += kappa * (fjik / np.sqrt(1 - fjik**2)) * DfdX
                    FILAy[i, j] += kappa * (fjik / np.sqrt(1 - fjik**2)) * DfdY
            #================================================================================
            
            if conn_node[i,j] != unconnected_const and conn_node[conn_node[i,j],j] != unconnected_const:            
            # if tempcon1 != unconnected_const and tempcon2 != unconnected_const:            
                tempcon1=conn_node[i,j]
                tempcon2=conn_node[conn_node[i,j],j]
                DfdX=0.0
                DfdY=0.0
                delxij=delyij=delxmj=delymj=0.0
                rij=rmj=0.0

                delxij = -(xpos[-1,i]-(xpos[-1, tempcon1 ]))	
                delyij = -(ypos[-1,i]-(ypos[-1, tempcon1 ]))
                rij = np.sqrt( pow(delxij,2)+pow(delyij,2) )
    
                delxmj = (xpos[-1,tempcon1]-( xpos[-1, tempcon2 ] ))
                delymj = (ypos[-1,tempcon1]-( ypos[-1, tempcon2 ] ))    
                rmj = np.sqrt( pow(delxmj,2)+pow(delymj,2) )
    
                cross_mji=delxij*delymj-delyij*delxmj
                fmji=np.sqrt(pow(cross_mji,2))/(rij*rmj)
                # THETAmji=asin(fmji);
    
                if(np.abs(cross_mji)>0.0):			
    			
                    DfdX=-delymj*cross_mji/(np.abs(cross_mji)*rij*rmj)-delxij*fmji/(pow(rij,2))
                    DfdY=delxmj*cross_mji/(np.abs(cross_mji)*rij*rmj)-delyij*fmji/(pow(rij,2))
    
                    # // trying the derivative of sin theta squared and not theta squared:
                    FILAx[i, j]=FILAx[i, j]+kappa*((fmji)/( np.sqrt(1-pow(fmji,2)) ))*DfdX
                    FILAy[i, j]=FILAy[i, j]+kappa*((fmji)/( np.sqrt(1-pow(fmji,2)) ))*DfdY
            #================================================================================

            if(conn_node[i,(j+3)%6] != unconnected_const and conn_node[conn_node[i,(j+3)%6],(j+3)%6] != unconnected_const):
                tempcon1=conn_node[i,(j+3)%6]
                tempcon2=conn_node[conn_node[i,(j+3)%6],(j+3)%6]

                DfdX=0.0
                DfdY=0.0
                delxik=delyik=delxlk=delylk=0.0
                rik=rlk=0.0

                delxik = -(xpos[-1,i]-(xpos[-1, tempcon1 ]))	
                delyik = -(ypos[-1,i]-(ypos[-1, tempcon1 ]))
                rik = np.sqrt( pow(delxik,2)+pow(delyik,2) )

                delxlk = (xpos[-1,tempcon1]-( xpos[-1, tempcon2 ] ))
                delylk = (ypos[-1,tempcon1]-( ypos[-1, tempcon2 ] ))
                rlk = np.sqrt( pow(delxlk,2)+pow(delylk,2) )

                cross_lki=-(delyik*delxlk-delxik*delylk)
                flki=np.sqrt(pow(cross_lki,2))/(rik*rlk)
# 				THETAlki=asin(flki);

                if(np.abs(cross_lki) > 0.0):
                    DfdX=-delylk*cross_lki/(np.abs(cross_lki)*rik*rlk)-delxik*flki/(pow(rik,2));
                    DfdY=delxlk*cross_lki/(np.abs(cross_lki)*rik*rlk)-delyik*flki/(pow(rik,2));

                    FILAx[i, j]=FILAx[i, j]+kappa*((flki)/( np.sqrt(1-pow(flki,2)) ))*DfdX;
                    FILAy[i, j]=FILAy[i, j]+kappa*((flki)/( np.sqrt(1-pow(flki,2)) ))*DfdY;

    
    return FILAx, FILAy                  # these are the bending forces along x and y

############################################################################################################################
# end of functions
############################################################################################################################

dont_throw_flag = 0                   # set to 1 if you do not want to throw the simulations where there maybe a completely disconnected node

# L = 16
L = 64
# L = 128

if L == 64:
    unconnected_const = 9999
if L == 128:
    unconnected_const = -1

num_pts = L*L

num_center = 1
# num_center = 2
# num_center = 3
# num_center = 4
# num_center = 5
# num_center = 10
# num_center = 15
# num_center = 20
# num_center = 25
# num_center = 30
# num_center = 35
# num_center = 40


# num_center = 2
# num_center = 3
# num_center = 4
# num_center = 5

#num_center = 7
#num_center = 12
#num_center = 15
# num_center = 20
#num_center = 25
# num_center = 30
#num_center = 35
#num_center = 40
#num_center = 50
#num_center = 60
#num_center = 70
#num_center = 80
#num_center = 100

#num_center = 15
#num_center = 20
#num_center = 30
#num_center = 40
#num_center = 45

num_dip = num_center*6
num = 5+1                                                                        # the total number of force steps
num = 10+1
#num = 15+1
#num = 10+1
#num = 50+1
#num = 100+1
#num = 1000+1

#loc = 135

pbond = 1
# pbond = 0.9
# pbond = 0.8
# pbond = 0.75
# pbond = 0.7

#pbond = 0.5
pbond = 0.55
# pbond = 0.6

#pbond = 0.61
#pbond = 0.63
#pbond = 0.65
#pbond = 0.7

mu = 1
mu_c = 1
# mu = 1e-3
# mu_c = 1e-3

tol = 2e-6
tol = 1.0e-6
tol = 1.0e-7
# tol = 1.0e-8

kappa = 1e-6
# kappa = 5e-6
kappa = 1e-5
#kappa = 2e-5
# kappa = 5e-5

# kappa = 1e-4

#kappa = 1e-2

rlen = 0.6    
rlen = 0.7
rlen = 0.8
rlen = 0.9
rlen_txt = "%.4f" % rlen

buckle_flag = 0

radial_flag = 0    # for inner and outer boundaries different; dipoles are hexagons

radial_inner_flag = 0   # for dipoles with just radial forces and inner boundary different from outer

inner_radial_arp_flag = 0       # for hex bonds

cluster_flag = 0

cluster_new_flag = 1

hex_flag = 0  # these are for cluster 100 simulations where the dipoles have all 12 bonds present - including the outer hex bonds
hex_rand_flag = 0
cluster_radial_flag = 1

if hex_flag == 1 or hex_rand_flag == 1 or cluster_radial_flag == 1:
    cluster_new_flag = 1

radial_only_flag = 0

force_flag = 0
if force_flag == 1:
    force_val = 0.01
    # force_val = 0.015
    rlen_txt = "%.9f" % force_val
    
if buckle_flag == 1:
    mu_c = 0.1
    mu_c = 0.5

srand_flag = 1
if srand_flag == 1:
    srand = 116
    # srand = 113
    # srand = 122

diff_flag = 1
if diff_flag == 1:
    diff_network1 = 111111
    # diff_network1 = 222222
    # diff_network1 = 333333
    # diff_network1 = 444444
    # diff_network1 = 555555
    # diff_network1 = 666666
    # diff_network1 = 777777
    # diff_network1 = 888888
    # diff_network1 = 999999
    # diff_network1 = 900000

    # diff_network2 = 101111         # 40
    # diff_network2 = 111111       #30, 40
    # diff_network2 = 121111
    # diff_network2 = 131111
    # diff_network2 = 141111
    # diff_network2 = 151111      # 5; 10; 15; 20; 30 ; 40
    # diff_network2 = 161111
    # diff_network2 = 171111
    # diff_network2 = 181111
    diff_network2 = 191111

    # diff_network2 = 202222
    # diff_network2 = 212222
    # diff_network2 = 222222
    # diff_network2 = 232222
    # diff_network2 = 242222      # 20
    # diff_network2 = 252222
    # diff_network2 = 262222
    # diff_network2 = 272222      # 30 
    # diff_network2 = 282222
    # diff_network2 = 292222

    # diff_network2 = 303333
    # diff_network2 = 313333
    # diff_network2 = 323333
    # diff_network2 = 333333
    # diff_network2 = 343333      # 30
    # diff_network2 = 353333      # 30; 40
    # diff_network2 = 363333
    # diff_network2 = 373333      # 30
    # diff_network2 = 383333      # 30
    # diff_network2 = 393333

    # diff_network2 = 404444
    # diff_network2 = 414444
    # diff_network2 = 424444
    # diff_network2 = 434444       # 30
    # diff_network2 = 444444       # 20 ; 40
    # diff_network2 = 454444       # 5; 10; 15; 20; 30; 40
    # diff_network2 = 464444       # 30 
    # diff_network2 = 474444       # 20 ; 30 ; 40 
    # diff_network2 = 484444
    # diff_network2 = 494444

    # diff_network2 = 505555
    # diff_network2 = 515555      # 20
    # diff_network2 = 525555      # 30 
    # diff_network2 = 535555      # 20; 30 ; 40
    # diff_network2 = 545555
    # diff_network2 = 555555      # 30; 40
    # diff_network2 = 565555
    # diff_network2 = 575555      # 20 ; 40 
    # diff_network2 = 585555
    # diff_network2 = 595555

    # diff_network2 = 606666      # 20
    # diff_network2 = 616666
    # diff_network2 = 626666      # 40 
    # diff_network2 = 636666
    # diff_network2 = 646666
    # diff_network2 = 656666
    # diff_network2 = 666666
    # diff_network2 = 676666      # 40 
    # diff_network2 = 686666
    # diff_network2 = 696666      # 30

    # diff_network2 = 707777
    # diff_network2 = 717777
    # diff_network2 = 727777
    # diff_network2 = 737777
    # diff_network2 = 747777      # 40 
    # diff_network2 = 757777
    # diff_network2 = 767777
    # diff_network2 = 777777
    # diff_network2 = 787777      # 40
    # diff_network2 = 797777
    
    # diff_network2 = 808888      # 5; 10; 15; 20; 30 ; 40 
    # diff_network2 = 818888
    # diff_network2 = 828888
    # diff_network2 = 838888      # 40 
    # diff_network2 = 848888
    # diff_network2 = 858888
    # diff_network2 = 868888
    # diff_network2 = 878888      # 10 ; 15
    # diff_network2 = 888888
    # diff_network2 = 898888

    # diff_network2 = 909999
    # diff_network2 = 919999     # 40
    # diff_network2 = 929999
    # diff_network2 = 939999     # 30 
    # diff_network2 = 949999
    # diff_network2 = 959999
    # diff_network2 = 969999     # 30 ; 40
    # diff_network2 = 979999
    # diff_network2 = 989999
    # diff_network2 = 999999

    # diff_network2 = 900000
    # diff_network2 = 900001     # 30
    # diff_network2 = 900002
    # diff_network2 = 900003
    # diff_network2 = 900004     # 5; 10; 15; 20 ; 30 ; 40 
    # diff_network2 = 900005
    # diff_network2 = 900006
    # diff_network2 = 900007     #  20 
    # diff_network2 = 900008
    # diff_network2 = 900009

#diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111',
#                     '222222,202222', '222222,212222', '222222,222222', '222222,232222', '222222,242222', '222222,252222', '222222,262222', '222222,272222', '222222,282222', '222222,292222',
#                     '333333,303333', '333333,313333', '333333,323333', '333333,333333', '333333,343333', '333333,353333', '333333,363333', '333333,373333', '333333,383333', '333333,393333',
#                    '555555,505555', '555555,515555', '555555,525555', '555555,535555', '555555,545555', '555555,555555', '555555,565555', '555555,575555', '555555,585555', '555555,595555',
#                     '666666,606666', '666666,616666', '666666,626666', '666666,636666', '666666,646666', '666666,656666', '666666,666666', '666666,676666', '666666,686666', '666666,696666',
#                     '777777,707777', '777777,717777', '777777,727777', '777777,737777', '777777,747777', '777777,757777', '777777,767777', '777777,777777', '777777,787777', '777777,797777',
#                     '888888,808888', '888888,828888', '888888,828888', '888888,838888', '888888,848888', '888888,858888', '888888,868888', '888888,878888', '888888,888888', '888888,898888',
#                     '999999,909999', '999999,929999', '999999,929999', '999999,939999', '999999,949999', '999999,959999', '999999,969999', '999999,979999', '999999,989999', '999999,999999',
#                     '900000,900000', '900000,900001', '900000,900002', '900000,900003', '900000,900004', '900000,900005', '900000,900006', '900000,900007', '900000,900008', '900000,900009']


    diff_network = str(diff_network1)+","+str(diff_network2)

auto_flag = 1
if auto_flag == 1:
    import sys
    if pbond != 1:
        diff_network = sys.argv[1]
        print("Diff network: ", diff_network)
        num_center = int(sys.argv[2])    
        num_dip = num_center*6    
        if srand_flag == 1:
            srand = int(sys.argv[3])
    else:
        num_center = int(sys.argv[1])    
        num_dip = num_center*6    
        if srand_flag == 1:
            srand = int(sys.argv[2])

if srand != 112 and num_center >= 30:
    dont_throw_flag = 1

pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
kappa_str = "%.2e" % kappa # to write the filename
if kappa == 0:
    kappa_str = "0.00e+00" # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"

folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'
if pbond == 1 and inner_radial_arp_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'      # for radial and HEX bonds
if pbond == 1 and radial_only_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/'
if pbond == 1 and force_flag == 1:
    folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_force/'

cent_dist_folder = folder

if (pbond != 1):
    ranseed = '667720,601210'    # -- network 1
#    ranseed = '666666,777777'
#    ranseed = '152240,193915'   # -- not being used for second project on contractions !!
#    ranseed = '167720,101210'    # -- network 2
#    ranseed = '712240,203915'   # -- network 3
#    ranseed = '537206,701210'   # -- network 4

#    ranseed = '399020,800099'   # -- network 4
#    ranseed = '507225,102270'   # -- network 4
#    ranseed = '367725,927820'

    if kappa == 2e-7:
        kappa_fname = 'kappa2_2e-7/'
    elif kappa == 5e-7:
        kappa_fname = 'kappa2_5e-7/'
    elif kappa == 1e-6:
        kappa_fname = 'kappa2_e-6/'
    elif kappa == 1e-5:
        kappa_fname = 'kappa2_e-5/'
    elif kappa == 2e-6:
        kappa_fname = 'kappa2_2e-6/'
    elif kappa == 5e-6:
        kappa_fname = 'kappa2_5e-6/'
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
    
    # if diff_flag == 0 and srand_flag == 1:
        # folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if diff_flag == 1 and srand_flag == 1 and radial_inner_flag == 0 and inner_radial_arp_flag == 0 and cluster_new_flag == 0 and force_flag == 0 and radial_only_flag == 0:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif radial_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_radial/'+ kappa_fname +ranseed+"/"
        cent_dist_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif radial_inner_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif inner_radial_arp_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"        # for hex bonds
        cent_dist_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
        print('folder set to: ', folder)
    elif radial_only_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif cluster_flag == 1:
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 0:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif cluster_new_flag == 1 and hex_flag == 1 and hex_rand_flag == 0 and cluster_radial_flag == 0:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 1 and cluster_radial_flag == 0:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 1:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    elif force_flag == 1:
        folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_force/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
        cent_dist_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_force/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
        
    else:        
        folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'+ kappa_fname +ranseed+"/"
        cent_dist_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'+ kappa_fname +ranseed+"/"+str(srand)+"/"

#hist_base = base + folder + "png/strain_hist/"

#------------------------------------------------------------------------------
# reading the forces on all nodes ---------------------------------------------
#------------------------------------------------------------------------------
force_all = []
for i in range(0,num-1):
    force_fname = base + folder + "txt/force/Force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i+1)+"_"+str(num-1)+"_force.txt"
#    if diff_flag == 1:
#        force_fname = base + folder + "txt/force/Force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i+1)+"_"+str(num-1)+"_force.txt"
    # if srand_flag == 1 and num_center > 1 and pbond == 1:
    if srand_flag == 1 and pbond == 1 and L == 64:
        force_fname = base + folder + "txt/force/Force_srand_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i+1)+"_"+str(num-1)+"_force.txt"
    elif srand_flag == 1 and L != 64 and pbond == 1:
        force_fname = base + folder + "txt/force/Force_srand_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i+1)+"_"+str(num-1)+"_force.txt"
    elif srand_flag == 1 and L != 64 and pbond != 1:
        force_fname = base + folder + "txt/force/Force_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i+1)+"_"+str(num-1)+"_force.txt"
    elif srand_flag == 1 and force_flag == 1:
        force_fname = base + folder + "txt/force/Force_srand_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i+1)+"_"+str(num-1)+"_force.txt"
        
    print('Force file:  ', force_fname)
    force = np.loadtxt(force_fname)
    force_all.append(force)

force_all = np.array(force_all)
force_all = -1.0*force_all   # mulitplying by -1 because the file has derivative of energy and not forces right now
force_mag = np.sqrt(np.square(force_all[:,:,0])+np.square(force_all[:,:,1]))

#------------------------------------------------------------------------------
# reading the boundary nodes --------------------------------------------------
#------------------------------------------------------------------------------
fname_boundary = base+folder+'txt/area/Boundary_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# if srand_flag == 1 and num_center > 1 and pbond == 1:
# if srand_flag == 1 and pbond == 1 and L == 64:
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    fname_boundary = base+folder+'txt/area/Boundary_nodes_srand_'+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
elif srand_flag == 1 and pbond == 1 and  L != 64:
    fname_boundary = base+folder+'txt/area/Boundary_nodes_srand_'+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
elif srand_flag == 1 and pbond != 1 and  L != 64:
    fname_boundary = base+folder+'txt/area/Boundary_nodes_'+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
# elif force_flag == 1:
    
print("circular boundary filename: ",fname_boundary)
boundary_data = np.loadtxt(fname_boundary)
boundary_nodes = boundary_data[:,1]
boundary_nodes = boundary_nodes.astype(int)

#------------------------------------------------------------------------------
# finding forces on the boundary nodes
#------------------------------------------------------------------------------

force_boundary_x = force_all[:,boundary_nodes,0]
force_boundary_y = force_all[:,boundary_nodes,1]

force_boundary_vec = []
for i in range(0,np.shape(force_boundary_x)[0]):
    force_boundary_vec.append([])
    for j in range(0,np.shape(force_boundary_x)[1]):
        force_boundary_vec[i].append([force_boundary_x[i,j], force_boundary_y[i,j]])

force_boundary_vec = np.array(force_boundary_vec)
#------------------------------------------------------------------------------

step_val = (1 - rlen)/(num-1)
rlen_steps = [float(1-i*step_val) for i in range(0,num)]
rlen_steps = rlen_steps[1:]

#force_outfname = base+folder+"txt/area/"+"bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+force_str+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_e_dip0_1000.txt"
# force_outfname = base+folder+"txt/area/"+"bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if diff_flag == 1:
#     force_outfname = base+folder+"txt/area/"+"bndry_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    
# print('here sum of all forces on BOUNDARY nodes are written to file : ',force_outfname)    
# heading = '# rlen    Magnitude force'
# fmt = '%12.5f', '%15.7e' 

# #np.savetxt(force_outfname, np.column_stack((rlen_steps, total_bndry_force)), header = heading, fmt = fmt)
# np.savetxt(force_outfname, np.column_stack((rlen_steps, force_boundary_x, force_boundary_y)), header = heading, fmt = fmt)


force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    force_outfname = base+folder+"txt/area/"+"bndry_node_force_srand_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    force_outfname = base+folder+"txt/area/"+"bndry_node_force_srand_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    force_outfname = base+folder+"txt/area/"+"bndry_node_force_srand_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

#if diff_flag == 1:
#    force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('here all individual forces (includes bending force component) on BOUNDARY nodes are written to file : ',force_outfname)    
heading = '# nodes       Fx           Fy'
fmt = '%12d', '%15.7e', '%15.7e' 

np.savetxt(force_outfname, np.column_stack((boundary_nodes, force_boundary_x[-1], force_boundary_y[-1])), header = heading, fmt = fmt)

##############################################################################################################################################################
# Finding the unbalanced force on the diole nodes except on the center node
##############################################################################################################################################################

##############################################################################################################################################################
# Plotting the boundary force histogram as a ratio of total boundary force
##############################################################################################################################################################
'''
bndry_force_ratio = boundary_force_mag[-1]/total_bndry_force[-1]

outlier_lim_abs = 5e-6
outlier_exlcluded_bf = np.array(boundary_force_mag[-1])
outlier_exlcluded_bf = outlier_exlcluded_bf[abs(outlier_exlcluded_bf) <= outlier_lim_abs] 
print('Outlier excluded boundary force: ',np.sum(outlier_exlcluded_bf))

fig = plt.figure(figsize=(8, 4), dpi=300)
#mean_pos_abhi_str = 'Mean = ' + "%.2e" % np.mean(mean_pos_abhi)
#max_lim = 0.02      # works for p = 0.55
max_lim = np.max(bndry_force_ratio)
min_lim = np.min(bndry_force_ratio)
#binwidth = (max_lim-0)/100.
binwidth = 0.005
binwidth_str = "%.2e" % binwidth
bins = np.arange(min_lim-binwidth,max_lim+binwidth,binwidth)

n,bins,patches = plt.hist(bndry_force_ratio, bins = bins, histtype = 'step', alpha = 1.0)#, color='blue')#, label = "CG " + mean_pos_abhi_str)
plt.yscale('log')
#plt.title('Positive Strains')
#plt.legend()
plt.xlabel('Normalilzed Boundary Forces', fontsize = 15)
plt.ylabel('Number of Nodes', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.subplots_adjust(top=0.92, left = 0.12, bottom = 0.15, right = 0.96)

plt.savefig(base+folder+"png/strain_hist/b_force_hist_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+binwidth_str+"_"+str(num-1)+".png")  # just plotting the last step data
plt.close()
plt.clf()

# scatter plot of boundary force vs node number
fig = plt.figure(figsize=(8, 4), dpi=300)
plt.scatter(np.arange(0,len(boundary_force_mag[-1])),boundary_force_mag[-1])
plt.xlabel('Node Number', fontsize = 15)
plt.ylabel('Boundary Forces', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.ylim([-5e-5,5e-5])
#plt.ylim([-5e-6,5e-6])
plt.subplots_adjust(top=0.92, left = 0.12, bottom = 0.15, right = 0.96)

plt.savefig(base+folder+"png/strain_hist/b_force_scatter_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+mu_str+"_"+mu_c_str+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+binwidth_str+"_"+str(num-1)+".png")  # just plotting the last step data
plt.close()
plt.clf()
'''

##############################################################################################################################################################
# calculating forces from strains :: so this accounds for stretching and compression forces only. 
# IMPORTANT::                                It does NOT account for bending forces ...!!!!!
##############################################################################################################################################################
# reading the connection array
fname_connect = base+folder+'txt/strain/Lattice_connect_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    fname_connect = base+folder+'txt/strain/Lattice_connect_srand_'+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
elif srand_flag == 1 and pbond == 1 and L != 64:
    fname_connect = base+folder+'txt/strain/Lattice_connect_srand_'+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
elif srand_flag == 1 and pbond != 1 and L != 64:
    fname_connect = base+folder+'txt/strain/Lattice_connect_'+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name

#if diff_flag == 1:
#    fname_connect = base+folder+'txt/strain/Lattice_connect_'+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
    
print("connection filename: ",fname_connect)
conn_data = np.loadtxt(fname_connect)
conn_node = conn_data[:,3:9]
conn_node = (np.rint(conn_node)).astype(int)

# reading the p = 1 connection array (so that we can find the connectivity of all nodes since we do not want to throw any data right now: for non-connected dipole nodes)
p1_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/'      # for radial and HEX bonds
if L == 64:
    fname_connect_p1 = base+p1_folder+'txt/strain/Lattice_connect_srand_'+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_1.00_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
if L != 64:
    fname_connect_p1 = base+p1_folder+'txt/strain/Lattice_connect_srand_'+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_1.00_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name

print("connection filename: ",fname_connect_p1)
conn_data_p1 = np.loadtxt(fname_connect_p1)
conn_node_p1 = conn_data_p1[:,3:9]
conn_node_p1 = (np.rint(conn_node_p1)).astype(int)


pos_folder = base + folder + "txt/strain/Lattice_" 
strain_pos_folder = base + folder + "txt/strain/Strain_Lattice_" 

all_data = []
strain_all_data = []
rlen_arr = []

for i in range(0,num):     # for all positions
    fname = pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"        # input file name
    strain_fname = strain_pos_folder + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"  # strain file
    # reading the rest-length array
    fname_rlen = base+folder+'txt/rlen/rlen_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"        # input file name


    if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
        fname = pos_folder + 'srand_' + str(srand) + "_" + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"        # input file name
        strain_fname = strain_pos_folder + 'srand_' + str(srand) + "_" + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"  # strain file
        # reading the rest-length array
        fname_rlen = base+folder+'txt/rlen/rlen_srand_' + str(srand) + "_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"        # input file name
    elif srand_flag == 1 and pbond == 1 and L != 64:
        fname = pos_folder +'srand_' + str(srand) + "_"+ str(L) + "_" + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"        # input file name
        strain_fname = strain_pos_folder + 'srand_' + str(srand) + "_" + str(L) + "_" + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"  # strain file
        # reading the rest-length array
        fname_rlen = base+folder+'txt/rlen/rlen_srand_' + str(srand) + "_" + str(L) + "_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"        # input file name
    elif srand_flag == 1 and pbond != 1 and L != 64:
        fname = pos_folder + str(L) + "_" + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"        # input file name
        strain_fname = strain_pos_folder + str(L) + "_" + str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+"_force.txt"  # strain file
        # reading the rest-length array
        fname_rlen = base+folder+'txt/rlen/rlen_' + str(L) + "_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"        # input file name

    print("rest-length filename: ",fname_rlen)
    print("Node position file name: ",fname)
    print("strain file name: ", strain_fname)

    data = np.loadtxt(fname)
#    xpos = data[:,1]
#    ypos = data[:,2]
    all_data.append(data)
    
    strain_data = np.loadtxt(strain_fname)
    strain_all_data.append(strain_data)

    rlen_temp = np.loadtxt(fname_rlen)
    rlen_arr.append(rlen_temp)

rlen_arr = np.array(rlen_arr)

all_data = np.array(all_data)
xpos = all_data[:,:,1]
ypos = all_data[:,:,2]

strain_all_data = np.array(strain_all_data)
bond_len = strain_all_data[:,:,1:7]
#bond_len = np.reshape(bond_len,(num,6*L*L))
#strain = strain - 1.0    # subtracting restlength from bond length to get strain
strain = []
for i in range(0,len(bond_len)):
#    del strain[i][rm_loc]
    strain.append([])
    for j in range(0,np.shape(bond_len)[1]):    # iterating over all nodes
        strain[i].append([])
        for k in range(0,np.shape(bond_len)[2]):     # iterating over all bonds
            if bond_len[i,j,k] != 0:
#                strain[i][j].append(bond_len[i,j,k] - 1.0)
                strain[i][j].append(bond_len[i,j,k] - rlen_arr[i,j,k])
            else:
                strain[i][j].append(np.nan)
    
strain = np.array(strain)

y_force_node = []
x_force_node = []
bndry_force_arr = []
bndry_node = []
for k in range(0,len(all_data)):
    y_force_node.append([])
    x_force_node.append([])
    bndry_force_arr.append([])
    
    for i in range(0,num_pts):
        row = int(i/L)
        col = i%L
    #for i in range(0,60):
#        row = int(boundary_nodes_sorted[i]/L)
#        column = boundary_nodes_sorted[i]%L
#        node_num = boundary_nodes_sorted[i]
        
        # bottom boundary without the edge nodes: only 0, 3, 4th and 5th bonds can play a role in perpendiv=cular forces
        bonds_list = [0,1,2,3,4,5]
    
        force_node_y = 0
        force_node_x = 0
        for bond in bonds_list:
            strain_val = strain[k][i][bond]    # fourth bond strain
            conn_val = conn_node[i][bond]     # node that is the fourth bond connnected to

            if np.isnan(strain_val) == False:
                x1 = all_data[k][i,1]
                y1 = all_data[k][i,2]
                x2 = all_data[k][conn_val,1]
                y2 = all_data[k][conn_val,2]
                angle_val = angle(x1,y1,x2,y2)
                if strain_val >= 0:
#                    force_bond = -1*strain_val*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
                    force_bond = strain_val*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
                elif strain_val < 0:
#                    force_bond = -1*strain_val*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary
                    force_bond = strain_val*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary
                force_bond_y = force_bond*math.sin(angle_val)
                force_bond_x = force_bond*math.cos(angle_val)
                
                force_node_y = force_node_y + force_bond_y  # adding bond forces to force on nodes
                force_node_x = force_node_x + force_bond_x  # adding bond forces to force on nodes
                                
        y_force_node[k].append(force_node_y)
        x_force_node[k].append(force_node_x)
        
        # bottom boundary without corner nodes
        if (row == 0) and (i != 0) and (i != L-1): 
            bndry_force_arr[k].append(force_node_y)
            if k == 0:
                bndry_node.append(i)
        # right boundary  without corner nodes
        if (col == L-1) and (i != L-1) and (i != L*L-1): 
            bndry_force_arr[k].append(-1.*force_node_x)
            if k == 0:
                bndry_node.append(i)
        # left boundary  without corner nodes
        if (col == 0) and (i != 0) and (i != L*(L-1)): 
            bndry_force_arr[k].append(force_node_x)
            if k == 0:
                bndry_node.append(i)
        # top boundary  without corner nodes
        if (row == L-1) and (i != L*(L-1)) and (i != L*L-1): 
            bndry_force_arr[k].append(-1.*force_node_x)
            if k == 0:
                bndry_node.append(i)
        
# net force on all nodes using x and y forces:
y_force_node = np.array(y_force_node)
x_force_node = np.array(x_force_node)


# writing all forces to output file
for i in range(0,len(x_force_node)):
    force_outfname = base+folder+"txt/force/"+"all_nodes_forces_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"
    if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
        force_outfname = base+folder+"txt/force/"+"all_nodes_forces_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"
    elif srand_flag == 1 and pbond == 1 and L != 64:
        force_outfname = base+folder+"txt/force/"+"all_nodes_forces_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"
    elif srand_flag == 1 and pbond != 1 and L != 64:
        force_outfname = base+folder+"txt/force/"+"all_nodes_forces_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"

#    if diff_flag == 1:
#        force_outfname = base+folder+"txt/force/"+"all_nodes_forces_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"
        
    print('here all individual forces on ALL nodes are written to file : ',force_outfname)    
    heading = '# node      Fx                   Fy'
    fmt = '%5d', '%20.7e', '%20.7e' 
    nodes_list = np.arange(0,L*L)
    np.savetxt(force_outfname, np.column_stack((nodes_list, x_force_node[i], y_force_node[i])), header = heading, fmt = fmt)

# writing boundary forces to output file for each node
for i in range(0,len(x_force_node)):
    force_outfname = base+folder+"txt/force/"+"bndry_nodes_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"
    if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
        force_outfname = base+folder+"txt/force/"+"bndry_nodes_force_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"
    elif srand_flag == 1 and pbond == 1 and L != 64:
        force_outfname = base+folder+"txt/force/"+"bndry_nodes_force_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"
    elif srand_flag == 1 and pbond != 1 and L != 64:
        force_outfname = base+folder+"txt/force/"+"bndry_nodes_force_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"

#    if diff_flag == 1:
#        force_outfname = base+folder+"txt/force/"+"bndry_nodes_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(i)+"_"+str(num-1)+".txt"

    print('here all individual forces on BOUNDARY nodes are written to file : ',force_outfname)    
    heading = '# node      Force'
    fmt = '%5d', '%20.7e' 
    nodes_list = np.arange(0,L*L)
    np.savetxt(force_outfname, np.column_stack((bndry_node, bndry_force_arr[i])), header = heading, fmt = fmt)

# writing boundary forces to output file for each rlen step
sum_bndry_forces = np.sum(bndry_force_arr, axis = 1)
sum_bndry_forces = sum_bndry_forces[1:]
force_outfname = base+folder+"txt/force/"+"bndry_force_sum_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    force_outfname = base+folder+"txt/force/"+"bndry_force_sum_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    force_outfname = base+folder+"txt/force/"+"bndry_force_sum_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    force_outfname = base+folder+"txt/force/"+"bndry_force_sum_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    
#if diff_flag == 1:
#    force_outfname = base+folder+"txt/force/"+"bndry_force_sum_"+str()+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('here all boundary forces summed together are written to file : ',force_outfname)    
heading = '# rest-length      Force'
fmt = '%0.5f', '%20.7e' 
np.savetxt(force_outfname, np.column_stack((rlen_steps, sum_bndry_forces)), header = heading, fmt = fmt)

##############################################################################################################################################################
# calculating forces from strains :: so this accounds for stretching and compression forces only. 
# IMPORTANT::                                It does NOT account for bending forces ...!!!!!
# -----------------------------------------------------------ENDS HERE---------------------------------------------------------------------------------
##############################################################################################################################################################

#------------------------------------------------------------------------------
# converting the x and y forces on boundary nodes to forces that are perpendicular to boundary or radial from the center
#------------------------------------------------------------------------------

# finding coordinates of center of the boundary
center = 2080
if L == 128:
    center = 8383
#center = 2016
center_x = xpos[:,center]
center_y = ypos[:,center]
center_vec = [[center_x[i],center_y[i]] for i in range(0,len(center_x))]

# coordinates of the boundary nodes
boundary_x = xpos[:,boundary_nodes]
boundary_y = ypos[:,boundary_nodes]

# radial vector from boundary to center
radial_vec = []
for i in range(0,np.shape(boundary_x)[0]):
    radial_vec.append([])
    for j in range(0,np.shape(boundary_x)[1]):
        radial_x = center_x[i] - boundary_x[i,j]
        radial_y = center_y[i] - boundary_y[i,j]
        
        radial_vec[i].append([radial_x,radial_y])

radial_vec = np.array(radial_vec)

# finding the projection of the boundary forces on the radial vector
temp_boundary_force = []
for i in range(0,np.shape(force_boundary_vec)[0]):
    perp_force = np.dot(force_boundary_vec[i], np.transpose(radial_vec[i+1])) / np.linalg.norm(radial_vec[i+1],axis=1)
    temp_boundary_force.append(perp_force)   

# temp_boundary_force has a size of (num-1,num boundary nodes, num boundary nodes) because the previous operation in not element wise
# so we will extract the diagonal terms from it that is the right calculation
perp_boundary_force = []    
for i in range(0,np.shape(temp_boundary_force)[0]):
    temp_force = temp_boundary_force[i]
    perp_boundary_force.append(temp_force.diagonal())
    
    # perp_boundary_force.append([])
    # for j in range(0,np.shape(temp_boundary_force)[1]):
    #     for k in range(0,np.shape(temp_boundary_force)[2]):
    #         if j == k:
    #             perp_boundary_force[i].append(temp_boundary_force[i,j,k])
    
    
perp_boundary_force = np.array(perp_boundary_force)

#------------------------------------------------------------------------------
# writing results to file :: these force include the boundary force !! WTF?
#------------------------------------------------------------------------------

# writing force on each boundary node
force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    force_outfname = base+folder+"txt/area/"+"bndry_node_force_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    force_outfname = base+folder+"txt/area/"+"bndry_node_force_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    force_outfname = base+folder+"txt/area/"+"bndry_node_force_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

#if diff_flag == 1:
#    force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('here all individual forces on BOUNDARY nodes are written to file : ',force_outfname)    
heading = '# nodes       Fx              Fy              Perperndicular F'
fmt = '%12d', '%15.7e', '%15.7e', '%15.7e' 

np.savetxt(force_outfname, np.column_stack((boundary_nodes, force_boundary_x[-1], force_boundary_y[-1], perp_boundary_force[-1])), header = heading, fmt = fmt)


# writing total perpendicular force on the boundary
force_outfname = base+folder+"txt/area/"+"bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    force_outfname = base+folder+"txt/area/"+"bndry_force_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    force_outfname = base+folder+"txt/area/"+"bndry_force_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    force_outfname = base+folder+"txt/area/"+"bndry_force_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    
#if diff_flag == 1:
#    force_outfname = base+folder+"txt/area/"+"bndry_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('here all individual forces on BOUNDARY nodes are written to file : ',force_outfname)    
heading = '# Rest length       Total Perperndicular F'
fmt = '%15d', '%15.7e'

np.savetxt(force_outfname, np.column_stack((rlen_steps, np.sum(perp_boundary_force, axis=1))), header = heading, fmt = fmt)


#******************************************************************************************************************
# calculating the FORCE DIPOLE (TRACE) MOMENT AT THE BOUNDARY
#******************************************************************************************************************

# arrays store values for each node
xx_moment_sum_array = force_boundary_vec[-1,:,0] * radial_vec[-1,:,0]    # fx dot rx
yy_moment_sum_array = force_boundary_vec[-1,:,1] * radial_vec[-1,:,1]    # fy dot ry
xy_moment_sum_array = force_boundary_vec[-1,:,0] * radial_vec[-1,:,1]  
yx_moment_sum_array = force_boundary_vec[-1,:,1] * radial_vec[-1,:,0]

# single values after summing over all the boundary nodes
xx_moment_sum = np.sum(xx_moment_sum_array)
yy_moment_sum = np.sum(yy_moment_sum_array)
xy_moment_sum = np.sum(xy_moment_sum_array)
yx_moment_sum = np.sum(yx_moment_sum_array)

# trace of dipole tensor: one single value
trace_moment = xx_moment_sum + yy_moment_sum     # this should be the dipole moment at the boundary
trace_moment_array = xx_moment_sum_array + yy_moment_sum_array     # dipole moment at each node

#total sum of dipole tensor: one single value
total_moment = xx_moment_sum + yy_moment_sum + xy_moment_sum + yx_moment_sum



# writing the moment on each node to file
force_outfname = base+folder+"txt/area/"+"bndry_node_moment_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    force_outfname = base+folder+"txt/area/"+"bndry_node_moment_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    force_outfname = base+folder+"txt/area/"+"bndry_node_moment_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    force_outfname = base+folder+"txt/area/"+"bndry_node_moment_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    

print('here all individual moments on BOUNDARY nodes are written to file : ',force_outfname)    
heading = '# nodes       Fx              Fy              Xpos_vec          Ypos_vec              x_moment           y_moment             f dot r'
fmt = '%12d', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e'  

np.savetxt(force_outfname, np.column_stack((boundary_nodes, force_boundary_x[-1], force_boundary_y[-1], radial_vec[-1,:,0], radial_vec[-1,:,1], xx_moment_sum_array, yy_moment_sum_array, trace_moment_array)), header = heading, fmt = fmt)

# writing the total dipole moment at the boundary is written to file
fname_dip_moment = base+folder+"txt/area/"+"bndry_node_dipole_moment_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    fname_dip_moment = base+folder+"txt/area/"+"bndry_node_dipole_moment_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
elif srand_flag == 1 and pbond == 1 and L != 64:
    fname_dip_moment = base+folder+"txt/area/"+"bndry_node_dipole_moment_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
elif srand_flag == 1 and pbond != 1 and L != 64:
    fname_dip_moment = base+folder+"txt/area/"+"bndry_node_dipole_moment_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
    
print('here total BOUNDARY dipole moment is written to file : ',force_outfname)    
heading = 'rlen       f dot r'
fmt = '%12.3f', '%15.7e'  

np.savetxt(fname_dip_moment, np.column_stack((rlen, trace_moment)), header = heading, fmt = fmt)


############################################################################################################################
# finding the bending forces on all nodes
############################################################################################################################
# but first need to find bending forces on all nodes
# force_bend_x = force_all[-1,:,0] - x_force_node[-1,:]        # total force - stretching force = bending force (lol)
# force_bend_y = force_all[-1,:,1] - y_force_node[-1,:]
# # writing the bending forces to file:
# bend_force_outfname = base+folder+"txt/area/"+"bend_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
#     bend_force_outfname = base+folder+"txt/area/"+"bend_force_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# elif srand_flag == 1 and pbond == 1 and L != 64:
#     bend_force_outfname = base+folder+"txt/area/"+"bend_force_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# elif srand_flag == 1 and pbond != 1 and L != 64:
#     bend_force_outfname = base+folder+"txt/area/"+"bend_force_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

# print('bending forces on each node written to file : ',bend_force_outfname)    
# heading = '# node       F_bend_x         F_bend_y'
# fmt = '%12d', '%15.7e', '%15.7e' 
# np.savetxt(bend_force_outfname, np.column_stack(( np.arange(0,num_pts), force_bend_x, force_bend_y)), header = heading, fmt = fmt)


force_bend_x_func, force_bend_y_func = calc_bending_forces()
force_bend_x_func_node = np.sum(force_bend_x_func, axis=1) 
force_bend_y_func_node = np.sum(force_bend_y_func, axis=1) 

# writing the newly calculated bending forces to file:
bend_force_test_outfname = base+folder+"txt/area/"+"bend_force_test_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    bend_force_test_outfname = base+folder+"txt/area/"+"bend_force_test_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    bend_force_test_outfname = base+folder+"txt/area/"+"bend_force_test_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    bend_force_test_outfname = base+folder+"txt/area/"+"bend_force_test_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('bending forces on each node (test) written to file : ',bend_force_test_outfname)    
heading = '# node       F_bend_x         F_bend_y'
fmt = '%12d', '%15.7e', '%15.7e' 
np.savetxt(bend_force_test_outfname, np.column_stack(( np.arange(0,num_pts), force_bend_x_func_node, force_bend_y_func_node)), header = heading, fmt = fmt)

############################################################################################################################

#******************************************************************************************************************
# calculating the LOCAL FORCE DIPOLE MOMENT 
#******************************************************************************************************************
# reading the unique dipole node array
fname_dipole = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    fname_dipole = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+ 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
elif srand_flag == 1 and pbond == 1 and L != 64:
    fname_dipole = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+ 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
elif srand_flag == 1 and pbond != 1 and L != 64:
    fname_dipole = base+folder+'txt/dip_nodes/Unique_dipole_nodes_'+ str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"        # input file name
    

#if diff_flag == 1:
#    fname_connect = base+folder+'txt/strain/Lattice_connect_'+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"__"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(0)+"_"+str(num-1)+"_force.txt"        # input file name
    
print("dipole node filename: ",fname_dipole)
dipole_data = np.loadtxt(fname_dipole)
dipole_nodes = dipole_data[:,1]
dipole_nodes = dipole_nodes.astype(int)

num_center_list = []
for i in range(0,len(dipole_nodes)):
    if i%7 == 0:
        num_center_list.append(dipole_nodes[i])
num_center_list = np.array(num_center_list)

center_force_x = force_all[-1,num_center_list,0]
center_force_y = force_all[-1,num_center_list,1]

center_force_mag = np.sqrt(np.square(center_force_x) + np.square(center_force_y))  # this just checks the net force and should be as close to zero as possible

# taking a shortcut by just using the central node of the dipole to find unbalanced forces
force_unbalanced = []
force_unbalanced_list = []
d_loc_short = []
dist_dip_node_list = []
for i in range(0,len(num_center_list)):
    force_unbalanced.append([])
    force_unbalanced[i] = 0
    d_loc_short.append([])
    dist_dip_node_list.append([])
    for j in range(0,6):
        bond_force = mu*strain[-1,num_center_list[i],j]   # force due to bond strain
        force_unbalanced[i] = force_unbalanced[i] + bond_force     # addingn unbalnced forces on the central node
        force_unbalanced_list.append(bond_force)
        # distance to each node from the central node
        if conn_node[num_center_list[i],j] != unconnected_const:
            dist_x = xpos[-1,num_center_list[i]] - xpos[-1,conn_node[num_center_list[i],j]]
            dist_y = ypos[-1,num_center_list[i]] - ypos[-1,conn_node[num_center_list[i],j]]
        elif dont_throw_flag == 1:
            dist_x = xpos[-1,num_center_list[i]] - xpos[-1,conn_node_p1[num_center_list[i],j]]
            dist_y = ypos[-1,num_center_list[i]] - ypos[-1,conn_node_p1[num_center_list[i],j]]

        dist_dip_node = np.sqrt(np.square(dist_x) + np.square(dist_y))
        dist_dip_node_list[i].append(dist_dip_node)  # list of dipole node distance from central node; 2d array format [num_center,6]
        #using the distance to find d_loc
        d_loc_short[i].append(dist_dip_node*bond_force)

force_unbalanced = np.array(force_unbalanced)
d_loc_short_sum = np.sum(d_loc_short)


# validating that the above short-cut is correct by finding unbalanced forces on the outer dipole nodes
y_force_unbalanced = []
x_force_unbalanced = []
strain_loc = []
test_force_list = []
dloc_j_list = []

j = 0
for i in range(0,len(dipole_nodes)):

    if i%7 == 0: 
        last_center = dipole_nodes[i]
        index = int(i/7)
        y_force_unbalanced.append([])
        x_force_unbalanced.append([])
        dip_nodes_outer = dipole_nodes[i+1:i+7]
        test_force_list.append([])
        dloc_j_list.append([])
        
        dip_mom = 0
                
    if i%7 != 0:

        # find the orientation of the spring connected to dipole node
        x_last_center = xpos[-1,last_center]
        y_last_center = ypos[-1,last_center]
        x_node = xpos[-1,dipole_nodes[i]]
        y_node = ypos[-1,dipole_nodes[i]]
        
        spring_vec = [x_node-x_last_center, y_node-y_last_center]
        
        angle_spring = angle(x_node, y_node, x_last_center, y_last_center)

        if force_flag == 1:            
            angle_val = angle(x_last_center, y_last_center, x_node, y_node)
            force_val_y = (num-1)*force_val*np.sin(angle_val)
            force_val_x = (num-1)*force_val*np.cos(angle_val)
            dip_mom += np.dot([force_val_x, force_val_y], spring_vec) #/ np.linalg.norm(spring_vec)
        
        bonds_list = [0,1,2,3,4,5]
    
        force_node_y = 0
        force_node_x = 0
        test_force_node = 0
        
        for bond in bonds_list:
            strain_val = strain[-1][dipole_nodes[i]][bond]    # fourth bond strain
            conn_val = conn_node[dipole_nodes[i]][bond]     # node that is the fourth bond connnected to
            
            if conn_val == last_center:
                # strain_active_force = 1. - (rlen + strain_val)
                strain_active_force = (rlen + strain_val) - 1.
                # strain_active_force = strain_val
                # print(dipole_nodes[i], bond, strain_active_force, strain_val)
                # print('here:: ', dipole_nodes[i], bond, strain_active_force, rlen, strain_val)
            else:
                strain_active_force = strain_val

            
            if np.isnan(strain_val) == False:# and conn_val != last_center:
                                
                x1 = all_data[-1][dipole_nodes[i],1]
                y1 = all_data[-1][dipole_nodes[i],2]
                x2 = all_data[-1][conn_val,1]
                y2 = all_data[-1][conn_val,2]
                angle_val = angle(x1, y1, x2, y2)

                if strain_val >= 0:
#                    force_bond = -1*strain_val*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
#                    force_bond = strain_val*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
                    force_bond = strain_active_force*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
                elif strain_val < 0:
#                    force_bond = -1*strain_val*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary
#                    force_bond = strain_val*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary
                    force_bond = strain_active_force*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary

                force_bond_y = force_bond*np.sin(angle_val)
                force_bond_x = force_bond*np.cos(angle_val)
                
                # if dipole_nodes[i] == 2079:
                #     print(bond, force_bond_x, force_bond_y, np.sqrt(np.square(force_bond_x) + np.square(force_bond_y)), force_bond)
                
                force_node_y = force_node_y + force_bond_y  # adding bond forces to force on nodes
                force_node_x = force_node_x + force_bond_x  # adding bond forces to force on nodes
                
                force_unbalanced_temp = force_bond*np.cos(angle_spring)
                # finding the projection of spring force along the radial outward vector                
                test_force = np.dot([force_bond_x, force_bond_y], spring_vec) / np.linalg.norm(spring_vec)
                # test_force = np.dot([force_bond_x, force_bond_y], spring_vec)
                
                # if i%7 == 1 and bond == 3:
                #     test_force = -1*test_force
                # if i%7 == 2 and bond == 4:
                #     test_force = -1*test_force
                # if i%7 == 3 and bond == 5:
                #     test_force = -1*test_force
                # if i%7 == 4 and bond == 0:
                #     test_force = -1*test_force
                # if i%7 == 5 and bond == 1:
                #     test_force = -1*test_force
                # if i%7 == 6 and bond == 2:
                #     test_force = -1*test_force
                
                # if dipole_nodes[i] == 2403:
                #     print(bond, [force_bond_x, force_bond_y])#, spring_vec)
                #     print(bond, force_bond, angle_val, force_unbalanced_temp, test_force, np.dot([force_bond_x, force_bond_y], spring_vec), np.linalg.norm(spring_vec))
                #     print(test_force)
                #     print(' ')
                
                test_force_node = test_force_node + test_force
                    
        test_force_list[index].append(test_force_node)
                    
        # adding bending forces to the total force on the node:
        # print(dipole_nodes[i], ' \t ' ,force_node_x, ' \t ' ,force_node_y, ' \t ' ,force_bend_x_func_node[dipole_nodes[i]], ' \t ' ,force_bend_y_func_node[dipole_nodes[i]])
        force_node_x = force_node_x + force_bend_x_func_node[dipole_nodes[i]]
        force_node_y = force_node_y + force_bend_y_func_node[dipole_nodes[i]]
        # print(dipole_nodes[i], ' \t ' ,force_node_x, ' \t ' ,force_node_y)
        
        # doing the dot product to find dipole moment: f_j <dot> r_ij
        dloc_j = np.dot([force_node_x, force_node_y], spring_vec)
        dloc_j_list[index].append(dloc_j)


        # y_force_unbalanced[index].append(force_node_y)
        # x_force_unbalanced[index].append(force_node_x)
        
#        y_force_unbalanced[i].append(force_node_y)
#        x_force_unbalanced[i].append(force_node_x)

        j += 1

test_force_list = np.array(test_force_list)
dloc_total = np.sum(dloc_j_list)

# y_force_unbalanced = np.array(y_force_unbalanced)
# x_force_unbalanced = np.array(x_force_unbalanced)

# force_mag_unbalanced_list = np.sqrt(np.square(x_force_unbalanced) + np.square(y_force_unbalanced))
# force_mag_unbalanced = np.sum(force_mag_unbalanced_list)

# finding dipole moment using unblanced forces on the outer nodes
# d_loc_unbalanced_list = dist_dip_node_list*force_mag_unbalanced_list
# d_loc_unbalanced = np.sum(d_loc_unbalanced_list)

d_loc_test_array = dist_dip_node_list*test_force_list
d_loc_test = np.sum(d_loc_test_array)
if force_flag == 1:
    d_loc_test = dip_mom




print('D_loc is : ', d_loc_test)
print('new D_loc is : ', dloc_total)

d_far_d_loc_ratio = trace_moment/d_loc_test

# writing the total local dipole moment is written to file
fname_dip_moment = base+folder+"txt/area/"+"dipole_moment_ratio_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    fname_dip_moment = base+folder+"txt/area/"+"dipole_moment_ratio_"+ 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    fname_dip_moment = base+folder+"txt/area/"+"dipole_moment_ratio_"+ 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    fname_dip_moment = base+folder+"txt/area/"+"dipole_moment_ratio_"+ str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    

print('here total LOCAL dipole moment is written to file : ',force_outfname)    
heading = 'rlen           D_loc             D_far         D_far/D_loc'
fmt = '%12.3f', '%15.7e', '%15.7e'  , '%15.7e'  

np.savetxt(fname_dip_moment, np.column_stack((rlen, d_loc_test, trace_moment, d_far_d_loc_ratio)), header = heading, fmt = fmt)

trace_moment_bndry_arr = xx_moment_sum_array + yy_moment_sum_array
index_pos = np.where(trace_moment_bndry_arr > 0)
pos_val_trace = np.sum(trace_moment_bndry_arr[index_pos])
print('the positive dfars add up to: ', pos_val_trace)

# 
#******************************************************************************************************************
#------------------------------------------------------------------------------
# making a scatter plot of boundary forces
#------------------------------------------------------------------------------
fname_plot = base+folder+"png/strain_hist/"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    fname_plot = base+folder+"png/strain_hist/"+ str(srand) + "_" 
elif L != 64:
    fname_plot = base+folder+"png/strain_hist/"+ str(srand) + "_"+ str(L) + "_" 
    

# normalization value using the total boundary force
total_bf = np.sum(perp_boundary_force, axis=1)
total_bf_scaled = perp_boundary_force[-1]/total_bf[-1]

# '''
# just commented to speed up things

# # plotting all the boundary forces on each node as a scatter plot
plt.figure(figsize=(8, 5), dpi=300)

plt.scatter(np.arange(0,np.shape(perp_boundary_force)[1]),perp_boundary_force[-1], marker = '.', s = 500, c = 'g', alpha=0.5, label = pbond_string)# Network 1')    

#plt.yscale('log')
#plt.tick_params(axis='x', which='minor', bottom=False, top=False, labelbottom=False)

plt.xlabel("Node number", fontsize = 15)
plt.ylabel("Boundary Force", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

#plt.title('$\widetilde{\kappa} = 10^{-6}$', fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.15,top=0.96,bottom=0.12)

plt.legend(loc = 'lower right')
plt.savefig(fname_plot+"bndry_force_scatter_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()

# # plotting all the normalized boundary forces on each node as a scatter plot
plt.figure(figsize=(8, 5), dpi=300)

plt.scatter(np.arange(0,np.shape(perp_boundary_force)[1]),total_bf_scaled, marker = '.', s = 500, c = 'g', alpha=0.5, label = pbond_string)# Network 1')    

#plt.yscale('log')

#plt.tick_params(axis='x', which='minor', bottom=False, top=False, labelbottom=False)

plt.xlabel("Node number", fontsize = 15)
plt.ylabel("Boundary Force / Total Boundary Force", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

#plt.title('$\widetilde{\kappa} = 10^{-6}$', fontsize = 15)
plt.subplots_adjust(right=0.96,left=0.15,top=0.96,bottom=0.12)

plt.legend(loc = 'lower right')
plt.savefig(fname_plot+"bndry_force_scatter_normalized_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()

# # limiting the boundary forces y axis
lim = 0.1
lim_str = "%.2e" % lim

plt.figure(figsize=(8, 5), dpi=300)

plt.scatter(np.arange(0,np.shape(perp_boundary_force)[1]),total_bf_scaled, marker = '.', s = 500, c = 'g', alpha=0.5, label = pbond_string)# Network 1')    

#plt.yscale('log')

#plt.tick_params(axis='x', which='minor', bottom=False, top=False, labelbottom=False)

plt.xlabel("Node number", fontsize = 15)
plt.ylabel("Boundary Force / Total Boundary Force", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.ylim([-lim,lim])

plt.subplots_adjust(right=0.96,left=0.15,top=0.96,bottom=0.12)

plt.legend(loc = 'lower right')
plt.savefig(fname_plot+"bndry_force_scatter_lim_normalized_"+lim_str+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()


# # plotting positive boundary forces on each node as a scatter plot on a log scale
yval = perp_boundary_force[-1]
yval_pos = []
yval_neg = []
for i in range(0,len(yval)):
    if yval[i] > 0:
        yval_pos.append(yval[i])
    elif yval[i] < 0:
        yval_neg.append(yval[i])

yval_pos = np.array(yval_pos)
yval_neg = np.array(yval_neg)


# plotting positive forces
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(np.arange(0,len(yval_pos)),yval_pos, marker = '.', s = 500, c = 'g', alpha=0.5, label = pbond_string)# Network 1')    
plt.yscale('log')
if pbond == 0.5 and (num_center == 1 or num_center == 3 or num_center == 5 or num_center == 7  or num_center == 10):
    plt.ylim([1e-11,1e-6])
elif pbond == 0.6 and (num_center == 1 or num_center == 3 or num_center == 5):# or num_center == 7  or num_center == 10):
    plt.ylim([1e-10,1e-4])
elif pbond == 0.6 and (num_center == 7  or num_center == 10):
    plt.ylim([1e-8,1e-2])
elif pbond == 0.7 and (num_center == 1):# or num_center == 3 or num_center == 5 or num_center == 7  or num_center == 10):
    plt.ylim([1e-8,5e-3])
elif pbond == 0.7 and (num_center == 3 or num_center == 5 or num_center == 7 or num_center == 10):
    plt.ylim([1e-7,2e-2])
#plt.ylim([1e-11,1e-5])
plt.title('Positive forces', fontsize = 15)
plt.xlabel("Node number", fontsize = 15)
plt.ylabel("Boundary Force", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.legend(loc = 'lower right')
plt.savefig(fname_plot+"bndry_force_scatter_positives_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()

# plotting negative forces
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(np.arange(0,len(yval_neg)),-1*yval_neg, marker = '.', s = 500, c = 'g', alpha=0.5, label = pbond_string)# Network 1')    
plt.yscale('log')

plt.title('Negative forces', fontsize = 15)
plt.xlabel("Node number", fontsize = 15)
plt.ylabel("Boundary Force", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
if pbond == 0.5 and (num_center == 1 or num_center == 3 or num_center == 5  or num_center == 7  or num_center == 10):
    plt.ylim([1e-11,1e-6])
elif pbond == 0.6 and (num_center == 1 or num_center == 3 or num_center == 5):
    plt.ylim([1e-10,1e-4])
elif pbond == 0.6 and (num_center == 7  or num_center == 10):
    plt.ylim([1e-8,5e-3])
elif pbond == 0.7 and (num_center == 1):# or num_center == 3 or num_center == 5 or num_center == 7  or num_center == 10):
    plt.ylim([1e-8,5e-3])
elif pbond == 0.7 and (num_center == 3 or num_center == 5 or num_center == 7  or num_center == 10):
    plt.ylim([1e-7,2e-2])

plt.legend(loc = 'lower right')
plt.savefig(fname_plot+"bndry_force_scatter_negatives_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()

# # plotting positive normalized boundary forces (nomralized by the total boundary force) on each node as a scatter plot on a log scale
yval_norm = total_bf_scaled
yval_pos_norm = []
yval_neg_norm = []
for i in range(0,len(yval)):
    if yval_norm[i] > 0:
        yval_pos_norm.append(yval[i])
    elif yval_norm[i] < 0:
        yval_neg_norm.append(yval[i])

yval_pos_norm = np.array(yval_pos_norm)
yval_neg_norm = np.array(yval_neg_norm)

# plotting positive forces
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(np.arange(0,len(yval_pos_norm)),yval_pos_norm, marker = '.', s = 500, c = 'g', alpha=0.5, label = pbond_string)# Network 1')    
plt.yscale('log')

plt.title('Positive forces', fontsize = 15)
plt.xlabel("Node number", fontsize = 15)
plt.ylabel("Boundary Force / Total Boundary Force", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.legend(loc = 'lower right')
plt.savefig(fname_plot+"bndry_force_scatter_positives_norm_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()

# plotting negative forces
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(np.arange(0,len(yval_neg_norm)),-1*yval_neg_norm, marker = '.', s = 500, c = 'g', alpha=0.5, label = pbond_string)# Network 1')    
plt.yscale('log')

plt.title('Negative forces', fontsize = 15)
plt.xlabel("Node number", fontsize = 15)
plt.ylabel("Boundary Force / Total Boundary Force", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.legend(loc = 'lower right')
plt.savefig(fname_plot+"bndry_force_scatter_negatives_norm_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()
# '''

###############################################################################
if pbond != 1:
    plot_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'+ kappa_fname +ranseed+"/"
    if diff_flag == 0 and srand_flag == 1:
        plot_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if diff_flag == 1 and srand_flag == 1 and radial_inner_flag == 0 and inner_radial_arp_flag == 0:
        plot_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if radial_flag == 1:
        plot_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_radial/'+ kappa_fname +ranseed+"/"
    if radial_inner_flag == 1:
        plot_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if radial_only_flag == 1:
        plot_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if cluster_flag == 1:
        plot_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 0:
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if cluster_new_flag == 1 and hex_flag == 1 and hex_rand_flag == 0 and cluster_radial_flag == 0:
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 1 and cluster_radial_flag == 0:
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"
    if cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 1:
        plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"

if pbond == 1 and inner_radial_arp_flag == 1:
    plot_folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp/png/'


# for hexagonal dipoles in networks with differing inner and outer boundaries
# manually plotting the data for circle_radial dipoles for N = 1, 5, 10
# x = [1,5,10]
# y1 = [6.5756805e-07, 2.3263706e-05, 7.2418898e-05]
# y2 = [7.6859217e-06, 1.2690729e-05, 5.2063915e-05]
# y3 = [-3.0512987e-07, 5.0379647e-05, 8.2555027e-05]
# y4 = [6.2506471e-06, 1.1678054e-05, 2.4499141e-05]

# plt.figure(figsize=(8, 5), dpi=300)
# plt.scatter(x,y1, marker = '.', s = 500, c = 'g', alpha=0.5, label = 'Network 1')    
# plt.scatter(x,y2, marker = '*', s = 50, c = 'b', alpha=0.5, label = 'Network 2')    
# plt.scatter(x,y3, marker = '>', s = 50, c = 'r', alpha=0.5, label = 'Network 3')    
# plt.scatter(x,y4, marker = 's', s = 50, c = 'k', alpha=0.5, label = 'Network 4')    

# #plt.title('Negative forces', fontsize = 15)
# plt.xlabel("Number of Dipoles", fontsize = 15)
# plt.ylabel("Boundary Force", fontsize = 15)
# plt.xticks(fontsize = 15)
# plt.yticks(fontsize = 15)

# plt.legend(loc = 'lower right')
# plt.savefig(plot_folder+"test.png")

# plt.clf()
# plt.close()


# for hexagonal dipoles in networks with differing inner and outer boundaries
# manually plotting the data for circle_radial dipoles for N = 1, 5, 10

# x = [1, 5, 10]
# # x = [1, 5, 10, 15, 20, 30, 40]
# if srand == 112:
#     y1 = [8.5316422e-06, 6.6953873e-06, 2.0660988e-05]
#     y2 = [9.7737342e-06, 1.1752456e-05, 2.0395879e-05]
#     y3 = [2.4290134e-05, 4.4942009e-05, 4.1483143e-03]
#     y4 = [6.7480717e-06, 1.6513861e-05, 2.4733604e-05]
#     #y5 = [2.2840178e-06, 1.7134595e-05, 3.9190374e-05]
#     #y6 = [2.2979636e-06, 1.1854762e-04, ]
    
#     y1 = [1.2286321e-05, 8.5397530e-05, 3.2087048e-04, 2.6854986e-04, 3.7777079e-02, 4.8578413e-02, 1.4331763e-01] 
#     y2 = [6.2455561e-06, 7.9587991e-05, 1.2560788e-04, 3.8911970e-03, 1.0523099e-02, 3.6765319e-02, 1.1435311e-01] 
#     y3 = [5.2805927e-05, 2.2244517e-04, 8.9267557e-05, 4.1653884e-04, 1.5981415e-02, 3.1998007e-02, 8.3256502e-02] 
#     y4 = [1.5157107e-05, 1.4000984e-04, 1.7058211e-04, 1.4506835e-04, 5.8977750e-04, 3.4563500e-02, 6.8732812e-02] 
    
# elif srand == 113:
#     y1 = [8.9193601e-06, 1.3427016e-05, 1.9100549e-05]
#     y2 = [9.6116900e-06, 1.5142519e-05, 4.2086230e-05]
#     y3 = [2.1116260e-05, 3.7035323e-05, 7.7391401e-05]
#     y4 = [7.2870567e-06, 1.0353375e-05, 3.1439133e-05]
# elif srand == 114:
#     y1 = [2.1964590e-06, 9.1792276e-06, 1.5537983e-05]
#     y2 = [1.2068628e-06, 1.5199709e-05, 3.4089756e-05]
#     y3 = [5.0780186e-06, 3.1929042e-05, 3.6857524e-05]
#     y4 = [2.0740104e-06, 1.0741810e-06, 1.4292047e-05]

# plt.figure(figsize=(8, 5), dpi=300)
# plt.scatter(x,y1, marker = '.', s = 500, c = 'g', alpha=0.5, label = 'Network 1')    
# plt.scatter(x,y2, marker = '*', s = 50, c = 'b', alpha=0.5, label = 'Network 2')    
# if srand != 112 or pbond == 0.61:
#     plt.scatter(x,y3, marker = '>', s = 50, c = 'r', alpha=0.5, label = 'Network 3')    
# plt.scatter(x,y4, marker = 's', s = 50, c = 'k', alpha=0.5, label = 'Network 4')    
# #if srand == 112:
# #    plt.scatter(x,y5, marker = 'o', s = 50, c = 'y', alpha=0.5, label = 'Network 5')    

# #plt.title('Negative forces', fontsize = 15)
# plt.xlabel("Number of Dipoles", fontsize = 15)
# plt.ylabel("Boundary Force", fontsize = 15)
# plt.xticks(fontsize = 15)
# plt.yticks(fontsize = 15)

# plt.legend(loc = 'upper left')
# plt.savefig(plot_folder+"test_p"+pbond_string+"_"+str(srand)+".png")

# plt.clf()
# plt.close()

# # for hexagonal dipoles in networks with differing inner and outer boundaries
# # manually plotting the data for circle_radial dipoles for N = 1, 5, 10  on log scale
# plt.figure(figsize=(8, 5), dpi=300)
# plt.scatter(x,y1, marker = '.', s = 500, c = 'g', alpha=0.5, label = 'Network 1')    
# plt.scatter(x,y2, marker = '*', s = 50, c = 'b', alpha=0.5, label = 'Network 2')    
# if srand != 112 or pbond == 0.61:
#     plt.scatter(x,y3, marker = '>', s = 50, c = 'r', alpha=0.5, label = 'Network 3')    
# plt.scatter(x,y4, marker = 's', s = 50, c = 'k', alpha=0.5, label = 'Network 4')    
# #if srand == 12:
# #    plt.scatter(x,y5, marker = 'o', s = 50, c = 'y', alpha=0.5, label = 'Network 5')    

# #plt.title('Negative forces', fontsize = 15)
# plt.xlabel("Number of Dipoles", fontsize = 15)
# plt.ylabel("Boundary Force", fontsize = 15)
# plt.xticks(fontsize = 15)
# plt.yticks(fontsize = 15)
# plt.yscale('log')
# #plt.xscale('log')
# plt.legend(loc = 'upper left')
# plt.savefig(plot_folder+"log_test_p"+pbond_string+"_"+str(srand)+".png")

# plt.clf()
# plt.close()


# # for hexagonal dipoles in networks with differing inner and outer boundaries
# # manually plotting the data for circle_radial dipoles for N = 1, 5, 10  on log scale
# plt.figure(figsize=(8, 5), dpi=300)
# plt.scatter(x,y1, marker = '.', s = 500, c = 'g', alpha=0.5, label = 'Network 1')    
# plt.scatter(x,y2, marker = '*', s = 50, c = 'b', alpha=0.5, label = 'Network 2')    
# if srand != 112 or pbond == 0.61:
#     plt.scatter(x,y3, marker = '>', s = 50, c = 'r', alpha=0.5, label = 'Network 3')    
# plt.scatter(x,y4, marker = 's', s = 50, c = 'k', alpha=0.5, label = 'Network 4')    
# #if srand == 12:
# #    plt.scatter(x,y5, marker = 'o', s = 50, c = 'y', alpha=0.5, label = 'Network 5')    

# #plt.title('Negative forces', fontsize = 15)
# plt.xlabel("Number of Dipoles", fontsize = 15)
# plt.ylabel("Boundary Force", fontsize = 15)
# plt.xticks(fontsize = 15)
# plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
# plt.legend(loc = 'upper left')
# plt.savefig(plot_folder+"log_log_test_p"+pbond_string+"_"+str(srand)+".png")

# plt.clf()
# plt.close()



############################################################################################################################
# Testing the cartesian stress calculation over individual nodes on the boundary : a test
############################################################################################################################
node1 = 223
# sigma_xx = 0
# sigma_xy = 0
# sigma_yx = 0
# sigma_yy = 0

sigma_xx_list = []
sigma_xy_list = []
sigma_yx_list = []
sigma_yy_list = []
sigma_rr_list = []
angle_val_list = []
strain_val_list = []

strain_valx_list = []
strain_valy_list = []

for i in range(0,len(boundary_nodes)):
    sigma_xx = 0
    sigma_xy = 0
    sigma_yx = 0
    sigma_yy = 0
    for j in range(0,6):
        conn_val = conn_node[boundary_nodes[i],j]
        if conn_val != unconnected_const:    # this will NOT work when L = 128!!!!!
            x1 = xpos[-1,boundary_nodes[i]]
            y1 = ypos[-1,boundary_nodes[i]]
            x2 = xpos[-1,conn_val]
            y2 = ypos[-1,conn_val]
            angle_val = angle(x1,y1,x2,y2)
            # angle_val_list.append(angle_val)
        
            dx = x2-x1
            dy = y2-y1
            
            strain_val = mu*strain[-1,boundary_nodes[i],j]        # technichally its force_val since we multiplied strain by mu
            # strain_val_list.append(strain_val)
            strain_val_x = strain_val*np.cos(angle_val)
            strain_val_y = strain_val*np.sin(angle_val)
            # strain_valx_list.append(strain_val_x)
            # strain_valy_list.append(strain_val_y)
        
            # using the strain forces only
            temp_sigma_xx = strain_val_x*dx
            temp_sigma_xy = strain_val_x*dy
            temp_sigma_yx = strain_val_y*dx
            temp_sigma_yy = strain_val_y*dy
                
            sigma_xx += temp_sigma_xx
            sigma_xy += temp_sigma_xy
            sigma_yx += temp_sigma_yx
            sigma_yy += temp_sigma_yy
    
    for k in range(0,3):  # bending part :: for collinear bonds!!
        conn_val_b1 = conn_node[boundary_nodes[i],k]
        conn_val_b2 = conn_node[boundary_nodes[i],k+3]
        if (conn_val_b1 != unconnected_const and conn_val_b2 != unconnected_const):    # this will NOT work when L = 128!!!!! both co-linear bonds need to be present to calculate the bending contribution        

            # 2nd node is the central node for bending calculations
            x2_bend = xpos[-1,boundary_nodes[i]]
            y2_bend = ypos[-1,boundary_nodes[i]]
            
            x1_bend = xpos[-1,conn_val_b1]
            y1_bend = ypos[-1,conn_val_b1]
            x3_bend = xpos[-1,conn_val_b2]
            y3_bend = ypos[-1,conn_val_b2]

            # radial vector from the center of the simulation. the center of the simulation does not change because the outer boundaries are always fixed
            dx1 = x1_bend-center_x[0]   # 0th is the unperturbed configuration
            dy1 = y1_bend-center_y[0]
            dx2 = x2_bend-center_x[0]
            dy2 = y2_bend-center_y[0]
            dx3 = x3_bend-center_x[0]
            dy3 = y3_bend-center_y[0]
            
            # testing with caluclated bending forces
            f1x = force_bend_x_func_node[conn_val_b1]
            f1y = force_bend_y_func_node[conn_val_b1]
            f2x = force_bend_x_func_node[boundary_nodes[i]]
            f2y = force_bend_y_func_node[boundary_nodes[i]]
            f3x = force_bend_x_func_node[conn_val_b2]
            f3y = force_bend_y_func_node[conn_val_b2]

            # prefactor of 1/3 comes from lammps docementation on stress/atom calculation
            bend_sigma_xx = (1/3.) * (dx1*f1x + dx2*f2x + dx3*f3x)
            bend_sigma_xy = (1/3.) * (dx1*f1y + dx2*f2y + dx3*f3y)
            bend_sigma_yx = (1/3.) * (dy1*f1x + dy2*f2x + dy3*f3x)
            bend_sigma_yy = (1/3.) * (dy1*f1y + dy2*f2y + dy3*f3y)
    
            sigma_xx += bend_sigma_xx
            sigma_xy += bend_sigma_xy
            sigma_yx += bend_sigma_yx
            sigma_yy += bend_sigma_yy

    # dividing by area
    sigma_xx = (1/(np.sqrt(3)*(1/2))) * sigma_xx
    sigma_xy = (1/(np.sqrt(3)*(1/2))) * sigma_xy
    sigma_yx = (1/(np.sqrt(3)*(1/2))) * sigma_yx
    sigma_yy = (1/(np.sqrt(3)*(1/2))) * sigma_yy

    sigma_xx_list.append(sigma_xx)
    sigma_xy_list.append(sigma_xy)
    sigma_yx_list.append(sigma_yx)
    sigma_yy_list.append(sigma_yy)
    
    x1 = xpos[-1,center]
    y1 = ypos[-1,center]
    x2 = xpos[-1,boundary_nodes[i]]
    y2 = ypos[-1,boundary_nodes[i]]
    angle_val = angle(x1,y1,x2,y2)
    
    sigma_rr = sigma_xx * np.square(np.cos(angle_val)) + sigma_yy * np.square(np.sin(angle_val)) + 2.0 * sigma_xy * np.cos(angle_val) * np.sin(angle_val) 
    sigma_rr_list.append(sigma_rr)

sigma_rr_list = np.array(sigma_rr_list)
# plotting radial stress vs radial boundary force
#find line of best fit
slope, intercept = np.polyfit(sigma_rr_list,perp_boundary_force[-1,:], 1)
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(sigma_rr_list,perp_boundary_force[-1,:], marker = 'o', s = 50, c = 'g', alpha=0.2, label = pbond_string)# Network 1')    
plt.plot(sigma_rr_list, slope*sigma_rr_list+intercept, label = '%0.4f *x + %0.2e' % (slope, intercept))
# plt.title('Positive forces', fontsize = 15)
plt.xlabel("Radial Force", fontsize = 15)
plt.ylabel("Radial Stress", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.legend(loc = 'best')
plt.subplots_adjust(top=0.92, left = 0.16, bottom = 0.15, right = 0.96)
plt.savefig(fname_plot+"radial_stress_radial_bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()


# writing the radial stress to file:
force_outfname = base+folder+"txt/displacement/"+"bndry_node_radial_stress_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    force_outfname = base+folder+"txt/displacement/"+"bndry_node_radial_stress_" + 'srand_'+ str(srand) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    force_outfname = base+folder+"txt/displacement/"+"bndry_node_radial_stress_" + 'srand_'+ str(srand) + "_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    force_outfname = base+folder+"txt/displacement/"+"bndry_node_radial_stress_" + str(L) + "_" +str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

#if diff_flag == 1:
#    force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('here all individual radial stresses on BOUNDARY nodes are written to file : ',force_outfname)    
heading = '# nodes       Sigma_rr              Perperndicular F'
fmt = '%12d', '%15.7e', '%15.7e'

np.savetxt(force_outfname, np.column_stack((boundary_nodes, sigma_rr_list, perp_boundary_force[-1])), header = heading, fmt = fmt)



# sigma_xx_list = np.squeeze(sigma_xx_list)
# sigma_xy_list = np.squeeze(sigma_xy_list)
# sigma_yx_list = np.squeeze(sigma_yx_list)
# sigma_yy_list = np.squeeze(sigma_yy_list)

# print(sigma_xx, sigma_xy, sigma_yx, sigma_yy)

# for node 223: sigma xx, xy, yx and yy: -4.24120392213691e-05 -4.973213217198923e-06 -4.9732132171989164e-06 -0.00012726621960954055
# for node 226: sigma xx, xy, yx and yy: -4.241203922136937e-05 4.973213217199147e-06 4.973213217199133e-06 -0.00012726621960954052

############################################################################################################################
# Testing the radial stress calculation over individual nodes on the boundary : a test
############################################################################################################################
# finding the angle from center

# x1 = xpos[-1,center]
# y1 = ypos[-1,center]
# x2 = xpos[-1,node1]
# y2 = ypos[-1,node1]
# angle_val = angle(x1,y1,x2,y2)

# sigma_rr = sigma_xx * np.square(np.cos(angle_val)) + sigma_yy * np.square(np.sin(angle_val)) + 2.0 * sigma_xy * np.cos(angle_val) * np.sin(angle_val) 
# print(sigma_rr)

# # finding the length constant in sigma_rr = (1/L) * radial_force
# loc_node1 = np.squeeze(np.where(boundary_nodes == node1))
# l_val = perp_boundary_force[-1,loc_node1]/sigma_rr
# print(l_val, perp_boundary_force[-1,loc_node1])

############################################################################################################################
# Testing the cartesian stress calculation over individual nodes on the INNER boundary : a test
# defining the nodes on inner boundary
radius = 13
dr = 0.5
inner_bndry_nodes = []
x2 = xpos[0,center]
y2 = ypos[0,center]
for i in range(0,num_pts):
    x1 = xpos[-1,i]
    y1 = ypos[-1,i]

    dx = x2-x1
    dy = y2-y1

    dist = np.sqrt(dx*dx + dy*dy)            

    if ((dist>radius-dr) and (dist<radius+dr)):
        inner_bndry_nodes.append(i)

# writing the inner boundary nodes to file:
inner_bndry_outfname = base+folder+"txt/area/"+"inner_bndry_nodes_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    inner_bndry_outfname = base+folder+"txt/area/"+"inner_bndry_nodes_srand_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    inner_bndry_outfname = base+folder+"txt/area/"+"inner_bndry_nodes_srand_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    inner_bndry_outfname = base+folder+"txt/area/"+"inner_bndry_nodes_srand_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

#if diff_flag == 1:
#    force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('here all inner BOUNDARY nodes are written to file : ',inner_bndry_outfname)    
heading = '# nodes       '
fmt = '%12d' 

np.savetxt(inner_bndry_outfname, inner_bndry_nodes, header = heading, fmt = fmt)


inner_bndry_force_x = force_all[:,inner_bndry_nodes,0]
inner_bndry_force_y = force_all[:,inner_bndry_nodes,1]

inner_bndry_force_vec = []
for i in range(0,np.shape(inner_bndry_force_x)[0]):
    inner_bndry_force_vec.append([])
    for j in range(0,np.shape(inner_bndry_force_x)[1]):
        inner_bndry_force_vec[i].append([inner_bndry_force_x[i,j], inner_bndry_force_y[i,j]])

inner_bndry_force_vec = np.array(inner_bndry_force_vec)


# finding the radial force on inner boundary nodes
center_x = xpos[:,center]
center_y = ypos[:,center]
center_vec = [[center_x[i],center_y[i]] for i in range(0,len(center_x))]

# coordinates of the boundary nodes
inner_bndry_nodes_x = xpos[:,inner_bndry_nodes]
inner_bndry_nodes_y = ypos[:,inner_bndry_nodes]

# radial vector from boundary to center
inner_bndry_radial_vec = []
for i in range(0,np.shape(inner_bndry_nodes_x)[0]):
    inner_bndry_radial_vec.append([])
    for j in range(0,np.shape(inner_bndry_nodes_x)[1]):
        radial_x = center_x[i] - inner_bndry_nodes_x[i,j]
        radial_y = center_y[i] - inner_bndry_nodes_y[i,j]
        
        inner_bndry_radial_vec[i].append([radial_x,radial_y])

inner_bndry_radial_vec = np.array(inner_bndry_radial_vec)

# finding the projection of the boundary forces on the radial vector
temp_boundary_force = []
for i in range(0,np.shape(inner_bndry_force_vec)[0]):
    perp_force = np.dot(inner_bndry_force_vec[i], np.transpose(radial_vec[i+1])) / np.linalg.norm(radial_vec[i+1],axis=1)
    temp_boundary_force.append(perp_force)   

# temp_boundary_force has a size of (num-1,num boundary nodes, num boundary nodes) because the previous operation in not element wise
# so we will extract the diagonal terms from it that is the right calculation
inner_perp_boundary_force = []    
for i in range(0,np.shape(temp_boundary_force)[0]):
    temp_force = temp_boundary_force[i]
    inner_perp_boundary_force.append(temp_force.diagonal())
    
    # perp_boundary_force.append([])
    # for j in range(0,np.shape(temp_boundary_force)[1]):
    #     for k in range(0,np.shape(temp_boundary_force)[2]):
    #         if j == k:
    #             perp_boundary_force[i].append(temp_boundary_force[i,j,k])
    
    
inner_perp_boundary_force = np.array(inner_perp_boundary_force)



# calculating stress
inner_sigma_xx_list = []
inner_sigma_xy_list = []
inner_sigma_yx_list = []
inner_sigma_yy_list = []
inner_sigma_rr_list = []
# inner_angle_val_list = []
# inner_strain_val_list = []

# inner_strain_valx_list = []
# inner_strain_valy_list = []

for i in range(0,len(inner_bndry_nodes)):
    inner_sigma_xx = 0
    inner_sigma_xy = 0
    inner_sigma_yx = 0
    inner_sigma_yy = 0
    for j in range(0,6):
        conn_val = conn_node[inner_bndry_nodes[i],j]
        if (conn_val != unconnected_const):    # this will NOT work when L = 128!!!!!        
            x1 = xpos[-1,inner_bndry_nodes[i]]
            y1 = ypos[-1,inner_bndry_nodes[i]]
            x2 = xpos[-1,conn_val]
            y2 = ypos[-1,conn_val]
            angle_val = angle(x1,y1,x2,y2)
            # angle_val_list.append(angle_val)
        
            dx = x2-x1
            dy = y2-y1
            
            strain_val = mu*strain[-1,inner_bndry_nodes[i],j]        # technichally its force_val since we multiplied strain by mu
            # strain_val_list.append(strain_val)
            strain_val_x = strain_val*np.cos(angle_val)
            strain_val_y = strain_val*np.sin(angle_val)
            # strain_valx_list.append(strain_val_x)
            # strain_valy_list.append(strain_val_y)
        
            temp_sigma_xx = strain_val_x*dx
            temp_sigma_xy = strain_val_x*dy
            temp_sigma_yx = strain_val_y*dx
            temp_sigma_yy = strain_val_y*dy
                
            inner_sigma_xx += temp_sigma_xx
            inner_sigma_xy += temp_sigma_xy
            inner_sigma_yx += temp_sigma_yx
            inner_sigma_yy += temp_sigma_yy
            
    for k in range(0,3):  # bending part :: for collinear bonds!!
        conn_val_b1 = conn_node[inner_bndry_nodes[i],k]
        conn_val_b2 = conn_node[inner_bndry_nodes[i],k+3]
        if (conn_val_b1 != unconnected_const and conn_val_b2 != unconnected_const):    # this will NOT work when L = 128!!!!! both co-linear bonds need to be present to calculate the bending contribution        

            # 2nd node is the central node for bending calculations
            x2_bend = xpos[-1,inner_bndry_nodes[i]]
            y2_bend = ypos[-1,inner_bndry_nodes[i]]
            
            x1_bend = xpos[-1,conn_val_b1]
            y1_bend = ypos[-1,conn_val_b1]
            x3_bend = xpos[-1,conn_val_b2]
            y3_bend = ypos[-1,conn_val_b2]

            # radial vector from the center of the simulation. the center of the simulation does not change because the outer boundaries are always fixed
            dx1 = x1_bend-center_x[0]   # 0th is the unperturbed configuration
            dy1 = y1_bend-center_y[0]
            dx2 = x2_bend-center_x[0]
            dy2 = y2_bend-center_y[0]
            dx3 = x3_bend-center_x[0]
            dy3 = y3_bend-center_y[0]
            
            # testing with caluclated bending forces
            f1x = force_bend_x_func_node[conn_val_b1]
            f1y = force_bend_y_func_node[conn_val_b1]
            f2x = force_bend_x_func_node[inner_bndry_nodes[i]]
            f2y = force_bend_y_func_node[inner_bndry_nodes[i]]
            f3x = force_bend_x_func_node[conn_val_b2]
            f3y = force_bend_y_func_node[conn_val_b2]

            # prefactor of 1/3 comes from lammps docementation on stress/atom calculation
            bend_sigma_xx = (1/3.) * (dx1*f1x + dx2*f2x + dx3*f3x)
            bend_sigma_xy = (1/3.) * (dx1*f1y + dx2*f2y + dx3*f3y)
            bend_sigma_yx = (1/3.) * (dy1*f1x + dy2*f2x + dy3*f3x)
            bend_sigma_yy = (1/3.) * (dy1*f1y + dy2*f2y + dy3*f3y)
    
            inner_sigma_xx += bend_sigma_xx
            inner_sigma_xy += bend_sigma_xy
            inner_sigma_yx += bend_sigma_yx
            inner_sigma_yy += bend_sigma_yy

    # dividing by area
    inner_sigma_xx = (1/(np.sqrt(3)*(1/2))) * inner_sigma_xx
    inner_sigma_xy = (1/(np.sqrt(3)*(1/2))) * inner_sigma_xy
    inner_sigma_yx = (1/(np.sqrt(3)*(1/2))) * inner_sigma_yx
    inner_sigma_yy = (1/(np.sqrt(3)*(1/2))) * inner_sigma_yy

    inner_sigma_xx_list.append(inner_sigma_xx)
    inner_sigma_xy_list.append(inner_sigma_xy)
    inner_sigma_yx_list.append(inner_sigma_yx)
    inner_sigma_yy_list.append(inner_sigma_yy)
    
    x1 = xpos[-1,center]
    y1 = ypos[-1,center]
    x2 = xpos[-1,inner_bndry_nodes[i]]
    y2 = ypos[-1,inner_bndry_nodes[i]]
    angle_val = angle(x1,y1,x2,y2)
    
    inner_sigma_rr = inner_sigma_xx * np.square(np.cos(angle_val)) + inner_sigma_yy * np.square(np.sin(angle_val)) + 2.0 * inner_sigma_xy * np.cos(angle_val) * np.sin(angle_val) 
    inner_sigma_rr_list.append(inner_sigma_rr)

inner_sigma_rr_list = np.array(inner_sigma_rr_list)
# plotting radial stress vs radial boundary force
#find line of best fit
inner_slope, inner_intercept = np.polyfit(inner_sigma_rr_list,inner_perp_boundary_force[-1,:], 1)
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(inner_sigma_rr_list,inner_perp_boundary_force[-1,:], marker = 'o', s = 50, c = 'g', alpha=0.2, label = pbond_string)# Network 1')    
plt.plot(inner_sigma_rr_list, inner_slope*inner_sigma_rr_list+inner_intercept, label = '%0.4f *x + %0.2e' % (inner_slope, inner_intercept))
# plt.title('Positive forces', fontsize = 15)
plt.xlabel("Radial Force", fontsize = 15)
plt.ylabel("Radial Stress", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.legend(loc = 'best')
plt.subplots_adjust(top=0.92, left = 0.16, bottom = 0.15, right = 0.96)
plt.savefig(fname_plot+"inner_radial_stress_radial_bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"__force.png")

plt.clf()
plt.close()

############################################################################################################################
# calculating the total and writing to file
############################################################################################################################
total_inner_stress = np.sum(inner_sigma_rr_list)

inner_stress_outfname = base+folder+"txt/area/"+"total_inner_stress_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    inner_stress_outfname = base+folder+"txt/area/"+"total_inner_stress_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    inner_stress_outfname = base+folder+"txt/area/"+"total_inner_stress_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    inner_stress_outfname = base+folder+"txt/area/"+"total_inner_stress_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

#if diff_flag == 1:
#    force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('here sum of inner stresses on inner boundary nodes are written to file : ',inner_stress_outfname)    
heading = '# nodes       Sigma_rr'
fmt = '%12d', '%15.7e'

np.savetxt(inner_stress_outfname, np.column_stack((len(inner_sigma_rr_list),total_inner_stress)), header = heading, fmt = fmt)


############################################################################################################################
# finding the stress in concentric circles around a central dipole in p = 1 network with N_d = 1
############################################################################################################################
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

# radius_ring_nodes = np.array(radius_ring_nodes)    # nested list of lists has different sizes so cannot be made into a numpy array

#==============================================================================
# calculating stress
ring_sigma_xx_list = []
ring_sigma_xy_list = []
ring_sigma_yx_list = []
ring_sigma_yy_list = []
ring_sigma_rr_list = []
# ring_angle_val_list = []
# ring_strain_val_list = []

# ring_strain_valx_list = []
# ring_strain_valy_list = []

for m in range(0,len(radius_ring_nodes)):
    ring_sigma_xx_list.append([])
    ring_sigma_xy_list.append([])
    ring_sigma_yx_list.append([])
    ring_sigma_yy_list.append([])
    ring_sigma_rr_list.append([])
    
    for i in range(0,len(radius_ring_nodes[m])):
        ring_sigma_xx = 0
        ring_sigma_xy = 0
        ring_sigma_yx = 0
        ring_sigma_yy = 0
        for j in range(0,6):
            conn_val = conn_node[radius_ring_nodes[m][i],j]
            if conn_val != unconnected_const:    # this will NOT work when L = 128!!!!!        
                x1 = xpos[-1,radius_ring_nodes[m][i]]
                y1 = ypos[-1,radius_ring_nodes[m][i]]
                x2 = xpos[-1,conn_val]
                y2 = ypos[-1,conn_val]
                angle_val = angle(x1,y1,x2,y2)
                # angle_val_list.append(angle_val)
            
                dx = x2-x1
                dy = y2-y1
                
                strain_val = mu*strain[-1,radius_ring_nodes[m][i],j]        # technichally its force_val since we multiplied strain by mu
                # strain_val_list.append(strain_val)
                strain_val_x = strain_val*np.cos(angle_val)
                strain_val_y = strain_val*np.sin(angle_val)
                # strain_valx_list.append(strain_val_x)
                # strain_valy_list.append(strain_val_y)
            
                # using the strains to induce forces
                temp_sigma_xx = strain_val_x*dx
                temp_sigma_xy = strain_val_x*dy
                temp_sigma_yx = strain_val_y*dx
                temp_sigma_yy = strain_val_y*dy

                # # testing using total force                        
                # temp_sigma_xx = force_all[-1,radius_ring_nodes[m][i],0]*dx
                # temp_sigma_xy = force_all[-1,radius_ring_nodes[m][i],0]*dy
                # temp_sigma_yx = force_all[-1,radius_ring_nodes[m][i],1]*dx
                # temp_sigma_yy = force_all[-1,radius_ring_nodes[m][i],1]*dy

                ring_sigma_xx += temp_sigma_xx
                ring_sigma_xy += temp_sigma_xy
                ring_sigma_yx += temp_sigma_yx
                ring_sigma_yy += temp_sigma_yy
                
        for k in range(0,3):  # bending part :: for collinear bonds!!
            conn_val_b1 = conn_node[radius_ring_nodes[m][i],k]
            conn_val_b2 = conn_node[radius_ring_nodes[m][i],k+3]
            if (conn_val_b1 != unconnected_const and conn_val_b2 != unconnected_const):    # this will NOT work when L = 128!!!!! both co-linear bonds need to be present to calculate the bending contribution        
    
                # 2nd node is the central node for bending calculations
                x2_bend = xpos[-1,radius_ring_nodes[m][i]]
                y2_bend = ypos[-1,radius_ring_nodes[m][i]]
                
                x1_bend = xpos[-1,conn_val_b1]
                y1_bend = ypos[-1,conn_val_b1]
                x3_bend = xpos[-1,conn_val_b2]
                y3_bend = ypos[-1,conn_val_b2]
    
                # radial vector from the center of the simulation. the center of the simulation does not change because the outer boundaries are always fixed
                dx1 = x1_bend-center_x[0]   # 0th is the unperturbed configuration
                dy1 = y1_bend-center_y[0]
                dx2 = x2_bend-center_x[0]
                dy2 = y2_bend-center_y[0]
                dx3 = x3_bend-center_x[0]
                dy3 = y3_bend-center_y[0]
                
                # testing with caluclated bending forces
                f1x = force_bend_x_func_node[conn_val_b1]
                f1y = force_bend_y_func_node[conn_val_b1]
                f2x = force_bend_x_func_node[radius_ring_nodes[m][i]]
                f2y = force_bend_y_func_node[radius_ring_nodes[m][i]]
                f3x = force_bend_x_func_node[conn_val_b2]
                f3y = force_bend_y_func_node[conn_val_b2]                
    
                # # testing with caluclated bending forces
                # f1x = force_all[-1,conn_val_b1,0]
                # f1y = force_all[-1,conn_val_b1,1]
                # f2x = force_all[-1,radius_ring_nodes[m][i],0]
                # f2y = force_all[-1,radius_ring_nodes[m][i],1]
                # f3x = force_all[-1,conn_val_b2,0]
                # f3y = force_all[-1,conn_val_b2,1]                

                # prefactor of 1/3 comes from lammps docementation on stress/atom calculation
                bend_sigma_xx = (1/3.) * (dx1*f1x + dx2*f2x + dx3*f3x)
                bend_sigma_xy = (1/3.) * (dx1*f1y + dx2*f2y + dx3*f3y)
                bend_sigma_yx = (1/3.) * (dy1*f1x + dy2*f2x + dy3*f3x)
                bend_sigma_yy = (1/3.) * (dy1*f1y + dy2*f2y + dy3*f3y)
        
                ring_sigma_xx += bend_sigma_xx
                ring_sigma_xy += bend_sigma_xy
                ring_sigma_yx += bend_sigma_yx
                ring_sigma_yy += bend_sigma_yy
    
        # dividing by area
        ring_sigma_xx = (1/(np.sqrt(3)*(1/2))) * ring_sigma_xx
        ring_sigma_xy = (1/(np.sqrt(3)*(1/2))) * ring_sigma_xy
        ring_sigma_yx = (1/(np.sqrt(3)*(1/2))) * ring_sigma_yx
        ring_sigma_yy = (1/(np.sqrt(3)*(1/2))) * ring_sigma_yy
    
        ring_sigma_xx_list[m].append(ring_sigma_xx)
        ring_sigma_xy_list[m].append(ring_sigma_xy)
        ring_sigma_yx_list[m].append(ring_sigma_yx)
        ring_sigma_yy_list[m].append(ring_sigma_yy)
        
        x1 = xpos[-1,center]
        y1 = ypos[-1,center]
        x2 = xpos[-1,radius_ring_nodes[m][i]]
        y2 = ypos[-1,radius_ring_nodes[m][i]]
        angle_val = angle(x1,y1,x2,y2)
        
        ring_sigma_rr = ring_sigma_xx * np.square(np.cos(angle_val)) + ring_sigma_yy * np.square(np.sin(angle_val)) + 2.0 * ring_sigma_xy * np.cos(angle_val) * np.sin(angle_val) 
        ring_sigma_rr_list[m].append(ring_sigma_rr)

# ring_sigma_rr_list = np.array(ring_sigma_rr_list)       # nested list of lists has different sizes so cannot be made into a numpy array

mean_ring_sigma = []
std_ring_sigma = []
for i in range(0,len(ring_sigma_rr_list)):
    mean_ring_sigma.append(np.mean(ring_sigma_rr_list[i]))
    std_ring_sigma.append(np.std(ring_sigma_rr_list[i]))

mean_ring_sigma = np.array(mean_ring_sigma)    
std_ring_sigma = np.array(std_ring_sigma)    
#===========================================================================================================================
# writing the stress on each ring to file
stress_ring_outfname = base+folder+"txt/area/"+"stress_ring_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    stress_ring_outfname = base+folder+"txt/area/"+"stress_ring_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    stress_ring_outfname = base+folder+"txt/area/"+"stress_ring_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    stress_ring_outfname = base+folder+"txt/area/"+"stress_ring_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

#if diff_flag == 1:
#    force_outfname = base+folder+"txt/area/"+"bndry_node_force_"+str(dip_ranseed)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('stress on rings of diff radius are written to file : ',stress_ring_outfname)    
heading = 'radius       mean stress       std stress'
fmt = '%12d', '%15.7e', '%15.7e' 

np.savetxt(stress_ring_outfname,  np.column_stack((radius_bins, mean_ring_sigma, std_ring_sigma)), header = heading, fmt = fmt)


# '''
# just commented to speed up things

# plotting the stress on each ring for a N=1 and p=1 case

#decay of stress
xaxis = np.linspace(radius_bins[0],radius_bins[-1],100)
y_fit = mean_ring_sigma[0]/np.power(xaxis,2)

plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(radius_bins, mean_ring_sigma, yerr = std_ring_sigma, linestyle = 'None', marker = '.', ms = 10, c = 'g', alpha=0.5, label = '$p=%s$' % (pbond_string))    
plt.plot(xaxis, y_fit, label = '$1/r^2$')
plt.xlabel("Radial distance", fontsize = 15)
plt.ylabel("Radial Stress", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc = 'best')
if L == 64:
    plt.savefig(base + folder+"png/radial_stress_vs_radial_dist_p_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".png")
if L != 64:
    plt.savefig(base + folder+"png/radial_stress_vs_radial_dist_p_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".png")

plt.clf()
plt.close()

cutoff_fit = np.where(xaxis >= 13)[0][0]
cutoff_data = np.where(radius_bins >= 13)[0][0]
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(radius_bins[cutoff_data:], mean_ring_sigma[cutoff_data:], yerr = std_ring_sigma[cutoff_data:], linestyle = 'None', marker = '.', ms = 10, c = 'g', alpha=0.5, label = '$p=%s$' % (pbond_string))    
plt.plot(xaxis[cutoff_fit:], y_fit[cutoff_fit:], label = '$1/r^2$')
plt.xlabel("Radial distance", fontsize = 15)
plt.ylabel("Radial Stress", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc = 'best')
if L == 64:
    plt.savefig(base + folder+"png/radial_stress_vs_radial_dist_p_outer_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".png")
if L != 64:
    plt.savefig(base + folder+"png/radial_stress_vs_radial_dist_p_outer_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".png")

plt.clf()
plt.close()


cutoff_fit = np.where(xaxis >= 13)[0][0]
cutoff_data = np.where(radius_bins >= 13)[0][0]
plt.figure(figsize=(8, 5), dpi=300)
plt.scatter(radius_bins[cutoff_data:-2], mean_ring_sigma[cutoff_data:-2], linestyle = 'None', marker = '.', s = 400, c = 'g', alpha=0.5, label = '$p=%s$' % (pbond_string))    
# plt.plot(xaxis[cutoff_fit:], y_fit[cutoff_fit:], label = '$1/r^2$')
plt.xlabel("Radial distance", fontsize = 15)
plt.ylabel("Radial Stress", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc = 'best')
if L == 64:
    plt.savefig(base + folder+"png/radial_stress_vs_radial_dist_p_outer_scatter_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".png")
if L != 64:
    plt.savefig(base + folder+"png/radial_stress_vs_radial_dist_p_outer_scatter_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".png")
plt.clf()
plt.close()

# plotting the stress on each ring for a N=1 and p=1 case
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(radius_bins, mean_ring_sigma, yerr = std_ring_sigma, linestyle = 'None', marker = '.', ms = 10, c = 'g', alpha=0.5, label = '$p=%s$' % (pbond_string))    
plt.plot(xaxis, y_fit, label = '$1/r^2$')
plt.xlabel("Radial distance", fontsize = 15)
plt.ylabel("Radial Stress", fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')
plt.legend(loc = 'best')
if L == 64:
    plt.savefig(base + folder+"png/radial_stress_vs_radial_dist_log_p_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".png")
if L != 64:
    plt.savefig(base + folder+"png/radial_stress_vs_radial_dist_log_p_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".png")
plt.clf()
plt.close()

# '''
#===========================================================================================================================



# '''
############################################################################################################################
# finding the mean distance between the dipole centers
############################################################################################################################
xpos_center_list = xpos[0,num_center_list]    
ypos_center_list = ypos[0,num_center_list]    

# Compute pairwise distances
dx = xpos_center_list[:, np.newaxis] - xpos_center_list[np.newaxis, :]
dy = ypos_center_list[:, np.newaxis] - ypos_center_list[np.newaxis, :]
dist_center_matrix = np.sqrt(dx**2 + dy**2)

i_indices, j_indices = np.triu_indices_from(dist_center_matrix, k=1)
# Extract corresponding node pairs and their distances
pairs = np.array(list(zip(num_center_list[i_indices], num_center_list[j_indices])))
pairwise_distances = dist_center_matrix[i_indices, j_indices]

mean_center_dist = np.mean(dist_center_matrix[np.triu_indices_from(dist_center_matrix, k=1)])
std_center_dist = np.std(dist_center_matrix[np.triu_indices_from(dist_center_matrix, k=1)])

# writing all the airwise distances to file
center_dist_outfname = base+cent_dist_folder+"means/"+"center_dist_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    center_dist_outfname = base+cent_dist_folder+"means/"+"center_dist_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    center_dist_outfname = base+cent_dist_folder+"means/"+"center_dist_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    center_dist_outfname = base+cent_dist_folder+"means/"+"center_dist_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

print('writing all the mEAN AND std pairwise distances to  file : ', center_dist_outfname)    
heading = 'center1          center2           dist'
fmt = '%12d', '%12d', '%15.7f'

if pbond !=1 and num_center != 1:
    np.savetxt(center_dist_outfname, np.column_stack((pairs, pairwise_distances)), header = heading, fmt = fmt)


############################################################################################################################
# finding the radius of gyration and convex hull area
############################################################################################################################
from scipy.spatial import ConvexHull

# Example: Replace this with your real data
# green_centers = np.array([[x1, y1], [x2, y2], [x3, y3], ...])
# Make sure it's an Nx2 array where N = number of green motifs
# Each row is (x, y) coordinates of a motif center

def compute_radius_of_gyration(centers):
    center_of_mass = np.mean(centers, axis=0)
    squared_distances = np.sum((centers - center_of_mass) ** 2, axis=1)
    radius_gyration = np.sqrt(np.mean(squared_distances))
    return radius_gyration

def compute_convex_hull_area(centers):
    if len(centers) < 3:
        return 0  # Area is zero if we have fewer than 3 points
    hull = ConvexHull(centers)
    return hull.area

# Example usage:
green_centers = np.column_stack((xpos_center_list,ypos_center_list))

radius_gyration = compute_radius_of_gyration(green_centers)
hull_area = compute_convex_hull_area(green_centers)

print(f"Radius of Gyration: {radius_gyration:.3f}")
print(f"Convex Hull Area: {hull_area:.3f}")


# writing all the airwise distances to file
outfname = base+cent_dist_folder+"means/"+"centers_quant_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if srand_flag == 1 and num_center > 1 and pbond == 1:
if srand_flag == 1 and pbond == 1 and L == 64 or force_flag == 1:
    outfname = base+cent_dist_folder+"means/"+"centers_quant_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond == 1 and L != 64:
    outfname = base+cent_dist_folder+"means/"+"centers_quant_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
elif srand_flag == 1 and pbond != 1 and L != 64:
    outfname = base+cent_dist_folder+"means/"+"centers_quant_"+str(srand)+"_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"

heading = 'radius_gyration    convex hull area'
fmt = '%15.7f', '%15.7f'

if pbond !=1 and num_center != 1:
    np.savetxt(outfname, np.column_stack((radius_gyration, hull_area)), header = heading, fmt = fmt)
