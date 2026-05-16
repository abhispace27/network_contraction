#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 22:49:54 2025
this code reads in all the simulation data for srand 113....122.
does not work for srand = 112 yet

this is for bending dominated networks only

it reads all the num_center, srand and diff_network values of dfar and energies in this order.
Then it finds the mean, median and, std and sem and writes them to a file.

when we use model 1, it also write the mean radial displacement of bending dominated networks to a master file for a specific N_d.

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

L = 16
L = 64
# L = 128
num_pts = L*L

ranseed = '667720,601210'    # -- network 1
    
cluster_new_flag = 1   # for results from the cluster
hex_flag = 0           # for dipoles with hexagonal bonds forced to be present
hex_rand_flag = 0      # for dipoles with hexagonal bonds randomly removed as per p value
cluster_radial_flag = 1   # for dipoles with just radial bonds

# num_center = 1
num_center_list = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
srand_list = [113, 114, 116, 117, 118, 119, 120, 121, 122]                         # 115 is BAD!!!!!

# num_dip = num_center*6
num = 5+1                                                                        # the total number of force steps
num = 10+1

pbond = 1
# pbond = 0.5
pbond = 0.55
# pbond = 0.6


mu = 1
mu_c = 1
tol = 2e-6
tol = 1.0e-6
tol = 1.0e-7

kappa = 1e-6
kappa = 1e-5
# kappa = 1e-4


if kappa == 1e-6:
    kappa_fname = 'kappa2_e-6/'
elif kappa == 1e-5:
    kappa_fname = 'kappa2_e-5/'
elif kappa == 5e-6:
    kappa_fname = 'kappa2_5e-6/'
elif kappa == 5e-5:
    kappa_fname = 'kappa2_5e-5/'
elif kappa == 1e-4:
    kappa_fname = 'kappa2_e-4/'

rlen = 0.9
rlen_txt = "%.4f" % rlen
rlen_txt_ring = "%.9f" % 0.1


def diff_network_func_kappa6(num_center, srand):
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    if L == 128:
        diff_network_list = ['111111,101111', '111111,111111']
        if srand == 116 and num_center == 80:
            diff_network_list = ['111111,101111']
            
    # if (srand == 114 and num_center == 10) or (srand == 120 and num_center == 10):
    #     diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']
    # elif (srand == 122 and num_center == 10):
    #     diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 113 and (num_center == 15 or num_center == 20 or num_center == 25 or num_center == 30 or num_center == 35 or num_center == 40)):
    #     diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,181111', '111111,191111']
    # # elif (srand == 114 and (num_center == 15 or num_center == 20 or num_center == 25 or num_center == 30 or num_center == 35 or num_center == 40)):
    # #     diff_network_list = ['111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 114 and (num_center == 30)):
    #     diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 116 and (num_center == 35 or num_center == 40)):
    #     diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 117 and (num_center >= 15)):
    #     diff_network_list = ['111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    #     # if (num_center == 30 or num_center == 35):
    #     #     diff_network_list = ['111111,111111', '111111,131111', '111111,141111', '111111,161111', '111111,171111', '111111,181111']
    #     # if (num_center == 40):
    #     #     diff_network_list = ['111111,131111', '111111,141111', '111111,171111', '111111,181111']
    #     if (num_center >= 30):
    #         diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 118 and (num_center >= 15)):
    #     diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 119 and (num_center >= 35)):
    #     diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 120 and (num_center >= 15)):
    #     diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']
    #     # if (num_center == 30):
    #     #     diff_network_list = ['111111,101111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']
    #     # if (num_center >= 35):
    #     #     diff_network_list = ['111111,121111', '111111,141111', '111111,161111', '111111,171111', '111111,191111']
    #     if (num_center >= 30):
    #         diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 121 and (num_center >= 15)):
    #     diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,181111', '111111,191111']
    #     if (num_center >= 25):
    #         # diff_network_list = ['111111,101111', '111111,131111', '111111,141111', '111111,161111', '111111,181111']
    #         diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 122 and (num_center == 15 or num_center == 20 or num_center == 25)):     # srand = 122 does not work for N >= 30 !!
    #     diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    # elif (srand == 122 and num_center >= 30):
    #     diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
    if L == 64 and cluster_radial_flag == 1:
        # print('here')
        if (srand == 114 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 117 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #161111 blows up...
        if (srand == 113 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    #-----------------------------------------------------------------------------------------------------

    if L == 64 and hex_flag == 1:
        if (srand == 114 and num_center == 20):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 114 and num_center == 25):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,191111']     #141111 blows up...
        if (srand == 114 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111']     #141111 blows up...
        if (srand == 114 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 116 and num_center == 20):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 116 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 117 and num_center == 10):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 117 and num_center == 20):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 117 and num_center == 25):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 117 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 118 and num_center == 25):
            diff_network_list = ['111111,101111', '111111,111111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 118 and num_center == 30):
            diff_network_list = ['111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 119 and num_center == 15):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']     #141111 blows up...
        if (srand == 119 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 120 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 120 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 120 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 121 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 122 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
    #-----------------------------------------------------------------------------------------------------

    if L == 64 and hex_rand_flag == 1:
        if (srand == 113 and num_center == 20):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 113 and num_center == 30):
            diff_network_list = ['111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 113 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 114 and num_center == 20):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 114 and num_center == 25):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 114 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111']     #141111 blows up...
        if (srand == 114 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111']     #141111 blows up...

        if (srand == 116 and num_center == 25):
            diff_network_list = ['111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 116 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,171111']     #141111 blows up...
        if (srand == 116 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']     #141111 blows up...

        if (srand == 117 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 117 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...

        if (srand == 118 and num_center == 25):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 118 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,121111', '111111,151111', '111111,161111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 118 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 118 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']     #141111 blows up...

        if (srand == 119 and num_center == 25):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 119 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']     #141111 blows up...
        if (srand == 119 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,191111']     #141111 blows up...

        if (srand == 120 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111']     #141111 blows up...

        if (srand == 121 and num_center == 20):
            diff_network_list = ['111111,101111', '111111,111111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 121 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,161111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 121 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 121 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111']     #141111 blows up...

        if (srand == 122 and num_center == 25):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']     #141111 blows up...
        if (srand == 122 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']     #141111 blows up...
        if (srand == 122 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,151111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 122 and num_center == 40):
            diff_network_list = ['111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
            
            
    diff_network_list = np.array(diff_network_list)
    return diff_network_list
    #-----------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------

def diff_network_func_kappa5(num_center, srand):
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
    if L == 64 and cluster_radial_flag == 1:
        # print('here')
        if (srand == 114 and num_center == 30):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 114 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
        if (srand == 116 and num_center == 40):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
    #-----------------------------------------------------------------------------------------------------

            
    diff_network_list = np.array(diff_network_list)
    return diff_network_list
    #-----------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------

def diff_network_func_kappa4(num_center, srand):
    diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
    if L == 64 and cluster_radial_flag == 1:
        # print('here')
        if (srand == 114 and num_center == 35):
            diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,161111', '111111,171111', '111111,191111']     #141111 blows up...
    #     if (srand == 114 and num_center == 40):
    #         diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,141111', '111111,151111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
    #     if (srand == 116 and num_center == 40):
    #         diff_network_list = ['111111,101111', '111111,111111', '111111,121111', '111111,131111', '111111,151111', '111111,161111', '111111,171111', '111111,181111', '111111,191111']     #141111 blows up...
    #-----------------------------------------------------------------------------------------------------

            
    diff_network_list = np.array(diff_network_list)
    return diff_network_list
    #-----------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------

pbond_string = "%.2f" % pbond # to write the filename
tol_str = "%.2e" % tol # to write the filename
kappa_str = "%.2e" % kappa # to write the filename
if kappa == 0:
    kappa_str = "0.00e+00" # to write the filename
    
mu_str = "%.4f" % mu # to write the filename
mu_c_str = "%.4f" % mu_c # to write the filename

base = "/home/abhinav/david/"


en_tot = []
en_stretch = []
en_compress = []
en_elast = []
en_bend = []
en_arp = []
en_ratio = []
colors = []

bndry_force = []
loc_dip_moment = []
far_dip_moment = []
ratio_dip_moment = []
radial_disp = []
# median_radial_disp = []
ring_stress = []

for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    num_dip = num_center*6
    
    bndry_force.append([])
    loc_dip_moment.append([])
    far_dip_moment.append([])
    ratio_dip_moment.append([])
    
    en_tot.append([])
    en_stretch.append([])
    en_compress.append([])
    en_elast.append([])
    en_bend.append([])
    en_arp.append([])
    en_ratio.append([])
    
    colors.append([])
    
    radial_disp.append([])
    # median_radial_disp.append([])
    ring_stress.append([])
    for j in range(0,len(srand_list)):
        srand = srand_list[j]
        if kappa == 1e-6:
            diff_network_list = diff_network_func_kappa6(num_center, srand)
        elif kappa == 1e-5:
            diff_network_list = diff_network_func_kappa5(num_center, srand)
        elif kappa == 1e-4:
            diff_network_list = diff_network_func_kappa4(num_center, srand)
                    
        bndry_force[i].append([])
        loc_dip_moment[i].append([])
        far_dip_moment[i].append([])
        ratio_dip_moment[i].append([])
        
        en_tot[i].append([])
        en_stretch[i].append([])
        en_compress[i].append([])
        en_elast[i].append([])
        en_bend[i].append([])
        en_arp[i].append([])
        en_ratio[i].append([])
        
        colors[i].append([])
        
        radial_disp[i].append([])
        # median_radial_disp[i].append([])
        ring_stress[i].append([])
        
        for k in range(0,len(diff_network_list)):
            diff_network = diff_network_list[k]
            
            folder = 'lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
            
            if cluster_new_flag == 1 and hex_flag == 1 and hex_rand_flag == 0 and cluster_radial_flag == 0:
                folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
                plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/plots/"
                out_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex/'+ kappa_fname +ranseed+"/means/"
            elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 1 and cluster_radial_flag == 0:
                folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
                plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/plots/"
                out_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds/'+ kappa_fname +ranseed+"/means/"
            elif cluster_new_flag == 1 and hex_flag == 0 and hex_rand_flag == 0 and cluster_radial_flag == 1:
                folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/"+str(srand)+"/"+diff_network+"/"
                plot_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/plots/"
                out_folder = 'cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial/'+ kappa_fname +ranseed+"/means/"
                    
            en_fname = base+folder+"energy/strain/"+"Lattice_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"    
            if L != 64:
                en_fname = base+folder+"energy/strain/"+"Lattice_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+"_force.txt"    
            
            en_temp = np.loadtxt(en_fname)
            en_tot[i][j].append(en_temp[-1,0])
            en_stretch[i][j].append(en_temp[-1,1])
            en_compress[i][j].append(en_temp[-1,2])
            en_bend[i][j].append(en_temp[-1,3])
            en_arp[i][j].append(en_temp[-1,4])

            en_elast_temp = en_temp[-1,1] + en_temp[-1,2]
            en_ratio_temp = en_elast_temp/en_temp[-1,3]
            en_ratio[i][j].append(en_ratio_temp)
            en_elast[i][j].append(en_elast_temp)

            if en_ratio_temp >= 1:
                colors[i][j].append('gray')
            else:
                colors[i][j].append('b')
                
                if cluster_radial_flag == 1:
                    # reading displacement values binned per ring in every simulation
                    disp_fname = base+folder+"txt/displacement/"+"bndry_node_radial_disp_mean_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
                    # heading = 'bin mid       Mean Radial dist        Meadian Radial dist        Std Radial dist'
                    disp_data = np.loadtxt(disp_fname)
                    radial_disp[i][j].append(disp_data[:,1])    # this is the mean displacement of all the nodes in a given ring
                    if i == 0 and j == 0 and k == 0: 
                        # print('i, j, k = ', i, j, k)
                        radial_bins = disp_data[:,0]
                
                    fname_ring_stress = base+folder+"txt/area/"+"stress_ring_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
                    data_ring_stress = np.loadtxt(fname_ring_stress)
                    ring_stress[i][j].append(data_ring_stress[:,1])

            force_fname = base+folder+"txt/area/"+"bndry_force_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
            if L != 64:
                force_fname = base+folder+"txt/area/"+"bndry_force_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
                
            bndry_force_temp = np.loadtxt(force_fname)
            bndry_force[i][j].append(bndry_force_temp[-1,1])
            
            #reading local dipole moment
            fname = base+folder+"txt/area/"+"dipole_moment_ratio_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
            if L != 64:
                fname = base+folder+"txt/area/"+"dipole_moment_ratio_"+str(L)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"    
        
            loc_data = np.loadtxt(fname)    
            loc_dip_moment_temp = loc_data[1]
            far_dip_moment_temp = loc_data[2]
            ratio_dip_moment_temp = loc_data[3]
            loc_dip_moment[i][j].append(loc_dip_moment_temp)
            far_dip_moment[i][j].append(far_dip_moment_temp)
            ratio_dip_moment[i][j].append(ratio_dip_moment_temp)
            

# radius_bins = radial_disp[0][0][0]

en_ratio_1d = []
en_elast_1d = []
en_bend_1d = []
en_tot_1d = []

colors_1d = []
all_dom_dip_mom_1d = []
radial_disp_1d = []
radial_stress_1d = []

all_dom_median_dip_mom = []
all_dom_mean_dip_mom = []
all_dom_std_dip_mom = []
all_dom_sem_dip_mom = []

# ratio of elastic and bending energies
median_en_ratio = []
mean_en_ratio = []
std_en_ratio = []
    
# total energies
median_en_tot = []
mean_en_tot = []
std_en_tot = []
sem_en_tot = []

# sum of stretching and compression energies
median_en_elast = []
mean_en_elast = []
std_en_elast = []
sem_en_elast = []

# bending energy
median_en_bend = []
mean_en_bend = []
std_en_bend = []
sem_en_bend = []

# radial displacement
median_radial_disp = []
mean_radial_disp = []
std_radial_disp = []
sem_radial_disp = []

# radial stress
median_radial_stress = []
mean_radial_stress = []
std_radial_stress = []
sem_radial_stress = []

num_sim = []   # number of simulations

for i in range(0,len(num_center_list)):
    num_center = num_center_list[i]
    num_dip = num_center*6
        
    all_dom_dip_mom_1d.append([item for sublist in far_dip_moment[i] for item in sublist])    
    if hex_flag == 1 and num_center_list[i] == 1:
        # # Find the index of the maximum value
        # max_index = all_dom_dip_mom_1d[i].index(max(all_dom_dip_mom_1d[i]))
        # # Remove the element at that index
        # del all_dom_dip_mom_1d[i][max_index]

        # Set percentile threshold (e.g., top 5%)
        percentile = 95
        threshold = np.percentile(all_dom_dip_mom_1d[i], percentile)        
        # Get indices of values above the threshold
        indices_to_remove = [i for i, val in enumerate(all_dom_dip_mom_1d[i]) if val > threshold]
        # Delete from the end to avoid index shifting
        for j in sorted(indices_to_remove, reverse=True):
            del all_dom_dip_mom_1d[i][j]
    
    if hex_flag == 1 and (num_center_list[i] == 2):
        # Find the index of the maximum value
        max_index = all_dom_dip_mom_1d[i].index(max(all_dom_dip_mom_1d[i]))
        # Remove the element at that index
        del all_dom_dip_mom_1d[i][max_index]

    median_dip_mom = np.median(all_dom_dip_mom_1d[i])
    mean_dip_mom = np.mean(all_dom_dip_mom_1d[i])
    std_dip_mom = np.std(all_dom_dip_mom_1d[i])
    all_dom_median_dip_mom.append(median_dip_mom)
    all_dom_mean_dip_mom.append(mean_dip_mom)
    all_dom_std_dip_mom.append(std_dip_mom)
    all_dom_sem_dip_mom.append(std_dip_mom/np.sqrt(len(all_dom_dip_mom_1d[i])))

    en_ratio_1d.append([item for sublist in en_ratio[i] for item in sublist])    
    en_bend_1d.append([item for sublist in en_bend[i] for item in sublist])    
    en_elast_1d.append([item for sublist in en_elast[i] for item in sublist])    
    en_tot_1d.append([item for sublist in en_tot[i] for item in sublist])    

    if hex_flag == 1 and num_center_list[i] == 1:
        # # Remove the element at that index
        # del en_ratio_1d[i][max_index]
        # del en_bend_1d[i][max_index]
        # del en_elast_1d[i][max_index]
        for j in sorted(indices_to_remove, reverse=True):
            del en_ratio_1d[i][j]
            del en_bend_1d[i][j]
            del en_elast_1d[i][j]
            del en_tot_1d[i][j]

    if hex_flag == 1 and num_center_list[i] == 2:
        # Remove the element at that index
        del en_ratio_1d[i][max_index]
        del en_bend_1d[i][max_index]
        del en_elast_1d[i][max_index]
        del en_tot_1d[i][max_index]

    median_en_ratio_temp = np.median(en_ratio_1d[i])
    mean_en_ratio_temp = np.mean(en_ratio_1d[i])
    std_en_ratio_temp = np.std(en_ratio_1d[i])
    median_en_ratio.append(median_en_ratio_temp)
    mean_en_ratio.append(mean_en_ratio_temp)
    std_en_ratio.append(std_en_ratio_temp)

    median_en_tot_temp = np.median(en_tot_1d[i])
    mean_en_tot_temp = np.mean(en_tot_1d[i])
    std_en_tot_temp = np.std(en_tot_1d[i])
    sem_en_tot_temp = std_en_tot_temp/len(en_tot_1d[i])
    median_en_tot.append(median_en_tot_temp)
    mean_en_tot.append(mean_en_tot_temp)
    std_en_tot.append(std_en_tot_temp)
    sem_en_tot.append(sem_en_tot_temp)

    median_en_elast_temp = np.median(en_elast_1d[i])
    mean_en_elast_temp = np.mean(en_elast_1d[i])
    std_en_elast_temp = np.std(en_elast_1d[i])
    sem_en_elast_temp = std_en_elast_temp/len(en_elast_1d[i])
    median_en_elast.append(median_en_elast_temp)
    mean_en_elast.append(mean_en_elast_temp)
    std_en_elast.append(std_en_elast_temp)
    sem_en_elast.append(sem_en_elast_temp)

    median_en_bend_temp = np.median(en_bend_1d[i])
    mean_en_bend_temp = np.mean(en_bend_1d[i])
    std_en_bend_temp = np.std(en_bend_1d[i])
    sem_en_bend_temp = std_en_bend_temp/len(en_bend_1d[i])
    median_en_bend.append(median_en_bend_temp)
    mean_en_bend.append(mean_en_bend_temp)
    std_en_bend.append(std_en_bend_temp)
    sem_en_bend.append(sem_en_bend_temp)
    
    if cluster_radial_flag == 1:
        # radial displacements
        radial_disp_1d.append([item for sublist in radial_disp[i] for item in sublist])    
        median_radial_disp_temp = np.median(radial_disp_1d[i], axis = 0)
        mean_radial_disp_temp = np.mean(radial_disp_1d[i], axis = 0)
        std_radial_disp_temp = np.std(radial_disp_1d[i], axis = 0)
        sem_radial_disp_temp = std_radial_disp_temp/np.sqrt(len(radial_disp_1d[i]))
        median_radial_disp.append(median_radial_disp_temp)
        mean_radial_disp.append(mean_radial_disp_temp)
        std_radial_disp.append(std_radial_disp_temp)
        sem_radial_disp.append(sem_radial_disp_temp)
        
        # writing the means for each dipole number to file
        outfname = base + out_folder + "radial_disp_srand_113_122_N_"+str(num_center)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
        heading = 'radius         mean disp         median disp       std disp         sem disp'
        fmt = '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e'
        np.savetxt(outfname, np.column_stack((radial_bins, mean_radial_disp_temp, median_radial_disp_temp, std_radial_disp_temp, sem_radial_disp_temp)), header = heading, fmt = fmt)

        # radial stress
        radial_stress_1d.append([item for sublist in ring_stress[i] for item in sublist])    
        median_radial_stress_temp = np.median(radial_stress_1d[i], axis = 0)
        mean_radial_stress_temp = np.mean(radial_stress_1d[i], axis = 0)
        std_radial_stress_temp = np.std(radial_stress_1d[i], axis = 0)
        sem_radial_stress_temp = std_radial_stress_temp/np.sqrt(len(radial_stress_1d[i]))
        median_radial_stress.append(median_radial_stress_temp)
        mean_radial_stress.append(mean_radial_stress_temp)
        std_radial_stress.append(std_radial_stress_temp)
        sem_radial_stress.append(sem_radial_stress_temp)

        # writing the means for each dipole number to file
        outfname = base + out_folder + "radial_stress_srand_113_122_N_"+str(num_center)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
        heading = 'radius         mean stress         median stress       std stress         sem stress'
        fmt = '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e'
        np.savetxt(outfname, np.column_stack((radial_bins, mean_radial_stress_temp, median_radial_stress_temp, std_radial_stress_temp, sem_radial_stress_temp)), header = heading, fmt = fmt)

    num_sim.append(len(all_dom_dip_mom_1d[i]))

    # plotting the energy ratios
    colors_1d.append([item for sublist in colors[i] for item in sublist])
    if hex_flag == 1 and num_center_list[i] == 1:
        # # Remove the element at that index
        # del colors_1d[i][max_index]
        for j in sorted(indices_to_remove, reverse=True):
            del colors_1d[i][j]

    if hex_flag == 1 and (num_center_list[i] == 2):
        # Remove the element at that index
        del colors_1d[i][max_index]

    plt.figure(figsize=(8, 5), dpi=300)
    plt.scatter(np.arange(0,len(en_ratio_1d[i])), en_ratio_1d[i], marker = '.', s = 500, c = colors_1d[i], alpha=0.5)#, label = pbond_string)# Network 1')    
    plt.xlabel("Simulation number", fontsize = 15)
    plt.ylabel("$E_{st}/E_{bend}$", fontsize = 15)
    # plt.xticks([0,25,50,75,100], fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.yscale('log')
    plt.subplots_adjust(right=0.96,left=0.18,top=0.94,bottom=0.12)
    plt.legend(loc = 'lower right')
    plt.savefig(plot_folder+"en_ratio_scatter_all_srand_"+str(num_center_list[i])+"_"+str(num_dip)+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
    plt.clf()
    plt.close()


# writing the number of bending dominated cases to file
num_sim_bend_dom = []
if cluster_radial_flag == 1:
    for i in range(0,len(num_center_list)):
        num_sim_bend_dom.append(len(radial_disp_1d[i]))
    outfname = base + out_folder + "bend_dom_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
    heading = 'N_d         # Bend dom Sims'
    fmt = '%12d', '%12d'
    np.savetxt(outfname, np.column_stack((num_center_list, num_sim_bend_dom)), header = heading, fmt = fmt)

#============================================================================================
# writing the dipole moments vs N to an output file
#============================================================================================
outfname = base + out_folder + "mean_dip_mom_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if bending_flag == 0:
#     outfname = base + folder + "all_dom_mean_dip_mom_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('here the alpha_m values are written to file : ',outfname)    
heading = 'N         mean dfar         median dfar       std dfar         sem dfar'
fmt = '%12d', '%15.7e', '%15.7e', '%15.7e', '%15.7e'
np.savetxt(outfname, np.column_stack((num_center_list, all_dom_mean_dip_mom, all_dom_median_dip_mom, all_dom_std_dip_mom, all_dom_sem_dip_mom)), header = heading, fmt = fmt)

#============================================================================================
# writing the mean energy values vs N to an output file
#============================================================================================
outfname = base + out_folder + "mean_energy_vs_N_values_srand_113_122_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
# if bending_flag == 0:
#     outfname = base + folder + "all_dom_mean_energy_vs_N_values_"+str(srand)+"_"+str(num_center)+"_"+str(num_dip)+"_"+pbond_string+"_"+tol_str+"_"+kappa_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+"_"+str(num-1)+".txt"
print('writing energy vs N to file : ',outfname)    
heading = '        N          #sim   mean el en      median el en    std el en       mean bend en    median bend en   std bend en    mean en ratio   median en ratio  std en ratio'
fmt = '%12d', '%12d', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e', '%15.7e'
np.savetxt(outfname, np.column_stack((num_center_list, num_sim, mean_en_elast, median_en_elast, std_en_elast, mean_en_bend, median_en_bend, std_en_bend, mean_en_ratio, median_en_ratio, std_en_ratio)), header = heading, fmt = fmt)


#========================================================================================================================================================================================
# plotting dfar vs energy ratios
#========================================================================================================================================================================================

# plotting each Nd case separately
for i in range(0,len(num_center_list)):
    plt.figure(figsize=(8, 5), dpi=300)
    plt.scatter(en_ratio_1d[i], all_dom_dip_mom_1d[i], marker = '.', s = 500, alpha=0.5, c = 'b')#, label = pbond_string)# Network 1')    
    plt.xlabel("$E_{st}/E_{bend}$", fontsize = 20)
    plt.ylabel("$D_{far}$", fontsize = 20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.subplots_adjust(right=0.98,left=0.14,top=0.98,bottom=0.14)
    # plt.legend(loc = 'lower right')
    plt.savefig(plot_folder+"dfar_vs_e_ratio_scatter_all_srand_"+str(num_center_list[i])+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
    plt.clf()
    plt.close()

    if hex_flag == 1:
        plt.figure(figsize=(8, 5), dpi=300)
        plt.scatter(en_elast_1d[i], all_dom_dip_mom_1d[i], marker = '.', s = 500, alpha=0.5, c = 'b')#, label = pbond_string)# Network 1')    
        plt.xlabel("$E_{st}$", fontsize = 20)
        plt.ylabel("$D_{far}$", fontsize = 20)
        plt.xscale('log')
        plt.yscale('log')
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)
        plt.subplots_adjust(right=0.98,left=0.14,top=0.98,bottom=0.14)
        # plt.legend(loc = 'lower right')
        plt.savefig(plot_folder+"dfar_vs_e_st_scatter_all_srand_"+str(num_center_list[i])+"_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
        plt.clf()
        plt.close()

# plotting all Nd cases together

colors_list = [
    "#1f77b4",  # muted blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # yellow-green
    "#17becf"   # cyan
]

markers_list = [
    "o",  # circle
    "s",  # square
    "^",  # triangle_up
    "v",  # triangle_down
    "D",  # diamond
    "X",  # x-filled
    "*",  # star
    "P",  # plus-filled
    "d",  # plus
    "h"   # x
]

plt.figure(figsize=(8, 5), dpi=300)
for i in range(0,len(num_center_list)):
    plt.scatter(en_ratio_1d[i], all_dom_dip_mom_1d[i], marker = markers_list[i], s = 100, alpha=0.4, c = colors_list[i], label = '$N_d = %d $' % num_center_list[i])# Network 1')    

if cluster_radial_flag == 1:
    if kappa == 1e-4:
        x_linspace = np.logspace(-3,-1.5,100) 
        line_fit1 = (np.power(x_linspace,2)/x_linspace[0])
        plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
        plt.text(1e-3, 1e-1, "$\\sim E_{ratio}^{2}$", fontsize = 15)

        x_linspace2 = np.logspace(-1.5,0,100) 
        line_fit2 = (0.05*np.power(x_linspace2,0.5)/x_linspace2[0])
        plt.plot(x_linspace2, line_fit2, c = 'b', ls = 'dashed') 
        plt.text(1e-1, 1, "$\\sim E_{ratio}^{0.5}$", fontsize = 15)

    elif kappa == 1e-5:
        x_linspace = np.logspace(-4,-2.5,100) 
        line_fit1 = (np.power(x_linspace,2)/x_linspace[0])
        plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
        plt.text(1e-4, 1e-2, "$\\sim E_{ratio}^{2}$", fontsize = 15)

        x_linspace2 = np.logspace(-2,1,100) 
        line_fit2 = (0.005*np.power(x_linspace2,0.5)/x_linspace2[0])
        plt.plot(x_linspace2, line_fit2, c = 'b', ls = 'dashed') 
        plt.text(1e-1, 1, "$\\sim E_{ratio}^{0.5}$", fontsize = 15)

    else:
        x_linspace = np.logspace(-5,-4,100) 
        line_fit1 = (np.power(x_linspace,2)/x_linspace[0])
        plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
        plt.text(6e-6, 5e-4, "$\\sim E_{ratio}^{2}$", fontsize = 15)

        x_linspace2 = np.logspace(-3,2,100) 
        line_fit2 = (0.0003*np.power(x_linspace2,0.5)/x_linspace2[0])
        plt.plot(x_linspace2, line_fit2, c = 'b', ls = 'dashed') 
        plt.text(1e-1, 5e-1, "$\\sim E_{ratio}^{0.5}$", fontsize = 15)

elif hex_rand_flag == 1:
    x_linspace = np.logspace(-4,-2,100) 
    line_fit1 = (np.power(x_linspace,1)/x_linspace[0])*1e-6
    plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
    plt.text(1e-3, 1e-6, "$\\sim E_{ratio}^{1}$", fontsize = 15)

    plt.ylim([1e-7,1e3])

plt.xlabel("$E_{st}/E_{bend}$", fontsize = 20)
plt.ylabel("$D_{far}$", fontsize = 20)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.14,top=0.97,bottom=0.14)
plt.legend(loc = 'best', fontsize = 12)
plt.savefig(plot_folder+"dfar_vs_e_ratio_scatter_all_srand_all_Nd_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()

# dfar vs bending energy
plt.figure(figsize=(8, 5), dpi=300)
for i in range(0,len(num_center_list)):
    plt.scatter(en_bend_1d[i], all_dom_dip_mom_1d[i], marker = markers_list[i], s = 100, alpha=0.4, c = colors_list[i], label = '$N_d = %d $' % num_center_list[i])# Network 1')    

if cluster_radial_flag == 1:
    if kappa == 1e-4:
        x_linspace = np.logspace(-5,-3,100) 
        line_fit1 = (np.power(x_linspace,1)/x_linspace[0])*0.003
        plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
        plt.text(2e-5, 1e-1, "$\\sim E_{bend}^{1}$", fontsize = 20)
    elif kappa == 1e-5:
        x_linspace = np.logspace(-6,-4,100) 
        line_fit1 = (np.power(x_linspace,1)/x_linspace[0])*0.0003
        plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
        plt.text(5e-6, 1e-2, "$\\sim E_{bend}^{1}$", fontsize = 20)
        
    else:        
        x_linspace = np.logspace(-7,-5,100) 
        line_fit1 = (np.power(x_linspace,1)/x_linspace[0])*0.00005
        plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
        plt.text(2e-7, 1e-3, "$\\sim E_{bend}^{1}$", fontsize = 20)

    # line_fit2 = (np.power(x_linspace,0.8)/x_linspace[0])*2e-6
    # plt.plot(x_linspace, line_fit2, c = 'r', ls = 'dashed') 
    # plt.text(1e-7, 5e-2, "$\\sim E_{bend}^{0.5}$", fontsize = 20)

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("$E_{bend}$", fontsize = 20)
plt.ylabel("$D_{far}$", fontsize = 20)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.15,top=0.97,bottom=0.16)
plt.legend(loc = 'best', fontsize = 12)
plt.savefig(plot_folder+"dfar_vs_e_bend_scatter_all_srand_all_Nd_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# dfar vs stret energy
plt.figure(figsize=(8, 5), dpi=300)
for i in range(0,len(num_center_list)):
    plt.scatter(en_elast_1d[i], all_dom_dip_mom_1d[i], marker = markers_list[i], s = 100, alpha=0.4, c = colors_list[i], label = '$N_d = %d $' % num_center_list[i])# Network 1')    

if cluster_radial_flag == 1:
    if kappa == 1e-4:
        x_linspace = np.logspace(-8,-2,100) 
        line_fit1 = np.power(x_linspace,0.5)*20
        plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
        plt.text(1e-7, 5e-2, "$\\sim E_{st}^{0.5}$", fontsize = 20)
    else:
        x_linspace = np.logspace(-11,-3,100) 
        line_fit1 = (np.power(x_linspace,0.5)/x_linspace[0])*2e-10
        plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
        plt.text(1e-7, 5e-2, "$\\sim E_{st}^{0.5}$", fontsize = 20)

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("$E_{st}$", fontsize = 20)
plt.ylabel("$D_{far}$", fontsize = 20)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.15,top=0.97,bottom=0.16)
plt.legend(loc = 'best', fontsize = 12)
plt.savefig(plot_folder+"dfar_vs_e_st_scatter_all_srand_all_Nd_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# dfar vs stret energy
plt.figure(figsize=(8, 5), dpi=300)
for i in range(0,len(num_center_list)):
    plt.scatter(en_tot_1d[i], all_dom_dip_mom_1d[i], marker = markers_list[i], s = 100, alpha=0.4, c = colors_list[i], label = '$N_d = %d $' % num_center_list[i])# Network 1')    

if cluster_radial_flag == 1:
    if kappa == 1e-4:
        x_linspace = np.logspace(-5,-2,100) 
        line_fit1 = (np.power(x_linspace,1)/x_linspace[0])*2e-4
        plt.text(1e-4, 8e-4, "$\\sim E_{tot}^{1}$", fontsize = 20)
    else:
        x_linspace = np.logspace(-7,-2,100) 
        line_fit1 = (np.power(x_linspace,1)/x_linspace[0])*2e-6
        plt.text(1e-4, 8e-4, "$\\sim E_{tot}^{1}$", fontsize = 20)
        
    plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks

plt.xlabel("$E_{tot}$", fontsize = 20)
plt.ylabel("$D_{far}$", fontsize = 20)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.15,top=0.97,bottom=0.16)
plt.legend(loc = 'best', fontsize = 12)
plt.savefig(plot_folder+"dfar_vs_e_tot_scatter_all_srand_all_Nd_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# plotting <E_st> vs Nd
plt.figure(figsize=(8, 5), dpi=300)
plt.errorbar(num_center_list, mean_en_elast, yerr = sem_en_elast, marker = 'o', ms = 15, alpha=1, c = 'g', ls = 'None', label = '$<E_{st}>$')# Network 1')    
plt.errorbar(num_center_list, mean_en_bend, yerr = sem_en_bend, marker = 'd', ms = 15, alpha=1, c = 'b', ls = 'None', label = '$<E_{bend}>$')# Network 1')    
plt.errorbar(num_center_list, mean_en_tot, yerr = sem_en_tot, marker = 's', ms = 15, alpha=0.3, c = 'r', ls = 'None', label = '$<E_{tot}>$')# Network 1')    
x_linspace = np.logspace(0,1,100) 
if kappa == 1e-4:
    line_fit1 = np.power(x_linspace,1)*1e-5
    plt.text(2, 1e-5, "$\\sim E_{bend}^{1}$", fontsize = 15)
    line_fit2 = np.power(x_linspace,1)*5e-7
    plt.text(3, 1e-6, "$\\sim E_{st}^{1}$", fontsize = 15)
elif kappa == 1e-5:
    line_fit1 = np.power(x_linspace,1)*1e-6
    plt.text(2, 5e-6, "$\\sim E_{bend}^{1}$", fontsize = 15)
    line_fit2 = np.power(x_linspace,1)*2e-9
    plt.text(2, 2e-9, "$\\sim E_{st}^{1}$", fontsize = 15)
else:
    line_fit1 = np.power(x_linspace,1)*1e-6
    plt.text(2, 5e-6, "$\\sim E_{bend}^{1}$", fontsize = 15)
    line_fit2 = np.power(x_linspace,1)*2e-10
    plt.text(2, 1e-9, "$\\sim E_{st}^{1}$", fontsize = 15)
    
plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
plt.plot(x_linspace, line_fit2, c = 'g', ls = 'dashed') 

ax = plt.gca()
# Increase tick width and length for major and minor ticks
ax.tick_params(axis='both', which='major', width=3, length=10)  # Major ticks
ax.tick_params(axis='both', which='minor', width=2, length=5)   # Minor ticks
plt.xlabel("Number of dipoles, $N_d$", fontsize = 20)
plt.ylabel("Energy", fontsize = 20)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.17,top=0.97,bottom=0.16)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(plot_folder+"e_st_and_e_bend_vs_Nd_scatter_all_srand_all_Nd_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()



if hex_flag == 1:
    plt.figure(figsize=(8, 5), dpi=300)
    for i in range(0,len(num_center_list)):
        plt.scatter(en_elast_1d[i], all_dom_dip_mom_1d[i], marker = markers_list[i], s = 100, alpha=0.4, c = colors_list[i], label = '$N_d = %d $' % num_center_list[i])# Network 1')            
    plt.xlabel("$E_{st}$", fontsize = 20)
    plt.ylabel("$D_{far}$", fontsize = 20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.subplots_adjust(right=0.98,left=0.14,top=0.97,bottom=0.14)
    plt.legend(loc = 'best', fontsize = 12)
    plt.savefig(plot_folder+"dfar_vs_e_st_scatter_all_srand_all_Nd_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
    plt.clf()
    plt.close()


# normalized version of dfar vs energy ratio:: normalized dfar only
plt.figure(figsize=(8, 5), dpi=300)
for i in range(0,len(num_center_list)):
    # norm_value = np.percentile(all_dom_dip_mom_1d[i], 5)
    norm_value = np.min(all_dom_dip_mom_1d[i])
    normalized_data = all_dom_dip_mom_1d[i] / norm_value
    # normalized_data = all_dom_dip_mom_1d[i] / np.median(all_dom_dip_mom_1d[i])
    # Option 1: Unique values only
    # third_lowest_unique = sorted(set(all_dom_dip_mom_1d[i]))[2]  # 0-based indexing
    # normalized_data = all_dom_dip_mom_1d[i] / third_lowest_unique
    
    plt.scatter(en_ratio_1d[i], normalized_data, marker = markers_list[i], s = 100, alpha=0.4, c = colors_list[i], label = '$N_d = %d $' % num_center_list[i])# Network 1')    
plt.xlabel("$E_{st}/E_{bend}$", fontsize = 20)
plt.ylabel("Normalized $D_{far}$", fontsize = 20)
# plt.xticks([0,25,50,75,100], fontsize = 15)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-1, 5e3])
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.14,top=0.98,bottom=0.14)
plt.legend(loc = 'best', fontsize = 12)
plt.savefig(plot_folder+"dfar_vs_e_ratio_scatter_all_srand_all_Nd_norm_dfar_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()


# normalized version of dfar vs energy ratio:: normalized dfar and energy ratio, both
plt.figure(figsize=(8, 5), dpi=300)
for i in range(0,len(num_center_list)):
    # norm_value = np.percentile(all_dom_dip_mom_1d[i], 5)
    norm_value = np.min(all_dom_dip_mom_1d[i])
    normalized_data = all_dom_dip_mom_1d[i] / norm_value
    # normalized_data = all_dom_dip_mom_1d[i] / np.median(all_dom_dip_mom_1d[i])
    # Option 1: Unique values only
    # third_lowest_unique = sorted(set(all_dom_dip_mom_1d[i]))[2]  # 0-based indexing
    # normalized_data = all_dom_dip_mom_1d[i] / third_lowest_unique
    norm_xvalue = np.min(en_ratio_1d[i])
    normalized_xdata = en_ratio_1d[i] / norm_value
    
    plt.scatter(normalized_xdata, normalized_data, marker = markers_list[i], s = 100, alpha=0.4, c = colors_list[i], label = '$N_d = %d $' % num_center_list[i])# Network 1')    

if cluster_radial_flag == 1:
    x_linspace = np.logspace(0,5,100) 
    line_fit1 = normalized_data[0]*np.power(x_linspace,0.5)/normalized_xdata[0]
    plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
    plt.text(1e4, 5e1, "$\\sim E_{ratio}^{0.5}$", fontsize = 20)
    # line_fit2 = normalized_data[0]*np.power(x_linspace,1)/normalized_xdata[0]
    # plt.plot(x_linspace, line_fit2, c = 'r', ls = 'dashed') 
    # plt.text(1e2, 1e3, "$\\sim E_{ratio}^{1}$", fontsize = 20)
    # line_fit3 = normalized_data[0]*np.power(x_linspace,0.67)/normalized_xdata[0]
    # plt.plot(x_linspace, line_fit3, c = 'g', ls = 'dashed') 
    # plt.text(1e2, 1e2, "$\\sim E_{ratio}^{1}$", fontsize = 20)
    plt.ylim([1e-1, 5e3])

if cluster_radial_flag == 1:
    # x_linspace = np.logspace(0,5,100) 
    # line_fit1 = normalized_data[0]*np.power(x_linspace,0.5)/normalized_xdata[0]
    # plt.plot(x_linspace, line_fit1, c = 'b', ls = 'dashed') 
    # plt.text(1e4, 5e1, "$\\sim E_{ratio}^{0.5}$", fontsize = 20)
    plt.xlim([1, 3e7])

plt.xlabel("Normalized $E_{st}/E_{bend}$", fontsize = 20)
plt.ylabel("Normalized $D_{far}$", fontsize = 20)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplots_adjust(right=0.98,left=0.14,top=0.98,bottom=0.14)
plt.legend(loc = 'best', fontsize = 12)
plt.savefig(plot_folder+"dfar_vs_e_ratio_scatter_all_srand_all_Nd_norm_dfar_Eratio_"+pbond_string+"_"+kappa_str+"_"+tol_str+"_"+rlen_txt+"_"+mu_str+"_"+mu_c_str+".png")
plt.clf()
plt.close()
