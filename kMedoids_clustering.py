#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 13:59:14 2020

@author: Lauri Neuvonen based on template by debbora
"""
# %% IMPORTING NECESSARY PACKAGES 

import numpy as np
from scipy import stats
import pandas as pd

# %% DEFINING FUNCTIONS

# ---------------------- CALCULATING DISTANCES --------------------------------

# 1) - Calculate Pearson Correlation between all policies (in a Pareto front) using the full policy vector (time, strength),
# saved in form ...
#    - Calculate Pearson distance out of given correlation:
#    d(x,y) = sqrt(0.5*(1-corr(x,y))) 
def CalcPearson(data_in):
    data = data_in.copy()
    
    corr = []
    dist = []
    [num_pol, n_t] = data.shape  # 1 time series per policy
    #print("num_pol: ", num_pol)
    
    #going thorugh all cell combinations
    for pol in range(0, num_pol):
        X = data[pol, :]   # get timeseries of cell 1
        corr_tmp = np.empty(num_pol)
        corr_tmp.fill(np.nan)
        dist_tmp = np.empty(num_pol)
        dist_tmp.fill(np.nan)
        for pol2 in range(0, num_pol):

            Y = data[pol2, :]  # timeseries cell 2
            use = np.logical_and(~np.isnan(X), ~np.isnan(Y)) # filters nan-containing data-pointsd
            # from use
            if np.sum(use) > 1: # if at least 2 points are available
                                # the distance is calculated
                corr_tmp[pol2] = stats.pearsonr(X[use], \
                                                         Y[use])[0]
                dist_tmp[pol2] = np.sqrt(0.5*(1 - \
                                              corr_tmp[pol2]))
                #print("corr_tmp: ", corr_tmp)
        corr.append(corr_tmp)
        dist.append(dist_tmp)
        #print("dist", dist)
        #print("len dist: ", len(dist))
        #print("dist 00: ", dist[0][0])

    return corr, dist

# --------------------------- k-Medoids Algorithm -----------------------------


# Definition of the k-Medoids algorithm:
# Step 1. k different objects are chosen as initial medoids by a greedy 
#         algorithm 
# Step 2. Each remaining object is associated with the medoid that is closest. 
#         The total costs are calculated as the sum over all squared distances 
#         between object and respective medoid.
# Step 3. For each pair of object and medoid, it is checked whether switching 
#         them (i.e. the normal object becoming medoid) would improve the 
#         clustering (i.e. decrease the total costs). After going through all
#         pairs, the switch yielding the biggest improvement is performed 
#         and step 3 is repeated. If none of the switches would yield an 
#         improvement the algorithm terminates.
    
# 0) Main part
def kMedoids(k, dist, start_medoids = None):
    num_pol = len(dist)
    terminated = False
    step = 0
    # Normally, the initial medoids are chosen thorugh a greedy algorithm, 
    # but if we wish to continue a run that had not yet terminated or for 
    # some other reason want to start with specific medoids these can be 
    # given to the function and will be used
    if start_medoids == None:
        medoids = GetInitialMedoids(k, dist, num_pol)
    else:
        medoids = start_medoids
    # get best cluster for these medoids and save the results of this step
    cluster, cost = GetCluster(k, dist, num_pol, medoids)
    # Now we iterated until either the maximal number of iteration is reached 
    # or no improvment can be found
    while terminated == False:
        # print(step)
        # going trough all possible switches and finding the best one
        new_cluster, new_cost, new_medoids = GetSwitch(k, dist, num_pol, medoids)
        # if the best possible switch actually improves the clusters we do the
        # switch and save the outcome of this step
        if new_cost < cost:
            cost, cluster, medoids = new_cost, new_cluster, new_medoids
            step += 1                         
            continue
        # if best switch is not improving the clustering, we print a 
        # corresponding message and terminate the algorithm
        # print("No improvement found")
        terminated = True
    return(cluster, medoids, cost)
    
# 1) Greedy algorithm for initial medoids
def GetInitialMedoids(k, dist, num_pol):
    medoids = []
    # for each medoid take the cell the is the is the best choice at this point
    # given the other mediods that are already chosen. 
    for l in range(1, k+1):
        best_cost = 0
        # for each new medoid we check each possible cell 
        for i in range(0, num_pol):

            # temporarily add this cell to medoids
            medoids_tmp = medoids.copy()
            medoids_tmp.append(i)
            # get new costs
            cluster, cost = \
                GetCluster(l, dist, num_pol, medoids_tmp)
            # if best choice so far (or first cell checked) remember
            if best_cost == 0:
                best_cost = cost
                best_medoid = i
            elif cost < best_cost:
                best_cost = cost
                best_medoid = i
        # add best choice to medoids
        medoids.append(best_medoid)
    return(medoids)

# 2) Subroutine to get clusters to given medoids:
def GetCluster(k, dist, num_pol, medoids):
    cluster = np.empty([num_pol])
    cluster.fill(np.nan)
    cl_dist = np.empty([num_pol])
    cl_dist.fill(np.nan)
    # loop over all grid cells
    for i in range(0, num_pol):

        # we index cluster by the position of the respective medoid
        # a medoid obviously belongs to its own cluster
        if i in medoids:
            cluster[i] = medoids.index(i)
            cl_dist[i] = 0
            continue
        # initilizing the best distance with 2,  as 1 is the worst possible
        # distance and we want something worse
        best_dist = 2
        # We check for each medoid how big the distance of the current grid
        # cell to that medoid is. If we found a better distance than the
        # current best_dist we update it and remember the medoid index
        for k in medoids:
            dist_tmp = dist[i][k]
            #print("dist_tmp: ", dist_tmp)
            if dist_tmp < best_dist:
                best_dist = dist_tmp
                best_med = medoids.index(k)
        # we then save the best distance and the corresponding cluster
        cluster[i] = best_med
        cl_dist[i] = best_dist
    # calculating the cost function: sum of all squared distances
    cost = np.nansum(cl_dist**2)
    return(cluster, cost)
                
# 3) Subroutine to get best change in medoids
def GetSwitch(k, dist, num_pol, medoids):
    new_cost = -1
    # loop over all grid cells
    for i in range(0, num_pol):

        # if the grid cell already is a cluster a switch makes no sense
        if i in medoids:
            continue
        # for each of the medoids we check what a switch would result in
        for k in medoids:
            # switching the medoids
            medoids_tmp = medoids[:]
            medoids_tmp[medoids_tmp.index(k)] = i
            # getting the new cluster
            cluster_tmp, cost_tmp = GetCluster(k, dist, num_pol, medoids_tmp)
            # updating if we found a better switch (or if this was the
            # first we tried)
            if cost_tmp < new_cost or new_cost == -1:
                new_cluster = cluster_tmp
                new_cost = cost_tmp
                new_medoids = medoids_tmp
    # returning best switch found (even if it is no improvement to current 
    # situation - this is checked after)
    return(new_cluster, new_cost, new_medoids)
