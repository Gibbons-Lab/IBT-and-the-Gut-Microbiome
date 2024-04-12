# -*- coding: utf-8 -*-
from os.path import expanduser
import sys
import numpy as np
import random
from multiprocessing import Pool
import time
import statsmodels.tsa.stattools as sta

mydir = expanduser("/proj/gibbons/kramos/github/IBT-and-the-Gut-Microbiome/IBM/")
GenPath = mydir + 'data/'

def headings():
    headings = 'sim,'
    headings += 'ct,'
    headings += 'Length,'
    headings += 'TotalAbundance,'
    headings += 'SpeciesRichness,'
    headings += 'Simpson,'
    headings += 'Simpson_e'

    return headings

def simpson_div(rad):
    
    ''' Simpson's diversity:
                D = 1 / Σni(ni-1) / N(N-1)
        Based on 1/D instead of 1 - D, since Simpson's evenness uses 1/D
    '''
    D = 0.0
    N = sum(rad)
    
    for x in rad: 
        D += (x/N)**2
    
    if D <= 0:
        print('Error:')
        print(rad)
        
    D = 1.0/D
    return D, D/len(rad) # this latter value is simpson's evenness
    

def run_model(sim, lgp, imr, rpr, length):
    start_time=time.time()
    seed = None
    N=0
    iteration = 0 #counter for each simulation
    xcoords=np.array([])
    community=np.array([])
    # Ns = []
    # Ss = []
    # Nstationary = 'No'
    # Sstationary = 'No'
    
    while iteration <= 10*length:  # end the simulation if the community has been maximized for more than 10 cycles
        processes=['immigration', 'reproduction', 'flow', 'emigration']
        np.random.RandomState(seed).shuffle(processes)

        for p in processes:
            if p=='flow':
                flow_start=time.time()
                xcoords[:] = xcoords+1
                flow_end=time.time()
                print("Flow time = {:.3f} seconds".format(flow_end - flow_start))

            if p=='immigration':
                im_start=time.time()
            # Immigration from a log-series distributed regional pool. The randomly drawn numbers are 
            # the species IDs. On average, there will be more 1's than 2's, more 2's than 3's, etc.
    
                immigrants = np.array(np.random.RandomState(seed).logseries(lgp, size=imr))
                community=np.concatenate((community,immigrants)) # Add the immigrants to the local community
                xcoords=np.concatenate((xcoords, np.array([0]*len(immigrants)))) # Make the immigrants enter at the inlet
        
                im_end=time.time()
                print("Immigration time = {:.3f} seconds".format(im_end- im_start))

            if p =='emigration':
                em_start = time.time()
                indices = xcoords <= length
                community = community[indices]
                xcoords = xcoords[indices]
                # community = np.delete(community, indices)
                # xcoords = np.delete(xcoords, indices)
                # indices = [i for i, v in enumerate(xcoords) if v > length] # find every individual that flowed out of bounds

                em_end = time.time()
                print("Emigration time = {:.3f} seconds".format(em_end - em_start))

            if p=='reproduction':
                rep_start = time.time()
                n = int(round(rpr * len(community))) # number of bouncing baby Bacilli 
                N = int(len(community))

                if n > 0:
                    # indices = np.random.default_rng(seed).choice(list(range(N)), size=n, replace=False) # choosing mothers at random
                    indices = random.sample(list(range(N)), n)
                    bbb = np.array(community)[indices] # Filter out non-selected indices
                    bbb_x = np.array(xcoords)[indices] # Filter out non-selected indices
                    community=np.concatenate((community, bbb))
                    xcoords=np.concatenate((xcoords, bbb_x))

                rep_end=time.time()
                print("Reproduction time = {:.3f} seconds".format(rep_end - rep_start))
            
            # if p=='death':
            #     n = int(round(dtr * len(community))) # number of members in the community that will die
            #     N = int(len(community))

            #     if n > 0:
            #         indices = np.random.RandomState(seed).choice(list(range(N)), size=n, replace=False) # choosing deaths at random
            #         community = np.delete(community, indices)
            #         xcoords = np.delete(xcoords, indices)
        
        iteration += 1
        
        N = len(community)
        # Ns.append(N)
        S=len(list(set(community)))
        # Ss.append(S)

        if iteration%100 == 0:
                it_start=time.time()
                unique, counts = np.unique(community, return_counts=True)
                sp_d = dict(zip(unique, counts))
                rad=list(sp_d.values())
                rad = sorted(rad)
                div, ev = simpson_div(rad)
                
                # Per sim, record abundance, species richness, simpson, evenness
                OUT = open(GenPath + 'simData.csv', 'a')
                outlist = [sim, iteration, length, N, S, div, ev]
                outlist = str(outlist).strip('[]')
                outlist = outlist.replace(" ", "")
                print(outlist, file=OUT)
                
                print('Instance:', sim, '|   N:', N, 'S:', S)
                it_end=time.time()
                print("Printing sim data time = {:.3f} seconds".format(it_end - it_start))

        # if len(Ns) > 1000:
        #     AugmentedDickeyFuller = sta.adfuller(Ns)
        #     val, pN = AugmentedDickeyFuller[0:2]
                
        #     AugmentedDickeyFuller = sta.adfuller(Ss)
        #     val, pS = AugmentedDickeyFuller[0:2]

        #     if pN > 0.05: 
        #         Nstationary = 'No'

        #     elif pN <= 0.05:
        #         Nstationary = 'Yes'
                    
        #     if pS > 0.05: 
        #         Sstationary = 'No'

        #     elif pS <= 0.05:
        #         Sstationary = 'Yes'
                
        #     Ns.pop(0) # ensuring that the adfuller test is conducted on N across 1K time steps
        #     Ss.pop(0) # ensuring that the adfuller test is conducted on S across 1K time steps
            
            
        #     if Nstationary == 'Yes' and Sstationary == 'Yes':
        #         # Sweet. Stationarity achieved. 
        #         print(iteration, ':  N=', N, '| S=', S, ' |  stationarity in N:', 
        #             Nstationary, ' |  stationarity in S:', Sstationary)
                
        #         print('Stationarity achieved in N and S across 1000 time steps.')
        #         end_time = time.time()
        #         print("Run time = {:.3f} seconds".format(end_time - start_time))
        #         break
                
            
    # Prompt user for parameters
    # sims=int(input("How many simulations do you want to run?\n"))
    # size=input("Do you want to run large- or small-scale simulations? Please enter: 'large' or 'small'\n").lower().strip()

sims=10
OUT = open(GenPath + "simData.csv", "w+")
h = headings()
print(h, file=OUT)
iterable=[]
for i in range(0, sims): # i=sim
    # Generate parameters for each simulation
    seed= i
    # lgp=0.999
    # imr=1
    # rpr=0.007
    # length=1000
    # dtr=0.007
    lgp=np.random.RandomState(seed).uniform(low=0.99, high=1.0, size = None)
    imr = np.random.RandomState(seed).randint(low=1,high=5)
    rpr=np.random.RandomState(seed).uniform(low=0.01, high=0.1, size = None) 
    # dtr=np.random.RandomState(seed).uniform(low=0.01, high=0.1, size=None)
    length = int(np.random.RandomState(seed).uniform(low=1, high=1000))
    iterable.append([i, lgp, imr, rpr, length])

pool= Pool(processes=10)
pool.starmap(run_model, iterable)

# ### Memory testing block
# OUT = open(GenPath + "simData.csv", "w+")
# h = headings()
# print(h, file=OUT)
# lgp=0.999
# sims = 0
# imr=1
# rpr=0.007
# length= 800
# dtr=0.007
# run_model(sims,lgp, imr, rpr, length)