# -*- coding: utf-8 -*-
from os.path import expanduser
import sys
import numpy as np
from multiprocessing import Pool
import numba as nb
from numba import jit
from numba.typed import List


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

@jit
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
    return D#, D/len(rad) # this latter value is simpson's evenness

@jit(nopython=True)
def run_model(sim, scale, all_sims):
    
    '''
    lgp: the log-series shape parameter; determines the structure of your regional pool
    imr: immigration rate. The number of individuals flowing in per time step. 
    rpr: General reproductive rate. The fraction of the community that reproduces each time step.
    dtr: Random chance of death. The fraction of the community that dies each time step. 
    
    '''
    lgp=np.random.uniform((0.99),(1.0))
    imr=np.random.randint((1),(3))
    rpr=np.random.uniform((0.01),(0.1)) 
    dtr=np.random.uniform((0.01), (0.1))

    if scale=='small':
        length = int(np.random.uniform((10), (20)))    
    else: 
        length = int(np.random.uniform((1), (1000)))

    iteration = 0 #counter for each simulation
    stop=10*length

    while iteration <= stop:
        if iteration == 0:
            RAD=np.empty(shape=(2,2), dtype='int16') # add to rad file
            # xcoords=np.empty(shape=1, dtype='int16')
            # community=np.empty(shape=1, dtype='int16')
        
        ''' Flow (directional passive dispersal) '''
        xcoords=RAD[:,0]
        xcoords=xcoords+1
        # xcoords[:] = [i+1 for i in xcoords] #not compatible with numba
        # optimize to check if they're going to be dead anyway

        ''' Immigration '''
        
        # Immigration from a log-series distributed regional pool. The randomly drawn numbers are 
        # the species IDs. On average, there will be more 1's than 2's, more 2's than 3's, etc.
        immigrants = np.random.logseries((lgp),(imr))
        # for i in range(len(immigrants)):
            # print(immigrants[i])
        # Add the immigrants to the local community

        community=RAD[:,1]

        # print(immigrants, community)
        community = np.append(community, immigrants)
        # for i in range(len(immigrants)):
            # community  = np.append(immigrants[i])
        # community.extend(immigrants) #not compatible with numba
        
        # Make the immigrants enter at the inlet
        xcoords=np.append(xcoords, [0]*len(immigrants))
        # xcoords.extend([0]*len(immigrants))

        ''' Emmigration''' 
        
        indices = [i for i, v in enumerate(xcoords) if v > length] # find every individual that flowed out of bounds
        
        for i in sorted(indices, reverse=True):
           community=np.delete(community, i)
           xcoords=np.delete(xcoords, i)
        
        ''' Reproduction '''
        
        n = int(round(rpr * len(community))) # number of bouncing baby Bacilli 
        N = int(len(community))

        indices = np.random.choice(np.arange(N), size=n, replace=False) # choosing mothers at random
        bbb = [ community[i] for i in indices]
        bbb_x = [ xcoords[i] for i in indices]
        # community.extend(bbb)
        xcoords=np.append(xcoords, bbb_x)
        community=np.append(community, bbb)
        # xcoords.extend(bbb_x)

        ''' Death '''

        n = int(round(dtr * len(community))) # number of members in the community that will die

        indices = np.random.choice(np.arange(N), size=n, replace=False) # choosing deaths at random
        
        for i in sorted(indices, reverse=True):
            community=np.delete(community, i)
            xcoords=np.delete(xcoords, i)


        iteration += 1
        RAD=np.vstack((xcoords, community)).T
        
        if iteration%10 == 0:
            print('Simulation:', sim, '|   N:', len(community), 'Length:', length)

        # Code that extends large community artifically unless run on small/special scale
        # if (scale=='small'): #only record this information if it is the last simulation
        #     rad =[] 
        #     sp_ls = list(set(community))
        #     for i in sp_ls:
        #         rad.append(community.count(i))
        #     rad = sorted(rad)
        #     div= simpson_div(rad)
            
        #     # Per sim, record abundance, species richness, simpson
        #     outlist = sim,length,div
        #     # outlist = str(outlist).strip('[]')
        #     # outlist = outlist.replace(" ", "")
        #     OUT[sim]=outlist
            
        # elif (scale=='large'): 
        #     big_com=community.copy()
        #     for i in community:
        #         big_com.extend([i]*100) #for each community member, add 100 more as representatives
        #     rad =[] 
        #     sp_ls = list(set(big_com))
        #     for i in sp_ls:
        #         rad.append(big_com.count(i))
        #     rad = sorted(rad)
        #     div = simpson_div(rad)
                
        #     # Per sim, record abundance, species richness, simpson, evenness
        #     outlist = sim,length*100,div
        #     # outlist = str(outlist).strip('[]')
        #     # outlist = outlist.replace(" ", "")
        #     OUT[sim]=outlist

    
# Prompt user for parameters
# sims=int(input("How many simulations do you want to run?\n"))
sims=1
# scale=input("Do you want to run large- or small-scale simulations? Please enter: 'large' or 'small'\n").lower().strip()
scale='large'
# Initialize the files simData (holds simpson diversity data) and RAD (holds species IDs and coordinates)
OUT = open(GenPath + "simData.csv", "w+")
h = headings()
print(h, file=OUT)
# RAD = np.memmap(GenPath +"RAD.mymemmap", mode = 'w+', shape = (sims,2), dtype=('int64')) # Creating an array with n=sims rows, two cols

# Run the simulations
pool = Pool(processes=1)
iterable=[]

# 1/3 of the simulations will be run on the normal scale to get the smallest two orders of magnitude
for i in range(0,sims):
    iterable.append([i, scale, sims])


pool.starmap(run_model, iterable)
