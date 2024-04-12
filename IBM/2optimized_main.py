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
# If I'm going to use numba, file I/O will only work on python, but I can separate the computations

def headings():
    headings = 'sim,'
    headings += 'ct,'
    headings += 'Length,'
    headings += 'TotalAbundance,'
    headings += 'SpeciesRichness,'
    headings += 'Simpson,'

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
def generate_params(scale):
     
    #  '''
    # lgp: the log-series shape parameter; determines the structure of your regional pool
    # imr: immigration rate. The number of individuals flowing in per time step. 
    # rpr: General reproductive rate. The fraction of the community that reproduces each time step.
    # dtr: Random chance of death. The fraction of the community that dies each time step. 
    
    # '''
    lgp=np.random.uniform((0.99), (1.0))
    imr=np.random.randint((1),(10))
    rpr=np.random.uniform((0.01),(0.1)) 
    dtr=np.random.uniform((0.01), (0.1))

    #Generate the length of the simulation
    if scale=='small':
        length = int(np.random.uniform((10), (20)))    
    else: 
        length = int(np.random.uniform((1), (1000)))

    # Simulation stopping point
    stop=10*length
    # params=[lgp, imr, rpr, dtr, length, scale, stop]
    return(lgp, imr, rpr, dtr, length, scale, stop)

@jit(nopython=True)
def procs(xcoords, community, params):
    xcoords=xcoords
    community=community
    outlist=[]
    lgp, imr, rpr, dtr, length, scale, stop = [i for i in params]
    ''' Flow (directional passive dispersal) '''
    xcoords[:] = [i+1 for i in xcoords] #not compatible with numba
    # optimize to check if they're going to be dead anyway

    ''' Immigration '''
    # Immigration from a log-series distributed regional pool. The randomly drawn numbers are 
    # the species IDs. On average, there will be more 1's than 2's, more 2's than 3's, etc.
    immigrants = np.random.logseries((lgp),(imr))
    
    for i in range(len(immigrants)):
        community.append(immigrants[i])

    xcoords.extend([0]*len(immigrants)) # Make the immigrants enter at the inlet

    ''' Emmigration''' 
    indices = [i for i, v in enumerate(xcoords) if v > length] # find every individual that flowed out of bounds
    
    for i in sorted(indices, reverse=True):
        del community[i]
        del xcoords[i]
    
    ''' Reproduction '''
    
    n = int(round(rpr * len(community))) # number of bouncing baby Bacilli 
    N = int(len(community))

    indices = np.random.choice(np.arange(N), size=n, replace=False) # choosing mothers at random
    bbb = [ community[i] for i in indices]
    bbb_x = [ xcoords[i] for i in indices]
    community.extend(bbb)
    xcoords.extend(bbb_x)

    ''' Death '''

    n = int(round(dtr * len(community))) # number of members in the community that will die

    indices = np.random.choice(np.arange(N), size=n, replace=False) # choosing deaths at random
    
    for i in sorted(indices, reverse=True):
        del community[i]
        del xcoords[i]

    return(xcoords, community)

def iterator(params, sim): # Will decide how many times to run procs, is going to connect python and numba
    lgp, imr, rpr, dtr, length, scale, stop = [i for i in params]
    params=params
    xcoords=RAD[sim, 0]
    community=RAD[sim, 1]
    for i in range(0,stop): #Iterate procs until stopping point has been reached
        xcoords, community = procs(xcoords, community, params )
        RAD.flush() #update RAD file after each process
        if i%100==0: #generate outlist
            print('Instance:', sim, ' iteration:', iteration, '|   N:', len(community), 'Length:', length)
            rad=[]
            sp_ls=list(set(community))
            for i in sp_ls:
                rad.append(community.count(i))
            rad=sorted(rad)
            div=simpson_div(rad)

            outlist=[sim, i, length, len(community), len(list(set(community))), div]
            outlist = str(outlist).strip('[]')
            outlist = outlist.replace(" ", "")
            print(outlist, file=OUT)

def run_model(sim, scale):
    # Generate params
    sim = sim
    scale=scale
    # Generate parameters for new set of simulations
    params=[generate_params(scale)]
    # Based on parameters, run the simulation
    iterator(params, sim)

#######################
# Prompt user for parameters
# sims=int(input("How many simulations do you want to run?\n"))
# size=input("Do you want to run large- or small-scale simulations? Please enter: 'large' or 'small'\n").lower().strip()
    
sims=100
scale='large'
   
# Generate the files I will be saving to    

OUT = open(GenPath + "simData.csv", "w+")
h = headings()
print(h, file=OUT)

RAD = np.memmap(GenPath +"RAD.mymemmap", mode = 'w+', shape = (sims,2), dtype=('int64')) # Creating an array with n=sims rows, two cols

# Run the simulations
pool = Pool(processes=10)
iterable=[]
for i in range(0,sims):
    iterable.append([i, scale])
pool.starmap(run_model, iterable)

        
            

