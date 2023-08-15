# -*- coding: utf-8 -*-
from os.path import expanduser
import sys
import numpy as np
from multiprocessing import Pool

mydir = expanduser(".../IBM/")
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
    

def run_model(sim, size, seed=None):
    
    '''
    lgp: the log-series shape parameter; determines the structure of your regional pool
    imr: immigration rate. The number of individuals flowing in per time step. 
    rpr: General reproductive rate. The fraction of the community that reproduces each time step.
    dtr: Random chance of death. The fraction of the community that dies each time step. 
    
    '''
    lgp=np.random.RandomState(seed).uniform(low=0.99, high=1.0, size = None)
    imr = np.random.RandomState(seed).randint(1,10)
    rpr=np.random.RandomState(seed).uniform(low=0.01, high=0.1, size = None) 
    dtr=np.random.RandomState(seed).uniform(low=0.01, high=0.1, size=None)
    
    if size=='small':
        length = np.random.RandomState(seed).uniform(low=10, high=20)
        
    else: 
        length = np.random.RandomState(seed).uniform(low=1, high=100)

    iteration = 0 #counter for each simulation
    
    while iteration < 10*length:
        if iteration == 0:
            RAD[sim] = [],[] # add to rad file
        
        ''' Flow (directional passive dispersal) '''
        
        xcoords = RAD[sim,0]
        xcoords[:] = [i+1 for i in xcoords]
        
        
        ''' Immigration '''
        
        # Immigration from a log-series distributed regional pool. The randomly drawn numbers are 
        # the species IDs. On average, there will be more 1's than 2's, more 2's than 3's, etc.
        immigrants = np.random.RandomState(seed).logseries(lgp, size=imr)  
        
        # Add the immigrants to the local community
        community=RAD[sim,1]
        community.extend(immigrants)
        
        # Make the immigrants enter at the inlet
        xcoords.extend([0]*len(immigrants))

        ''' Emmigration''' 
        
        indices = [i for i, v in enumerate(xcoords) if v > length] # find every individual that flowed out of bounds
        
        for i in sorted(indices, reverse=True):
            del community[i]
            del xcoords[i]
        
        ''' Reproduction '''
        
        n = int(round(rpr * len(community))) # number of bouncing baby Bacilli 
        N = int(len(community))

        indices = np.random.RandomState(seed).choice(list(range(N)), size=n, replace=False) # choosing mothers at random
        bbb = [ community[i] for i in indices]
        bbb_x = [ xcoords[i] for i in indices]
        community.extend(bbb)
        xcoords.extend(bbb_x)

        ''' Death '''

        n = int(round(dtr * len(community))) # number of members in the community that will die

        indices = np.random.RandomState(seed).choice(list(range(N)), size=n, replace=False) # choosing deaths at random
        
        for i in sorted(indices, reverse=True):
            del community[i]
            del xcoords[i]


        iteration += 1
        RAD.flush()

         
        
        if iteration%10 == 0:
            print('Instance:', sim, ' iteration:', iteration, '|   N:', len(community), 'Length:', length)
            
            # Code that differs by the 'big' and 'small' parameter so I can artificially extend the community if I want it to be bigger
            if size == 'small':
                rad =[] 
                sp_ls = list(set(community))
                for i in sp_ls:
                    rad.append(community.count(i))
                rad = sorted(rad)
                div, ev = simpson_div(rad)
                
                # Per sim, record abundance, species richness, simpson, evenness
                OUT = open(GenPath + 'simData.csv', 'a')
                outlist = [sim, iteration, length, len(community), len(list(set(community))), div, ev]
                outlist = str(outlist).strip('[]')
                outlist = outlist.replace(" ", "")
                print(outlist, file=OUT)
                
            else: 
                big_com=community.copy()
                for i in community:
                    big_com.extend([i]*100) #for each community member, add 100 more as representatives
                rad =[] 
                sp_ls = list(set(big_com))
                for i in sp_ls:
                    rad.append(big_com.count(i))
                rad = sorted(rad)
                div, ev = simpson_div(rad)
                
                # Per sim, record abundance, species richness, simpson, evenness
                OUT = open(GenPath + 'simData.csv', 'a')
                outlist = [sim, iteration, length*100, len(big_com), len(list(set(big_com))), div, ev]
                outlist = str(outlist).strip('[]')
                outlist = outlist.replace(" ", "")
                print(outlist, file=OUT)

    
if __name__ == '__main__':
    # Prompt user for parameters
    sims=int(input("How many simulations do you want to run?\n"))
    size=input("Do you want to run large- or small-scale simulations? Please enter: 'large' or 'small'\n").lower().strip()
    
    OUT = open(GenPath + "simData.csv", "w+")
    h = headings()
    print(h, file=OUT)
    
    RAD = np.memmap(GenPath +"RAD.mymemmap", mode = 'w+', shape = (sims,2), dtype=object) # Creating an array with n=sims rows, two cols
    pool = Pool(processes=20)
    iterable=[]
    for i in range(0,sims):
        iterable.append([i, size])
    pool.starmap(run_model, iterable)
