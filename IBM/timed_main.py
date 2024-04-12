# -*- coding: utf-8 -*-
from os.path import expanduser
import sys
import numpy as np
from multiprocessing import Pool
import time

mydir = expanduser("/proj/gibbons/kramos/github/IBT-and-the-Gut-Microbiome/IBM/")
GenPath = mydir + 'data/'

def headings():
    headings = 'sim,'
    headings += 'sim_clone,'
    headings += 'ct,'
    headings += 'Length,'
    headings += 'TotalAbundance,'
    headings += 'SpeciesRichness,'
    headings += 'Simpson,'
    headings += 'Simpson_e'

    return headings

# def del_list_numpy(l, id_to_del):
#     arr = np.array(l, dtype='int32')
#     return list(np.delete(arr, id_to_del))

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
    

def run_model(sim, sim_clone, size, lgp, imr, rpr, dtr):
    
    iteration = 0 #counter for each simulation
    stopping_point=10*length
    strikes = 0
    seed= sim_clone

    while (iteration <= stopping_point) and (strikes <=10): # end the simulation if the community has been maximized for more than 10 cycles
        if iteration == 0:
            N=0
            # initialize RAD file for this sim
            # RAD will be a 3D array and will look like the following
            # RAD = [ [[x1, x2, ..., xn],[s1, s2, ..., sn]], ... ]

            RAD=np.zeros(shape=(10,2,0), dtype='int16')
            # initialize x coordinates and community
            xcoords=[]
            community=[]
            # RAD=np.empty(shape=(10,2))
            
            # RAD = np.memmap(GenPath +"RAD.mymemmap", mode = 'w+', shape = (10,2), dtype=object) # Creating an array with n=sims rows, two cols
        N=int(len(community))
        S=len(set(np.unique(community)))
        start_time=time.time()
        
        ''' Flow (directional passive dispersal) '''
        # xcoords = RAD[sim,0]      
        xcoords[:] = [i+1 for i in xcoords]
        # xcoords= xcoords + 1
        end_time=time.time()
        print("Flow={:.3f} seconds".format(end_time-start_time))
        
        ''' Immigration '''
        start_time=time.time()
        # Immigration from a log-series distributed regional pool. The randomly drawn numbers are 
        # the species IDs. On average, there will be more 1's than 2's, more 2's than 3's, etc.
        if N<=10000: # Only present new immigrations if community is not at its max
            immigrants = np.array(np.random.RandomState(seed).logseries(lgp, size=imr))
            
            # Add the immigrants to the local community
            # community=RAD[sim,1]
            community=np.concatenate((community,immigrants))
            
            # Make the immigrants enter at the inlet
            # xcoords=np.append(xcoords, np.array([0]*len(immigrants)))
            xcoords=np.concatenate((xcoords, np.array([0]*len(immigrants))))
            end_time=time.time()
            print('Immigration={:.3f} seconds'.format(end_time-start_time))
        else: 
            strikes += 1

        ''' Emmigration''' 
        start_time=time.time()
        indices = [i for i, v in enumerate(xcoords) if v > length] # find every individual that flowed out of bounds
        
        community = np.delete(community, indices)
        xcoords = np.delete(xcoords, indices)

        # for i in sorted(indices, reverse=True):
        #     del community[i]
        #     del xcoords[i]

        end_time=time.time()
        print('Emmigration={:.3f} seconds'.format(end_time-start_time))

        ''' Reproduction '''
        start_time=time.time()
        if N<=100000: # Only present new immigrations if community is not at its max
            n = int(round(rpr * len(community))) # number of bouncing baby Bacilli 
            N = int(len(community))

            if n > 0:
                indices = np.random.RandomState(seed).choice(list(range(N)), size=n, replace=False) # choosing mothers at random
                # bbb = np.array(community[i] for i in indices)
                # bbb_x = np.array(xcoords[i] for i in indices)
                # Filter out non-selected indices
                bbb = np.array(community)[indices]
                bbb_x = np.array(xcoords)[indices]
                community=np.concatenate((community, bbb))
                xcoords=np.concatenate((xcoords, bbb_x))
                end_time=time.time()
                print('Reproduction={:.3f} seconds'.format(end_time-start_time))
        else: 
            strikes += 1

        ''' Death '''
        start_time=time.time()
        n = int(round(dtr * len(community))) # number of members in the community that will die

        indices = np.random.RandomState(seed).choice(list(range(N)), size=n, replace=False) # choosing deaths at random
        
        community = np.delete(community, indices)
        xcoords = np.delete(xcoords, indices)

        # for i in sorted(indices, reverse=True):
        #     del community[i]
        #     del xcoords[i]
        
        end_time=time.time()
        print('Death={:.3f} seconds'.format(end_time-start_time))
        iteration += 1

        start_time=time.time()
        # RAD[sim_clone]=community
        # RAD[sim_clone]=np.array(list(zip(xcoords, community)))
        end_time=time.time()
        print('Update RAD {:.3f} seconds'.format(end_time-start_time))
         
        if iteration%50 == 0:
            start_time=time.time()
            # Code that differs by the 'big' and 'small' parameter so I can artificially extend the community if I want it to be bigger
            if size == 'small':
                rad =[] 
                sp_ls = list((np.unique(community)))
                for i in sp_ls:
                    rad.append(community.count(i))
                rad = sorted(rad)
                div, ev = simpson_div(rad)
                
                # Per sim, record abundance, species richness, simpson, evenness
                OUT = open(GenPath + 'simData.csv', 'a')
                outlist = [sim, iteration, length, len(community), len(list(np.unique(community))), div, ev]
                outlist = str(outlist).strip('[]')
                outlist = outlist.replace(" ", "")
                print(outlist, file=OUT)
                
            else: 
                big_com=community.copy()
                for i in community:
                    big_com=np.append(big_com, [i]*100) #for each community member, add 100 more as representatives
                if iteration==stopping_point:
                    rad =[] 
                    sp_ls = list((np.unique(big_com)))
                    for i in sp_ls:
                        rad.append(big_com.count(i))
                    rad = sorted(rad)
                    div, ev = simpson_div(rad)
                else: 
                    div,ev=float('nan'), float('nan')
                
                # Per sim, record abundance, species richness, simpson, evenness
                OUT = open(GenPath + 'simData.csv', 'a')
                outlist = [sim, sim_clone, iteration, length*100, len(big_com), len(list(np.unique(big_com))), div, ev]
                outlist = str(outlist).strip('[]')
                outlist = outlist.replace(" ", "")
                print(outlist, file=OUT)
                end_time=time.time()
                print('Div metrics={:.3f} seconds'.format(end_time-start_time))
            
            print('Instance:', sim, ' iteration:', iteration, '|   N:', len(community), 'Length:', length)
            
# Prompt user for parameters
# sims=int(input("How many simulations do you want to run?\n"))
# size=input("Do you want to run large- or small-scale simulations? Please enter: 'large' or 'small'\n").lower().strip()

sims=10
size=1
OUT = open(GenPath + "simData.csv", "w+")
h = headings()
print(h, file=OUT)

for i in range(0, sims): # i=sim
    # Generate parameters for each simulation
    seed= i
    lgp=np.random.RandomState(seed).uniform(low=0.99, high=1.0, size = None)
    imr = np.random.RandomState(seed).randint(low=1,high=10)
    rpr=np.random.RandomState(seed).uniform(low=0.01, high=0.1, size = None) 
    dtr=np.random.RandomState(seed).uniform(low=0.01, high=0.1, size=None)
    length = int(np.random.RandomState(seed).uniform(low=700, high=1000))
    print("Generated params")
    #Generate clone IDs for each simulation
    for j in range(0, 10): # j=sim_clone
        if j==0:
             iterable=[]
             print('Started iterable')
        iterable.append([i, j, size, lgp, imr, rpr, dtr])
        if j==9:
            print('starting pool')
            pool= Pool(processes=10)
            pool.starmap(run_model, iterable)
    print(iterable)


# # RAD = np.memmap(GenPath +"RAD.mymemmap", mode = 'w+', shape = (sims,2), dtype=object) # Creating an array with n=sims rows, two cols
# pool = Pool(processes=1)
# iterable=[]
# for i in range(0,sims):
#     iterable.append([i, size, RAD])
# pool.starmap(run_model, iterable)