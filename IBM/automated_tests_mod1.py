import os 
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import statsmodels.tsa.stattools as sta
import random
import time
from multiprocessing import Pool
from os.path import expanduser

mydir = expanduser("/proj/gibbons/kramos/github/IBT-and-the-Gut-Microbiome/IBM/")
GenPath = mydir + 'data/small_tests_50-100/'

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

def run_model(sim, alpha, im, r, length, path):
    ################# DECLARE SOME STUFF #########################################
    COMM = np.empty(shape=1, dtype='uint16')   # community list
    COORDS = np.empty(shape=1, dtype='uint16') # coordinates ... indices correspond to COMM
    Ns = np.empty(shape=1, dtype='uint16')     # list to track total abundance
    Ss = np.empty([], dtype='uint16')     # list to track taxa richness

    Nstationary = 'No' # will indicate whether stationarity in N has been reached
    Sstationary = 'No' # will indicate whether stationarity in S has been reached


    ################# SIMULATE SOME STUFF #########################################

    start_time = time.time()
    time_step = 0
    while time_step == time_step: # WTH? An infinite loop? Is this guy nuts?
        

        ################# SIMULATE SOME PROCESSES #########################################
        
        processes = ['immigration', 'reproduction', 'flow']
        seed = None
        random.seed(seed) #initialize a new seed per time step
        # The order in which immigration, reproduction, and flow/emigration occur
        # is randomized in each time step ... to prevent systemic artifacts
        random.shuffle(processes)
        
        for p in processes:
            if p == 'immigration':
                # immigration from a log-series metacommunity ... as in Hubbell 2001
                immigrants = np.random.RandomState(seed).logseries(alpha, size = im).astype('uint16')
                COMM = np.concatenate((COMM, immigrants))
                # individuals enter at the upstream edge
                COORDS = np.concatenate((COORDS, np.array([0]* im, dtype = 'uint16'))) 
                
            elif p == 'reproduction':
                n = int(round(r * len(COMM))) # number of bouncing baby Bacilli 
                N = int(len(COMM))

                if n > 0:
                    indices = random.sample(list(range(N)), n) # choosing mothers at random
                    bbb = np.array(COMM)[indices] # Filter out non-selected indices
                    bbb_x = np.array(COORDS)[indices] # Filter out non-selected indices
                    COMM = np.concatenate((COMM, bbb))
                    COORDS=np.concatenate((COORDS, bbb_x))
                
            elif p == 'flow':
                COORDS = COORDS + 1
                indices = COORDS <= length
                COMM = COMM[indices]
                COORDS = COORDS[indices]
        
        
        ################# CALCULATE SOME STUFF #########################################
        
        N = len(COMM)
        Ns = np.append(Ns, N).astype('uint16')
        S = len(list(set(COMM)))
        Ss = np.append(Ss, S).astype('uint16')
        
        
        ################# CHECK FOR STATIONARITY #########################################
        
        # Once we've burned through 1000 time steps, check for stationarity in N and S
        # stationarity is like a fluctuating equilibrium
        
        if len(Ns) > 1000:
            '''
            Augmented Dickey-Fuller test. Tests for stationarity (i.e., serial autocorrelation).
            https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.adfuller.html
            
            We'll test for stationarity with respect to N and S.
            '''
            
            AugmentedDickeyFuller = sta.adfuller(Ns)
            val, pN = AugmentedDickeyFuller[0:2]
                
            AugmentedDickeyFuller = sta.adfuller(Ss)
            val, pS = AugmentedDickeyFuller[0:2]

            if pN > 0.05: 
                Nstationary = 'No'

            elif pN <= 0.05:
                Nstationary = 'Yes'
                    
            if pS > 0.05: 
                Sstationary = 'No'

            elif pS <= 0.05:
                Sstationary = 'Yes'
                
            Ns = np.delete(Ns, 0) #Ns.pop(0) # ensuring that the adfuller test is conducted on N across 1K time steps
            Ss = np.delete(Ss, 0) #Ss.pop(0) # ensuring that the adfuller test is conducted on S across 1K time steps
            
            if Nstationary == 'Yes' and Sstationary == 'Yes': # Sweet. Stationarity achieved. 
                unique, counts = np.unique(COMM, return_counts=True)
                sp_d = dict(zip(unique, counts))
                rad=list(sp_d.values())
                rad = sorted(rad)
                div, ev = simpson_div(rad)
                
                # Per sim, record abundance, species richness, simpson, evenness
                OUT = open(path, 'a+')
                outlist = [sim, time_step, length, N, S, div, ev]
                outlist = str(outlist).strip('[]')
                outlist = outlist.replace(" ", "")

                print(time_step, ':  N=', N, '| S=', S, ' |  stationarity in N:', 
                    Nstationary, ' |  stationarity in S:', Sstationary)
                print(outlist, file=OUT)
                print('Stationarity achieved in N and S across 1000 time steps.')
                end_time = time.time()
                print("Run time = {:.3f} seconds".format(end_time - start_time))

                break # Kill the infinite loop
                
                
        ################# PRINT SOME STUFF AFTER EACH TIME STEP ##################################
        
        print(time_step, ':  N=', N, '| S=', S, ' |  stationarity in N:', 
            Nstationary, ' |  stationarity in S:', Sstationary)
        
        time_step += 1

################# Code to run the model ######################

im_rates = [1, 25, 50]
im_rate_desc = ['low_im', 'med_im', 'high_im']
r_rates = [0.001, 0.01, 0.1]
r_rate_desc = ['low_r', 'med_r', 'high_r']

for i in range(0,3): # for each possible immigration rate 
    seed = None
    alpha = 0.9999
    iterable = []
    im = im_rates[i]
    im_desc = im_rate_desc[i]
    
    for j in range(0,3): # for each possible reproduction rate
        r = r_rates[j]
        r_desc = r_rate_desc[j]
    
        sims=300
        path = GenPath + im_desc + '_' + r_desc + ".csv"
        OUT = open(path, "w")
        h = headings()
        print(h, file = OUT)
        
        for i in range(sims): # generate iterable for each sim with the parameters selected
            length = np.random.RandomState(seed).uniform(low=np.log(50), high=np.log(100))
            iterable.append([i, alpha, im, r, np.exp(length), path])
        
        pool= Pool(processes=10)
        pool.starmap(run_model, iterable) 




