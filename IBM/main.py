import os 
import numpy as np
import statsmodels.tsa.stattools as sta
import random
import time
from multiprocessing import Pool
from os.path import expanduser

mydir = expanduser("/proj/gibbons/kramos/github/IBT-and-the-Gut-Microbiome/IBM/")
GenPath = mydir + 'data/'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

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

def run_model(sim, lgp, im, r, length, flow_rate):
    ################# DECLARE VARIABLES/LISTS #########################################

    COMM = np.empty(shape=1, dtype='uint16')   # community list
    COORDS = np.empty(shape=1, dtype='uint16') # coordinates ... indices correspond to COMM
    Ns = np.empty(shape=1, dtype='uint16')     # list to track total abundance
    Ss = np.empty([], dtype='uint16')     # list to track taxa richness

    Nstationary = 'No' # will indicate whether stationarity in N has been reached
    Sstationary = 'No' # will indicate whether stationarity in S has been reached


    ################# SIMULATION START #########################################

    start_time = time.time()
    time_step = 0
    while time_step == time_step: 
        

        ################# PROCESSES #########################################
        
        processes = ['immigration', 'reproduction', 'flow']
        seed = None
        random.seed(seed) #initialize a new seed per time step
        random.shuffle(processes)
        # The order in which immigration, reproduction, and flow/emigration occur
        # is randomized in each time step ... to prevent systemic artifacts
        
        for p in processes:
            if p == 'immigration':
                immigrants = np.random.RandomState(seed).logseries(lgp, size = im).astype('uint16') # immigration from a log-series metacommunity ... as in Hubbell 2001
                COMM = np.concatenate((COMM, immigrants))
                COORDS = np.concatenate((COORDS, np.array([0]* im, dtype = 'uint16'))) # individuals enter at the upstream edge
                
            elif p == 'reproduction':
                n = int(round(r * len(COMM))) # number of bouncing baby Bacilli 
                N = int(len(COMM))

                if n > 0:
                    indices = random.sample(list(range(N)), n) # choosing mothers at random
                    bbb = np.array(COMM)[indices] # Filter out non-selected indices
                    bbb_x = np.array(COORDS)[indices] # Filter out non-selected indices
                    COMM = np.concatenate((COMM, bbb))
                    COORDS=np.concatenate((COORDS, bbb_x))
                
            elif p == 'flow': # add one unit to each individual's current x-coord
                COORDS = COORDS + flow_rate
                indices = COORDS <= length # remove individuals who have flown out of the system
                COMM = COMM[indices]
                COORDS = COORDS[indices]
        
        
        ################# CALCULATIONS #########################################
        
        N = len(COMM)
        Ns = np.append(Ns, N).astype('uint16')
        S = len(list(set(COMM)))
        Ss = np.append(Ss, S).astype('uint16')
        
        
        ################# CHECK FOR STATIONARITY #########################################
        
        # Once we've burned through 1000 time steps, check for stationarity in N and S
        
        if len(Ns) > 1000:
            '''
            Augmented Dickey-Fuller test. Tests for stationarity (i.e., serial autocorrelation).
            https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.adfuller.html
            
            We'll test for stationarity with respect to N and S.
            '''
            
            AugmentedDickeyFuller = sta.adfuller(Ns)
            val, pN = AugmentedDickeyFuller[0:2]
            nobs = AugmentedDickeyFuller[3]

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
                
            Ns = np.delete(Ns, 0) # ensuring that the adfuller test is conducted on N across 1K time steps
            Ss = np.delete(Ss, 0) # ensuring that the adfuller test is conducted on S across 1K time steps
            
            if Nstationary == 'Yes' and Sstationary == 'Yes': # Stationarity achieved. 
                unique, counts = np.unique(COMM, return_counts=True)
                sp_d = dict(zip(unique, counts))
                rad=list(sp_d.values())
                rad = sorted(rad)
                div, ev = simpson_div(rad)
                
                # Per sim, record abundance, species richness, simpson, evenness
                OUT = open(GenPath + file_name, 'a+')
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


sims=10
file_name = "test.csv"
OUT = open(GenPath + file_name, "w")
h = headings()
print(h, file = OUT)

iterable = []

for i in range(0, sims): # i=sim
    seed = None
    # Generate parameters for each simulation
    lgp = 0.9999 # log-series parameter. # The closer it is to 1.0, the more species you'll get from a sample 
    low_im = 1 # immigration rate = individuals inflowing per time step
    med_im = 25
    high_im = 50
    low_r = 0.001 # per capita probability of reproduction in a given time step 
    med_r = 0.01
    high_r = 0.1
    flow_rate = 1 # units of distance each individual flows downstream per time step
    # Decide which parameters to move forward with
    im=med_im
    r=med_r
    length = np.random.RandomState(seed).uniform(low=np.log(1), high=np.log(1000)) # This parameter determines the scale of your simulations
    iterable.append([i, lgp, im, r, np.exp(length), flow_rate])

pool= Pool(processes=10)
pool.starmap(run_model, iterable)



