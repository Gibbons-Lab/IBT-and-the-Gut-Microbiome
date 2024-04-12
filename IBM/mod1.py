import numpy as np
import statsmodels.tsa.stattools as sta
import random
import time

################# DECLARE SOME STUFF #########################################

COMM = []   # community list
COORDS = [] # coordinates ... indices correspond to COMM
Ns = []     # list to track total abundance
Ss = []     # list to track taxa richness

alpha = 0.999 # log-series parameter. 
              # The closer it is to 1.0, the more species you'll get from a sample
              
flow_rate = 1 # units of distance each individual flows downstream per time step
length = 1000 
im = 1      # immigration rate = individuals inflowing per time step
r = 0.007   # per capita probability of reproduction in a given time step
            # Note: small changes in r produce big changes in N

Nstationary = 'No' # will indicate whether stationarity in N has been reached
Sstationary = 'No' # will indicate whether stationarity in S has been reached


################# SIMULATE SOME STUFF #########################################

start_time = time.time()
time_step = 0
while time_step == time_step: # WTH? An infinite loop? Is this guy nuts?
	

    ################# SIMULATE SOME PROCESSES #########################################
    
    processes = ['immigration', 'reproduction', 'flow']
    # The order in which immigration, reproduction, and flow/emigration occur
    # is randomized in each time step ... to prevent systemic artifacts
    random.shuffle(processes)
    
    for p in processes:
        if p == 'immigration':
            # immigration from a log-series metacommunity ... as in Hubbell 2001
            COMM.extend(np.random.logseries(alpha, size=im).tolist()) 
            # individuals enter at the upstream edge
            COORDS.extend([0] * im) 
            
        elif p == 'reproduction':
            # create a binary (0/1) list for which individual will reproduce
            # 1's and 0's are chosen at random via the per capita probability of
            # reproduction in a given time step
            ls = np.random.binomial(1, r, size=len(COMM)).tolist()
            
            # find the indices with 1's (i.e., individuals that will reproduce)
            indices = [i for i, value in enumerate(ls) if value == 1]
            
            # reproduce the individuals with 1's
            daughters = [COMM[i] for i in indices]
            # daughters get the same coordinates as their mothers
            daughter_coords = [COORDS[i] for i in indices]
            
            # add the daughters to the community
            COMM.extend(daughters)
            COORDS.extend(daughter_coords)
            
        elif p == 'flow':
            indices = COORDS <= length
            COMM = COMM[indices]
            COORDS = COORDS[indices]
    
    ################# CALCULATE SOME STUFF #########################################
    
    N = len(COMM)
    Ns.append(N)
    S = len(list(set(COMM)))
    Ss.append(S)
    
    
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
            
        Ns.pop(0) # ensuring that the adfuller test is conducted on N across 1K time steps
        Ss.pop(0) # ensuring that the adfuller test is conducted on S across 1K time steps
        
        
        if Nstationary == 'Yes' and Sstationary == 'Yes':
            # Sweet. Stationarity achieved. 
            
            print(time_step, ':  N=', N, '| S=', S, ' |  stationarity in N:', 
                  Nstationary, ' |  stationarity in S:', Sstationary)
            
            print('Stationarity achieved in N and S across 1000 time steps.')
            end_time = time.time()
            print("Run time = {:.3f} seconds".format(end_time - start_time))
            
            break # Kill the infinite loop
            
    ################# PRINT SOME STUFF AFTER EACH TIME STEP ##################################
    
    print(time_step, ':  N=', N, '| S=', S, ' |  stationarity in N:', 
          Nstationary, ' |  stationarity in S:', Sstationary)
    
    time_step += 1