import numpy as np
import math
import random

for D in range (13,14):
    t = 60000 #number of steps that will be considered
    #D = 3 #dimension of the hypercube
    #n = 2**D #number of nodes

    '''we will do several trials in orden to estimate a mean of the time that the time the particle take to arrive to the opposite node.
    This mean time will be compared with the hitting time we estimated in previous program and should be similar'''
    n_trials = 1000 #number of trials we will do
    t_trial = np.zeros(n_trials) #the position i will save the time that the particle take to arrive to the opposite node in the trial i

    for trial in range(0,n_trials): 
        p_dir = np.zeros(D) #probability of moving along direction i
        dir = np.arange(D)   #number of direccions. We identify direccion x as dir[0], direction y as dir[1], etc.
        #s = generate_hypercube_nodes(D)



        x = np.zeros(D) #position vector. Initially is on node (0,..,0)

        #define de values of p_dir:
        #p_dir = generate_random_probabilities(D)
        p_dir = np.ones(D)/D

        T=0 #time variable

        while (T < t) and (not np.all(x == 1)):
            dir_t = np.random.choice(dir,p=p_dir) #direction of movement at time step T
            x[dir_t] = np.abs(x[dir_t] - 1)
            T = T+1
            #print( 'x[',T,']=',x )
        t_trial[trial] = T
        print('trial =', trial, 'T = ',T)

    #print(t_trial)
    print('D = ',D,', t = ',t,', mean time = ',np.mean(t_trial), 'standard deviation = ', np.std(t_trial))
