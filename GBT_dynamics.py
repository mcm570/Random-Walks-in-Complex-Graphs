import numpy as np
import random
import matplotlib.pyplot as plt
import myconfig

def tau_tree_analytic(depth):
    return 8*(2**(depth)-1)*(1-1/(2**depth)) # analytic formula for the hitting time of a binary tree

def simulation_tree():
    depth_arr = np.arange(1,11)
    t_arr = np.empty_like(depth_arr)
    for depth in depth_arr:
        t = 60000 #number of steps that will be considered

        n_trials = 5000 #number of trials we will do
        t_trial = np.array([])

        for trial in range(0,n_trials): 
            k = 0
            T = 0 #time variable

            probability = 2/3

            while (T< t) and (k < 2*depth):
                if k==0: p=1
                elif k<depth: p = probability
                elif k==depth: p=1/2
                else: p=1-probability
                if (random.random() < p):
                    k = k + 1
                else:
                    k = k - 1 if k > 0 else 0
                T = T + 1

            if (k == 2*depth): t_trial = np.append(t_trial, T)

        print(len(t_trial)/n_trials)
        print('depth = ',depth, 't = ',t,', mean time = ',np.mean(t_trial), ' analytic = ', tau_tree_analytic(depth), 'standard deviation = ', np.std(t_trial))

        t_arr[depth-1] = np.mean(t_trial)

    #plot
    plt.figure(1)
    plt.scatter(depth_arr, t_arr, label="simulation")
    # plt.scatter(depth_arr, tau_tree_analytic(depth_arr), label="analytic")
    plt.xlabel("depth")
    plt.ylabel("mean time")
    # plt.savefig("hitting_time_tree.png", bbox_inches='tight', dpi=300)
    
    plt.figure(2)
    plt.yscale('log')
    plt.scatter(depth_arr, t_arr, label="simulation")
    # plt.scatter(depth_arr, tau_tree_analytic(depth_arr), label="analytic")
    plt.xlabel("depth")
    plt.ylabel("mean time")
    # plt.savefig("hitting_time_tree_log.png", bbox_inches='tight', dpi=300)
    # plt.legend()
    plt.show()

simulation_tree()
