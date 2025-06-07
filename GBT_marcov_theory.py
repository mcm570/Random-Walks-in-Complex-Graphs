import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.optimize import curve_fit
import math
import random 

def generate_random_probabilities(D):
    probabilities = np.random.rand(D)  # generate D random numbers between 0 and 1
    probabilities /= probabilities.sum()  # Normalize
    return probabilities


def P_simple_binary_tree(P,n,p_dir):
    #Define the values of the P matrix for the simple binary tree: 
    for i in range(0,n):
        for j in range(0,n):
            if j == 2*i+1: 
                P[i][j] = p_dir[0]
                P[j][i] = p_dir[0]
            elif j == 2*i+2:
                P[i][j] = p_dir[1]
                P[j][i] = p_dir[1]


def normalize_rows(matrix):
    row_sums = matrix.sum(axis=1)
    normalized_matrix = matrix/row_sums[:, np.newaxis]
    return normalized_matrix

# --- Función modelo: exponencial ---
def exponential(t, A, lamb):
    return A * np.exp(-lamb * t)


def p_0_graphics(p_0_ex,D,tau,steps):    #this will plot p_0 in function of time and D
    '''
        # --- Ajuste ---
        params, covariance = curve_fit(exponential, steps, p_0_ex, p0=(1e-4, 1e-4))
        A_fit, lambda_fit = params
        param_errors = np.sqrt(np.diag(covariance))

    '''
    # configuración del gráfico
    plt.xlabel('t [steps]', fontdict={'fontsize': 15})
    plt.ylabel(f'$p_0$', fontdict={'fontsize': 15})
    plt.title(f'Probability of arriving to the opposite node in function of the number of steps with D={D}', fontdict={'fontsize': 15, 'fontweight': 'bold'})
    plt.plot(p_0_ex, linewidth=0.5 ,markersize=0.5)
    #plt.plot(steps, exponential(steps, *params), 'r--', label=f'Ajuste:  A = {A_fit:.2e} ± {param_errors[0]:.1e}, λ = {lambda_fit:.2e} ± {param_errors[1]:.1e}')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # Añadir una línea vertical en t = T
    plt.axvline(x=tau, color='r', linestyle='--', linewidth=2, label=r'$\tau$')
    # Línea vertical en t=t_min
    t_min = 400
    plt.axvline(x=t_min, color='b', linestyle='--', linewidth=2, label=r'$t_{min}$')

    plt.legend(fontsize=15,loc='upper right')
    plt.show()
    plt.savefig(f"tree_D={D}.png", dpi=300, bbox_inches='tight')



D = 3 #depth of the simple binary tree
n1 = np.sum(2 ** np.arange(D + 1)) #number of nodes in simple binary tree
n2 = 2*np.sum(2 ** np.arange(D)) + 2**D #number of nodes in the glued binary tree
P = np.zeros((n2,n2)) #transition matrix
p_dir = np.ones(2)/2 #probability of moving to the right or to the left with equal probability
#p_dir = generate_random_probabilities(2) #probability of moving to the right or to the left with random probabilities
#p_dir = np.arange(1,3)

'''
Define de superior submatrix which is equal to the transition matrix of a simple binary tree 
'''
P_simple_binary_tree(P,n1,p_dir) 

'''
In function of this matrix, lets define the other values for the glued binary tree
'''
P[n1:n2 , n1:n2] = P[0:n2-n1 , 0:n2-n1] 

P[ 2**D-1 : n1 , n1:n2 ] = P[ 2**D-1 : n1 , 0:n2-n1 ] 
P[ n1:n2 , 2**D-1 : n1] = P[ 2**D-1 : n1 , n1:n2 ].T

''' 
Impose the condition of first arrival in the transition matriz
'''
P[n1,:] = 0
P[n1,n1] = 1

P = normalize_rows(P)


''' 
Lets estimate the hitting time
'''
t = 600 #number of steps that will be considered
v = np.zeros((t,n2)) #vector of probabilities of node i at time step t
v[0,0] = 1 #we impose the particle begins at the first node
#p_0 = np.zeros(t) #probability that the particle hits the final vertex for the 1st time at step t or less'
p_0_ex = np.zeros(t) #probability that the particle hits the final vertex for the 1st time at exactly t steps
steps = np.arange(0, t)

file = f"tau_vector_D={D}_t={t}.txt"
with open(file, 'w') as file:

    for T in range (1,t):
        v[T,:] = v[T-1,:]@P
        #p_0[T] = v[T,n-1] #probability of arriving to (0') in T or less steps
        p_0_ex[T] = v[T,n1] - v[T-1,n1]  #probability of arriving to (0') in exactly T steps
        #print ('t=',T,'p_0_ex[T]=', p_0_ex[T])
        if int(T%1000) == 0:
            print('T = ',T)
        file.write(f"{T} {T*p_0_ex[T]} \n")

        tau = np.dot(steps,p_0_ex)
    file.write(f"tau = {tau}\n")
    
        

tau_minus =  np.dot(steps[0:t-1],p_0_ex[0:t-1])
error = (tau - tau_minus)/tau
#print(P)

#print('p_0_ex = ', p_0_ex)
print('D = ',D,'t = ',t,'tau = ',tau, 'error = ', error)
p_0_graphics(p_0_ex,D,tau,steps)

'''
file = "tau_vector.txt"
with open(file, 'w') as file:
    file.write(f"{(steps*p_0_ex).T}")
'''