import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import math
import random 


def generate_hypercube_nodes(D):
    """
    Genera un array de dimensión D x 2^D, donde cada columna representa un nodo
    en un hipercubo de dimensión D.

    Parámetros:
    D (int): Dimensión del hipercubo

    Retorna:
    np.ndarray: Matriz de tamaño (D, 2^D) con cada columna como un nodo
    """
    n = 2 ** D  # Número total de nodos
    s = np.zeros((n, D), dtype=int)

    for i in range(n):
        binary_representation = np.binary_repr(i, width=D)
        s[i,:] = [int(bit) for bit in binary_representation]

    return s



def displacement(s1, s2):
    '''
    This function returns a boolean variable that show if s1 and s2 are consecutive nodes. If they are consecutive, the function also 
    returns the coordinate where the particle has moved.
    That is, if s1 = [0,0,2] and s2 = [0,0,3] then the fuction will return 2, which is the coordinated that have changed. 
    '''
    diff = s2 - s1
    nonzero = np.nonzero(diff)[0]  # índices donde hay diferencia

    if len(nonzero) == 1 and abs(diff[nonzero[0]]) == 1:
        return nonzero[0], True  # desplazamiento válido en un único eje
    else:
        return None, False  # no hay desplazamiento unitario válido
    

'''
lets create a function that will define randomly the vector p_dir[D] whose elements are 
the probabilites of moving along direction i. 
'''
def generate_random_probabilities(D):
    probabilities = np.random.rand(D)  # generate D random numbers between 0 and 1
    probabilities /= probabilities.sum()  # Normalize
    return probabilities


def p_0_graphics(p_0_ex,D,tau):    #this will plot p_0 in function of time and D
    # configuración del gráfico
    plt.xlabel('t [steps]', fontdict={'fontsize': 15})
    plt.ylabel(f'$p_0$', fontdict={'fontsize': 15})
    plt.title(f'Probability of arriving to the opposite node in function of the number of steps with D={D}', fontdict={'fontsize': 15, 'fontweight': 'bold'})
    plt.plot(p_0_ex, linewidth=0.5 ,markersize=0.5)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # Añadir una línea vertical en t = T
    plt.axvline(x=tau, color='r', linestyle='--', linewidth=2, label=r'$\tau$')
    # Línea vertical en t=t_min
    t_min = 30
    plt.axvline(x=t_min, color='b', linestyle='--', linewidth=2, label=r'$t_{min}$')
    plt.legend(fontsize=15,loc='upper right')
    plt.show()
    plt.savefig(f"hypercube_D={D}.png", dpi=300, bbox_inches='tight')


#for D in range(1,5):
D=3
#for t in range (10000, )
t = 100 #number of steps that will be considered
#D = 1 #dimension of the hypercube
n = 2**D #number of nodes
P = np.zeros((n,n)) #transition matrix
p_dir = np.zeros(D) #probability of moving along direction i
v = np.zeros((t,n)) #vector of probabilities of node i at time step t
v[0,0] = 1 #we impose the particle begins at the first vertex/node
p_0 = np.zeros(t) #probability that the particle hits the final vertex for the 1st time at step t or less'
p_0_ex = np.zeros(t) #probability that the particle hits the final vertex for the 1st time at step t or less'
steps = np.arange(0, t)

'''lets define an array s that will contains the position vector of each node.
That is, for D=2, s is a matrix of dimensions 4x2 that contains the nodes positions {(00),(01),(10),(11)}
Each row contains the vector of node i
'''
s = generate_hypercube_nodes(D)

#define de values of p_dir:
#p_dir = generate_random_probabilities(D)
p_dir = np.ones(D)/D #generate equal probabilities for all directions

#lets define the values of matrix P
'''
we will work with a loop that will iterate over all values of the matrix P and will define its values. The values of matrix P are non-zero 
only for elements that represent a movement between consecutive nodes, so we will define a funtion 'displacement' that receives the 
position of 2 nodes s1, s2 and will return two values: a boolean which says if these nodes are consecutive, and 
in the case they are consecutive, it will also return the position/coordinate/direction where the displacement takes place. 
'''

for i in range (0,n):
    for j in range (0,n):
        dir_movement, consecutive = displacement(s[i,:],s[j,:])
        if consecutive == True: 
            P[i,j] = p_dir[dir_movement]

print('Transition matrix: ', P)

'''
Lets upload the vector s[n,:] that contains the probabilities of moving from node (11..1)
They all have to be 0 because we are considering the FIRST arrival to this node, so when 
particle arrive to this node it can't moves anymore
That is, P[n,:] =[0,0,...,0,1]
'''
P[n-1,:] = 0
P[n-1,n-1] = 1

'''print(P)
print(p_dir)'''
        

'''
Now we finally estimate the hitting time
'''

file = f"tau_vector_D={D}_t={t}.txt"
with open(file, 'w') as file:
    for T in range (1,t):
        v[T,:] = v[T-1,:]@P
        p_0[T] = v[T,n-1] #probability of arriving to (1..1) in T or less steps
        p_0_ex[T] = v[T,n-1] - v[T-1,n-1]  #probability of arriving to (1..1) in exactly T steps
        #print ('t=',T,'p_0_ex[T]=', p_0_ex[T])
        '''in order to avoid an innecessary amount of steps, I will add a condition such that if the probability of 
        arriving to the opposite node in exactly T steps is too small (in particular, I will say less than 0.0001), then 
        the loop will break 
        because the next terms will not either contribute to the sum at all. This condition is added to the condition T>D because
        it is the minimum number of steps that are required to arrive to the opposite node
        if (T>D) and  (p_0_ex[T] + p_0_ex[T-1]) < 0.0001: 
            break
        '''
        if int(T%100) == 0:
            print('T = ',T,'p_0_ex[T] + p_0_ex[T-1] = ', p_0_ex[T]+p_0_ex[T-1])
        file.write(f"{T} {T*p_0_ex[T]} \n")

    tau = np.dot(steps,p_0_ex)
    file.write(f"tau = {tau}\n")


print('D = ',D,'t = ',t,'tau = ',tau)
p_0_graphics(p_0_ex,D,tau)

#print('steps=',steps,'p_0=',p_0,'tau=',tau)