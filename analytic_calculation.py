'''
Analytic calculation of Hitting time tau
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial, binom
from scipy.optimize import curve_fit
import myconfig

def tau_cube(dim):
    k_arr = np.arange(0,dim,1)
    return 2**(dim-1)*np.sum(1/(binom(dim-1, k_arr)))

def tau_tree(depth):
    return 8*(2**depth - 1)*(1 - 1/(2**depth))

def m_cube(dim):
    return 2**(dim-1)*dim

def m_tree(depth):
    return 4*(2**depth - 1)

def R_cube(dim):
    k_arr = np.arange(0,dim,1)
    return np.sum(1/(binom(dim-1, k_arr)))

def R_tree(depth):
    return 2*(1-1/2**depth) 

def n_cube(dim):
    return 2**dim

def n_tree(depth):
    return (2**(depth+1) - 2) + 2**depth



def general_plot_function_tau(x_arr_cube, x_arr_tree, d_arr, x_label="x"):

    tau_cube_arr = [tau_cube(d) for d in d_arr]
    tau_tree_arr = [tau_tree((d)/2) for d in d_arr if d % 2 == 0]

    tau_tree_arr = np.array(tau_tree_arr)
    tau_cube_arr = np.array(tau_cube_arr)[:len(tau_tree_arr)]
    
    x_arr_cube = np.array(x_arr_cube)[:len(tau_tree_arr)]
    x_arr_tree = np.array(x_arr_tree)

    #lin fit linear function
    fit_cube_lin = np.polyfit(x_arr_cube[:5], tau_cube_arr[:5], 1)
    fit_tree_lin = np.polyfit(x_arr_tree[:5], tau_tree_arr[:5], 1)
    x_fit_cube_lin = np.linspace(0, x_arr_cube[8], 100)
    x_fit_tree_lin = np.linspace(0, x_arr_tree[8], 100)
    y_fit_cube_lin = fit_cube_lin[0] * x_fit_cube_lin + fit_cube_lin[1]
    y_fit_tree_lin = fit_tree_lin[0] * x_fit_tree_lin + fit_tree_lin[1]

    print("LinFit cube: m=",fit_cube_lin[0], "offset=", fit_cube_lin[1])
    print("LinFit tree: m=",fit_tree_lin[0], "offset=", fit_tree_lin[1])

    def fit_func(x, m, offset):
        return m * x + offset
    #lin fit power function
    fit_cube_pow, cov_cube = curve_fit(fit_func, np.log(x_arr_cube[10:]), np.log(tau_cube_arr[10:]))
    fit_tree_pow, cov_tree = curve_fit(fit_func, np.log(x_arr_tree[10:]), np.log(tau_tree_arr[10:]))
    # fit_cube_pow, cov_cube = np.polyfit(np.log(x_arr_cube[10:]), np.log(tau_cube_arr[10:]), 1)
    # fit_tree_pow, cov_tree = np.polyfit(np.log(x_arr_tree[10:]), np.log(tau_tree_arr[10:]), 1)
    x_fit_cube_pow = np.linspace(x_arr_cube[10], max(x_arr_cube), 100)
    x_fit_tree_pow = np.linspace(x_arr_tree[10], max(x_arr_tree), 100)
    # y_fit_cube = np.exp(fit_cube[1]) * x_fit_cube**fit_cube[0]
    # y_fit_tree = np.exp(fit_tree[1]) * x_fit_tree**fit_tree[0]
    y_fit_cube_pow = np.exp(fit_cube_pow[1]) * x_fit_cube_pow**fit_cube_pow[0]
    y_fit_tree_pow = np.exp(fit_tree_pow[1]) * x_fit_tree_pow**fit_tree_pow[0]

    print("PowFit cube: m=",np.exp(fit_cube_pow[1]),r"$\pm$", np.sqrt(cov_cube[1,1]), "exponent=", fit_cube_pow[0],r"$\pm$", np.sqrt(cov_cube[0,0]))
    print("PowFit tree: m=",np.exp(fit_tree_pow[1]),r"$\pm$", np.sqrt(cov_tree[1,1]), "exponent=", fit_tree_pow[0],r"$\pm$", np.sqrt(cov_tree[0,0]))


    plt.figure(0, figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.scatter(x_arr_cube[:6], tau_cube_arr[:6], label='Cube', color='#1f77b4', alpha=1, s=20, zorder=5)
    plt.scatter(x_arr_tree[:6], tau_tree_arr[:6], label='Tree', color='#ff7f0e', alpha=1, s=16, zorder=5)
    plt.plot(x_fit_cube_lin, y_fit_cube_lin, label='Linear fits', color="#C32DFE", alpha=0.5, ls="--", linewidth=1.5, zorder=6)
    plt.plot(x_fit_tree_lin, y_fit_tree_lin, color='#C32DFE', alpha=0.5, ls="--", linewidth=1.5, zorder=6)
    plt.xlim(0, max(x_arr_cube[:6]) * 1.5)
    plt.ylim(min(min(y_fit_cube_lin[:6]), min(y_fit_tree_lin[:6]))*1.1, max(tau_cube_arr[:6]) * 1.1)
    plt.grid(which="major", linestyle='-', linewidth=0.5, color="lightgray",zorder=0)
    plt.xlabel(x_label)
    plt.ylabel(r'Hitting Time $\tau_{0, \text{end}}$')
    plt.legend(loc='upper left')

    plt.subplot(1, 2, 2)
    plt.scatter(x_arr_cube, tau_cube_arr, label='Cube', color='#1f77b4', alpha=1, s=20, zorder=5)
    plt.scatter(x_arr_tree, tau_tree_arr, label='Tree', color='#ff7f0e', alpha=1, s=16, zorder=5)
    plt.plot(x_fit_cube_pow, y_fit_cube_pow, label='Power function fits', color="#C32DFE", alpha=0.9, linewidth=1.5, zorder=6)
    plt.plot(x_fit_tree_pow, y_fit_tree_pow, color='#C32DFE', alpha=0.9, linewidth=1.5, zorder=6)

    # plt.scatter(x_arr_tree, tau_cube_arr[:len(tau_tree_arr)]/tau_tree_arr, label='Cube/Tree')
    # plt.ylabel(r'Hitting Time $\tau_{0, \text{end}}$')
    plt.xlabel(x_label)
    plt.xscale("log")
    plt.yscale("log")

    plt.xticks(np.logspace(0,15, 4))
    plt.xticks(np.logspace(0,15, 16),"", minor=True)
    plt.yticks(np.logspace(0,15, 4))
    plt.yticks(np.logspace(0,15, 16),"", minor=True)

    plt.grid(which="major", linestyle='-', linewidth=0.5, color="lightgray",zorder=0)
    plt.grid(which="minor", linestyle=':', linewidth=0.5, color="lightgray", zorder=0)
    # plt.xlim(0, min(max(x_arr_cube), max(x_arr_tree)) * 1.1)
    # plt.ylim(0, min(max(tau_cube_arr), max(tau_tree_arr)) * 1.1)
    # plt.title('Hitting Time vs ' + x_label)
    plt.legend(loc='upper left')
    # plt.show()


d_table = np.arange(1,14,1)
tau_cube_arr = [tau_cube(d) for d in d_table]
# tau_tree_arr = [tau_tree((d)/2) for d in d_table if d % 2 == 0]

print("d_table:", d_table)
print("tau_cube_arr:", tau_cube_arr)

d_arr = np.arange(1, 100, 1)
n_arr_cube = [n_cube(d) for d in d_arr]
n_arr_tree = [n_tree((d)/2) for d in d_arr if d % 2 == 0]
general_plot_function_tau(n_arr_cube, n_arr_tree, d_arr, x_label='Number of nodes')
# plt.savefig("plots/tau_vs_n.png", bbox_inches='tight', dpi=300)

# d-sweep
# L=2d -> d=(L)/2
d_arr_cube = np.arange(1, 40, 1)
d_arr_tree = np.arange(2, 40, 2)

tau_cube_arr = [tau_cube(d) for d in d_arr_cube]
tau_tree_arr = [tau_tree((d)/2) for d in d_arr_tree]
# tau_tree_arr = [tau_tree(d) for d in d_arr_tree]


# Plotting log scale
plt.figure(2)
plt.scatter(d_arr_cube, tau_cube_arr, label='Cube')
plt.yscale('log')
plt.scatter(d_arr_tree, tau_tree_arr, label='Tree')
# plt.xlabel('dimension d (cube) | length L= 2d+1 (binary tree)')
plt.ylabel(r'Hitting Time $\tau_{0, \text{end}}$')
plt.xlabel(r'Length $L$ (d-cube: $L=d$, GBT: $L=2\gamma$)')
plt.legend(loc='upper left')
# plt.savefig("plots/tau_vs_d.png", bbox_inches='tight', dpi=300)
# plt.legend()

# Plotting m with one plot linear
m_cube_arr = [m_cube(d) for d in d_arr_cube]
m_tree_arr = [m_tree((d)/2) for d in d_arr_tree]
plt.figure(4)
plt.scatter(d_arr_cube, m_cube_arr, label='Cube')
plt.scatter(d_arr_tree, m_tree_arr, label='Tree')
plt.ylabel(r'Number of links $m$')
plt.xlabel(r'Length $L$ (d-cube: $L=d$, GBT: $L=2\gamma$)')
plt.yscale('log')
plt.legend(loc='upper left')
# plt.savefig("plots/m_vs_d.png", bbox_inches='tight', dpi=300)

# Plotting n

n_cube_arr = [n_cube(d) for d in d_arr_cube]
n_tree_arr = [n_tree((d)/2) for d in d_arr_tree]
plt.figure(6)
plt.scatter(d_arr_cube, n_cube_arr, label='Cube')
plt.scatter(d_arr_tree, n_tree_arr, label='Tree')
plt.ylabel(r'Number of nodes $N$')
plt.xlabel(r'Length $L$ (d-cube: $L=d$, GBT: $L=2\gamma$)')
plt.yscale('log')
plt.legend(loc='upper left')
# plt.savefig("plots/nodes_vs_d.png", bbox_inches='tight', dpi=300)

# Plotting resistance with one plot linear and one plot log-log combined in with subplots

R_cube_arr = [R_cube(d) for d in d_arr_cube]
R_tree_arr = [R_tree((d)/2) for d in d_arr_tree]
plt.figure(5)
plt.scatter(d_arr_cube, R_cube_arr, label='Cube')
plt.scatter(d_arr_tree, R_tree_arr, label='Tree')
plt.ylabel(r'Resistance $R_{0, \text{end}}$')
plt.xlabel(r'Length $L$ (d-cube: $L=d$, GBT: $L=2\gamma$)')

plt.legend(loc='right')
# plt.savefig("plots/R_vs_d.png", bbox_inches='tight', dpi=300)
plt.show()