"""Solve
u_t + (u^2 / 2)_x = 0 on [0, 5]
with fixed boundary conditions.
"""
import numpy as np
import matplotlib.pyplot as plt


def gaussian(x):
    u = np.exp(-100 * (x - 0.25)**2)
    return u

def step(x):
    u = np.zeros(x.shape)
    for j in range(len(x)):
        if (x[j] >= 0.6) and (x[j] <= 0.8):
            u[j] = 1.0
    return u

def exact_rarefraction(...):
    # IMPLEMENT
    # .....
    return u

T = 2.0
gamma = 0.95
nx = 128


x, hx = np.linspace(0, 5, nx, endpoint=False, retstep=True)
# Ghost cell mask: pretend first and last DoF is a ghost cell
mask =  np.ones(len(x), dtype=bool) 
mask[:1] = mask[-1:] = False
# Indexing arrays
K = np.arange(0, nx)   # 0, ..., nx-1
Km1 = np.roll(K, 1)    # nx-1, 0, 1, ..., nx-2
Kp1 = np.roll(K, -1)   # 1, ..., nx

ht = hx * gamma
nt = int(np.ceil(T/ht))
ht = T/nt

print('T = %g' % T)
print('tsteps = %d' % nt)
print('    hx = %g' % hx)
print('    ht = %g' % ht)
print('lambda = %g' % gamma)

u = step(x)
u0 = u.copy()

def f(u):
    return u**2/2

def fprime(u):
    return u

for n in range(1, nt+1):

    # add code here
    # ...
    
    # Which u values need to be updated?
    # u = u - ht/hx * (flux[K]-flux[Km1])  
    # uexact = exact_rarefraction(x, time, u0)

    # Plot Computed and exact solution 
    time = n * ht
    if abs(time-1.) < ht/2 or abs(time-2) < ht/2.: 
        plt.title('t=%g, i=%g' % (n * ht, n))
        plt.plot(x[mask], u[mask], 'r-', linewidth=1, label='approximate')
        #plt.plot(x[mask], uexact[mask], '-.', linewidth=3, label='exact')
        plt.legend(); plt.show()
