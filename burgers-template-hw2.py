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

def exact_rarefraction(x,time,u0):
    # IMPLEMENT
    # .....
    return u0

T = 2.0
gamma = 0.95
nx = 128
nx_ghost = 2 # number of ghost cells on each side

x, hx = np.linspace(0, 5, nx + 2*nx_ghost, endpoint=False, retstep=True)
xx = np.linspace(0, 5, 1000, endpoint=False)
# Ghost cell mask
mask =  np.ones(len(x), dtype=bool) 
mask[:2] = mask[-2:] = False

ht = hx * gamma
nt = int(np.ceil(T/ht))
ht = T/nt

print('T = %g' % T)
print('tsteps = %d' % nt)
print('    hx = %g' % hx)
print('    ht = %g' % ht)
print('lambda = %g' % gamma)

K = np.arange(1, nx+2)   # 1, ..., nx+1
Km1 = K-1   # 0,..., nx
Kp1 = K+1   # 2, ..., nx+2

u = step(x)
u0 = u.copy()

def f(u):
    return u**2/2

def fprime(u):
    return u

for n in range(1, nt+1):

    # add code here
    # ...
    
    # u[:] = u - ht/hx * (flux[K]-flux[Km1])  
    # uexact = exact_rarefraction(x, time, u0)

    # Plot Computed and exact solution 
    time = n * ht
    if abs(time-1.) < ht/2 or abs(time-2) < ht/2.: 
        plt.title('t=%g, i=%g' % (n * ht, n))
        plt.plot(x[mask], u[mask], 'r-', linewidth=1, label='approximate')
        #plt.plot(x[mask], uexact[mask], '-.', linewidth=3, label='exact')
        plt.legend(); plt.show()
