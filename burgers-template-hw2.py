"""Solve
u_t + (u^2 / 2)_x = 0 on [0, 5]
with periodic boundary conditions.
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


T = 2.0
gamma = 0.95
nx = 128

x, hx = np.linspace(0, 5, nx, endpoint=False, retstep=True)
xx = np.linspace(0, 5, 1000, endpoint=False)

ht = hx * gamma
nt = int(np.ceil(T/ht))
ht = T/nt

print('T = %g' % T)
print('tsteps = %d' % nt)
print('    hx = %g' % hx)
print('    ht = %g' % ht)
print('lambda = %g' % gamma)

K = np.arange(0, nx)   # 0, ..., nx-1
Km1 = np.roll(K, 1)    # nx-1, 0, 1, ..., nx-2
Kp1 = np.roll(K, -1)   # 1, ..., nx

u = step(x)
u0 = u.copy()

plt.ion()
fig, ax = plt.subplots()
ax.plot(x, u0, 'r-', linewidth=1)
uline, = ax.plot(x, u, '-', linewidth=3)
txt = ax.text(1.0 / 3.0, 0.8, 't=%g, i=%g' % (0.0, 0), fontsize=16)

def f(u):
    return u**2/2

def fprime(u):
    return u

for n in range(1, nt+1):

    # add code here
    # ...

    u[:] = u - ht/hx * (flux[K]-flux[Km1])

    uline.set_ydata(u)
    txt.set_text('t=%g, i=%g' % (n * ht, n))
    ax.axis([x.min(), x.max(), -0.25, 1.25])
    fig.canvas.draw()
    fig.canvas.flush_events()

plt.show(block=True)
