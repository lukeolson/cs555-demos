"""
u_t + (c(x) u)_x = 0 on [a,b]
with u(b) = u(a)
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


def c(x):
    return 1.0 + 0.0 * x
    #return 2.0 + np.sin(2 * np.pi * x)

T = 1.0
gamma = 0.95
nx = 128

x, hx = np.linspace(0, 1, nx, endpoint=False, retstep=True)
xx = np.linspace(0, 1-hx, 1000)

ht = hx * gamma / c(x).max()
nt = int(np.ceil(T/ht))
ht = T/nt

print('T = %g' % T)
print('tsteps = %d' % nt)
print('    hx = %g' % hx)
print('    ht = %g' % ht)
print('lambda = %g' % gamma)

J = np.arange(0, nx)  # all vertices
Jm1 = np.roll(J, 1)
Jp1 = np.roll(J, -1)

u = gaussian(x) + step(x)
u0 = u.copy()

plt.ion()
fig, ax = plt.subplots()
ax.plot(x, u0, 'r-', linewidth=1)
uline, = ax.plot(x, u, '-', linewidth=3)
txt = ax.text(1.0 / 3.0, 0.8, 't=%g, i=%g' % (0.0, 0), fontsize=16)


method = 'constant'

for n in range(1, nt+1):

    ci = c(x[Jp1])  # at x = i + 1/2

    flux = ci * (u[J] + u[Jp1]) - np.abs(ci) * (u[Jp1] - u[J])
    flux /= 2
    u[J] = u[J] - (ht / hx) * (flux[J] - flux[Jm1])

    uline.set_ydata(u)
    txt.set_text('t=%g, i=%g' % (n * ht, n))
    ax.axis([0, 1, -0.25, 1.25])
    fig.canvas.draw()
    fig.canvas.flush_events()

plt.show(block=True)
