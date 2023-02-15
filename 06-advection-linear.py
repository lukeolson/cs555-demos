"""
u_t + u_x = 0 on [a,b]
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

method = 'constant'
method = 'linear'
#sigtype = 'LaxW'
#sigtype = 'BW'
sigtype = 'Fromm'
#sigtype = 'minmod'
#sigtype = 'superbee'

T = 1.0
gamma = 0.95
nx = 128

x, hx = np.linspace(0, 1, nx, endpoint=False, retstep=True)
xx = np.linspace(0, 1-hx, 1000)

ht = hx * gamma
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

plt.ion()
fig, ax = plt.subplots()

uline, = ax.plot(x, u, '-', linewidth=3)
txt = ax.text(1.0 / 3.0, 0.8, 't=%g, i=%g' % (0.0, 0), fontsize=16)

def minmod(a, b):
    c = 0 * a
    for i in range(len(a)):
        if a[i] * b[i] > 0:
            if abs(a[i]) <= abs(b[i]):
                c[i] = a[i]
            else:
                c[i] = b[i]
        else:
            c[i] = 0.0
    return c


def maxmod(a, b):
    c = 0 * a
    for i in range(len(a)):
        if a[i] * b[i] > 0:
            if abs(a[i]) >= abs(b[i]):
                c[i] = a[i]
            else:
                c[i] = b[i]
        else:
            c[i] = 0.0
    return c


TV = []
for n in range(0, nt):

    if method == 'constant':
        u[J] = u[J] - (ht / hx) * (u[J] - u[Jm1])

    if method == 'linear':

        if sigtype == 'LaxW':
            sig = (u[Jp1] - u[J]) / hx

        if sigtype == 'BW':
            sig = (u[J] - u[Jm1]) / hx

        if sigtype == 'Fromm':
            sig = (u[Jp1] - u[Jm1]) / (2.0 * hx)

        if sigtype == 'minmod' or sigtype == 'superbee':
            sigm = (u[J] - u[Jm1]) / hx
            sigp = (u[Jp1] - u[J]) / hx

            if sigtype == 'minmod':
                sig = minmod(sigm, sigp)

            if sigtype == 'superbee':
                sig1 = minmod(sigm, 2 * sigp)
                sig2 = minmod(2 * sigm, sigp)
                sig = maxmod(sig1, sig2)

        u[J] = u[J] - (ht / hx) * (u[J] - u[Jm1])\
                    - (ht / (2.0 * hx)) * (sig[J] - sig[Jm1]) * (hx - ht)
    uline.set_ydata(u)
    txt.set_text('t=%g, i=%g' % ((n + 1) * ht, n))
    ax.axis([0, 1, -0.25, 1.25])

    fig.canvas.draw()
    fig.canvas.flush_events()
plt.show(block=True)
