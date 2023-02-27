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

#method = 'constant'

method = 'linear'
sigtype = 'LaxW'
#sigtype = 'BW'
#sigtype = 'Fromm'
#sigtype = 'minmod'
#sigtype = 'superbee'

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

K = np.arange(0, nx)  # all vertices
Km1 = np.roll(K, 1)
Kp1 = np.roll(K, -1)

u = gaussian(x) #+ step(x)

plt.ion()
fig, axs = plt.subplots(nrows=2)
ax = axs[0]
axtvd = axs[1]

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

def calctv(u):
    """assumes periodic"""
    return np.abs(u[1:] - u[:-1]).sum() + abs(u[0] - u[-1])

uline, = ax.plot(x, u, '-', linewidth=3)
txt = ax.text(1.0 / 3.0, 0.8, 't=%g, i=%g' % (0.0, 0), fontsize=16)

tv = np.full(nt, np.nan)
tv[0] = calctv(u)
tvdline, = axtvd.plot(tv, 'o', color='tab:red', markerfacecolor='w', ms=3)

for n in range(1, nt+1):

    if method == 'constant':
        ci = c(x[Kp1])  # at x = i + 1/2
        fi = 0 * ci   # at x = i + 1/2
        pos = np.where(ci >= 0)[0]
        fi[pos] = ci[pos] * u[K[pos]]
        neg = np.where(ci < 0)[0]
        fi[neg] = ci[neg] * u[Kp1[neg]]

        u[K] = u[K] - (ht / hx) * (fi[K] - fi[Km1])

    if method == 'linear':

        if sigtype == 'LaxW':
            sig = (u[Kp1] - u[K]) / hx

        if sigtype == 'BW':
            sig = (u[K] - u[Km1]) / hx

        if sigtype == 'Fromm':
            sig = (u[Kp1] - u[Km1]) / (2.0 * hx)

        if sigtype == 'minmod' or sigtype == 'superbee':
            sigm = (u[K] - u[Km1]) / hx
            sigp = (u[Kp1] - u[K]) / hx

            if sigtype == 'minmod':
                sig = minmod(sigm, sigp)

            if sigtype == 'superbee':
                sig1 = minmod(sigm, 2 * sigp)
                sig2 = minmod(2 * sigm, sigp)
                sig = maxmod(sig1, sig2)

        u[K] = u[K] - (c(x[K]) * ht / hx) * (u[K] - u[Km1])\
                    - (c(x[K]) * ht / (2.0 * hx)) * (sig[K] - sig[Km1])\
                                                  * (hx - c(x[K]) * ht)
    uline.set_ydata(u)
    txt.set_text('t=%g, i=%g' % ((n + 1) * ht, n))
    ax.axis([0, 1, -0.25, 1.25])

    tv[n-1] = calctv(u)
    tvdline.set_ydata(tv)
    axtvd.axis([0, nt+1, -0.25, 5])

    fig.canvas.draw()
    fig.canvas.flush_events()
plt.show(block=True)
