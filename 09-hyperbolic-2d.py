import numpy as np
import matplotlib.pyplot as plt

nx = 128  # number of cells
ny = nx  # assume equal

# domain on -1 to 1
x1d, hx = np.linspace(-1, 1, nx, endpoint=False, retstep=True)
hy = hx
x1d += hx/2
x, y = np.meshgrid(x1d, x1d)

def u0a(x, y):
    u0 = 0 * x
    for i in range(nx):
        for j in range(ny):
            if x[i, j] >= 0.1 and x[i, j] <= 0.6\
                    and y[i, j] >= -0.25 and y[i, j] <= 0.25:
                u0[i, j] = 1.0
    return u0


def u0b(x, y):
    u0 = 0 * x
    for i in range(nx):
        for j in range(ny):
            r = np.sqrt((x[i, j] + 0.45)**2 + y[i, j]**2)
            if r < 0.35:
                u0[i, j] = 1.0 - r / 0.35
    return u0


def u0c(x, y):
    return np.exp(-100 * (x**2 + y**2))


b1 = 1.0
b2 = 1.0

K = np.arange(0, nx)
Kp1 = np.roll(K, -1)
Km1 = np.roll(K, 1)

u = u0a(x, y)
#+ u0b(x, y)
#u = u0c(x, y)

fig, ax = plt.subplots()
# um = np.ma.masked_where(np.fabs(u) < 0.001, u)
# uplt = plt.pcolormesh(xhalf, yhalf, u, shading='gouraud')
uplt = ax.pcolor(u)
txt = ax.text(1.0 / 3.0, 0.8, 't=%g, i=%g' % (0.0, 0), fontsize=16, color='w')
plt.axis('equal')
vtkplot = True
noplot = True

if vtkplot:
    from pyvisfile.vtk import write_structured_grid
    mesh = np.rollaxis(np.dstack((x, y)), 2)
    write_structured_grid('output/test0.vts', mesh,
                          point_data=[('u', u[np.newaxis, :, :])])

method = 'DCU'
T = 2.0
lmbda = 0.95
ht = hx * lmbda / 1.0
tsteps = int(T/ht)
nt = tsteps
for i in range(0, nt):
    print("time step %d of %d" % (i, nt))
    if method == 'DCU':
        Fm = b1 * u[Km1][:, K]
        Fp = b1 * u[K][:, K]
        Gm = b2 * u[K][:, Km1]
        Gp = b2 * u[K][:, K]
        u = u[K][:, K] - (ht / hx) * (Fp - Fm) - (ht / hy) * (Gp - Gm)

    if method == 'DCU2step':
        Fm = b1 * u[Km1][:, K]
        Fp = b1 * u[K][:, K]
        Gm = b2 * u[K][:, Km1]
        Gp = b2 * u[K][:, K]
        u = u[K][:, K] - (ht / hx) * (Fp - Fm)
        u = u[K][:, K] - (ht / hy) * (Gp - Gm)

    if method == 'CTU':
        Fm = b1 * u[Km1][:, K]
        Fp = b1 * u[K][:, K]
        Gm = b2 * u[K][:, Km1]
        Gp = b2 * u[K][:, K]
        u = u[K][:, K] - (ht / hx) * (Fp - Fm)\
                       - (ht / hy) * (Gp - Gm)\
                       + (ht**2 * b1 * b2 / (hx * hy)) * (u[K][:, K] - u[K][:, Km1] - u[Km1][:, K] + u[Km1][:, Km1])

    if vtkplot:
        write_structured_grid('output/test%d.vts' % (i+1), mesh,
                              point_data=[('u', u[np.newaxis, :, :])])
    else:
        um = u
        uplt.set_array(um.ravel())
        fig.canvas.draw()
        fig.canvas.flush_events()

plt.show(block=True)
err = np.abs(u0c(x, y) - u)
l2err = np.sqrt(hx * hy * np.sum(err**2))
print(l2err)
