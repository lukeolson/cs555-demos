"""DG for Burgers with limiting."""
import sys
import numpy as np
from scipy.special import legendre
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
import quadpy

#  An outline of this code is a follows:
#      1. problem defition
#      2. mesh setup
#      3. definition of f, the numerical flux, and the right hand side
#      4. the slop limiter
#      5. some random plotting routines
#      4. the main time stepping loop

# --------------------------------------------------------------------
# 1
# --------------------------------------------------------------------
# Define the problem
# Periodic boundary conditions
xa, xb = (0, 1)        # domain
nel = 100              # number of elements
p = 3                  # poly order
applylimiter = True    # should we slope limit?
T = 0.5                # length of time
ht = 0.001            # time step size


def uinit(z):
    """Set the initial condition.

    z : ndarray
        Assume z is p+1 x nel
    """
    uu = np.zeros_like(z)
    for i in range(uu.shape[1]):
        if np.all(z[:, i] <= 0.5):
            # make sure the entire cell is less than 0.5
            uu[:, i] = 1.0
    J = np.where(z <= 0.25)
    uu[J] = 4 * z[J]
    return uu


# --------------------------------------------------------------------
# 2
# --------------------------------------------------------------------
# Set the mesh points
x = np.linspace(xa, xb, nel + 1)
h = x[1] - x[0]

# Construct the local quadrature points
scheme = quadpy.c1.gauss_lobatto(p + 1)
qx, qw = scheme.points, scheme.weights

# Full mesh, all nodes
xx = np.zeros((p + 1, nel))
for e in range(nel):
    xx[:, e] = x[e] + (h / 2) * (qx + 1)


# Vandermonde / mass
def nlegendre(qorder):
    """Compute the Legendre with unite L2."""
    return legendre(qorder) / np.sqrt(2/(2*qorder+1))


V = np.zeros((p+1, p+1))      # Vandermonde
dV = np.zeros((p+1, p+1))     # derivative Vandermonde
for q in range(p+1):
    poly = nlegendre(q)
    dpoly = poly.deriv()
    V[:, q] = poly(qx)
    dV[:, q] = dpoly(qx)
Vinv = np.linalg.inv(V)
M = Vinv.T @ Vinv            # Mass
Minv = (2/h) * V @ V.T       # Mass inverwse
S = M @ dV @ Vinv            # S (see text)
D = dV @ Vinv                # Derivative

n = (p+1)*nel                # total number of points

#   v-left-ext      v-right-ext
# ---|-------------|-----------|
#     ^left-int   ^right-int

# set the vectorized indices
I = np.arange(nel)           # index i
Im1 = np.roll(I, 1)          # index i-1
Ip1 = np.roll(I, -1)         # index i+1

# --------------------------------------------------------------------
# 3
# --------------------------------------------------------------------
# f: function in the conservation law u_t + f(u)_x = 0
#
# nflux: numerical flux
#
# L: right hand side


def f(u):
    """Calculate u**2 / 2."""
    return u**2 / 2


def df(u):
    """Calcualte derivative, u."""
    return u


def nflux(u_int, u_ext, n_int, n_ext):
    """Lax-Fredrichs flux.

    u : nel x 1
    """
    avgf = (f(u_int) + f(u_ext)) / 2
    jumpu = n_int * u_int + n_ext * u_ext
    C = np.max([np.abs(df(u_int)), np.abs(df(u_ext))])
    return avgf + (C/2) * jumpu


def L(u):  # noqa: N802
    """Right-hand side.

    dudt = M^-1 S.T @ f(u) - M^-1 f*(xright)[0,0,...,1] + f*(xleft)
    """
    flux = np.zeros_like(u)
    u_left_int = u[0, I]
    u_left_ext = u[p, Im1]
    u_right_int = u[p, I]
    u_right_ext = u[0, Ip1]
    flux[0, :] = nflux(u_left_int, u_left_ext, -1, 1)     # left
    flux[p, :] = -nflux(u_right_int, u_right_ext, 1, -1)  # right
    Lu = Minv @ S.T @ f(u) + Minv @ flux
    return Lu


# --------------------------------------------------------------------
# 4
# --------------------------------------------------------------------
def minmodv(a, b, c):
    """Vector-based minmod.

    minmod(a,b,c) = { sign(a) min(|a|, |b|, |c|) if sign(a)=sign(b)=sign(c)
                    {                            else 0
    """
    a = a.ravel()
    b = b.ravel()
    c = c.ravel()
    I = np.where(np.abs(np.sign(a) + np.sign(b) + np.sign(c))
                 == 3                          # all three values are the same
                 )[0]
    mm = np.zeros_like(a)
    mm[I] = np.sign(a[I]) * np.vstack((np.abs(a[I]),
                                       np.abs(b[I]),
                                       np.abs(c[I]))).min(axis=0)
    return mm


def slopelimiter(u):
    """Slopelimiters."""
    uhat = Vinv @ u                      # make it modal

    u1 = uhat.copy()
    u1[2:, :] = 0                        # strip to linear
    u1 = V @ u1                          # make it nodal
    ux = (u1[-1, :] - u1[0, :]) / h      # get the slope for each element
    ux = ux.ravel()

    uhat[1:, :] = 0                      # strip to constant (mean)
    ubar = (V @ uhat)[0, :].ravel()      # convert to nodal, and take the left value
    ubarm1 = ubar[Im1]                   # -1
    ubarp1 = ubar[Ip1]                   # +1

    newslope = minmodv(ux, (ubarp1 - ubar)/h, (ubar - ubarm1)/h)

    xbar = (xx[-1, :] + xx[0, :]) / 2    # midpoint of each element

    # make everything p+1 x nel
    bshape = (p+1, nel)                  # broadcast shape
    xcenter = xx - xbar                  # (broadcast) center x to the element
    ubar = np.broadcast_to(ubar, bshape)
    newslope = np.broadcast_to(newslope, bshape)

    # construct the new linear
    newu = ubar + xcenter * newslope
    return newu


def blanklimiter(u):
    """Apply no limiter."""
    return 1.0 * u


# --------------------------------------------------------------------
# 5
# --------------------------------------------------------------------
def plotnodal(xx, u, ax=None):
    """Plot all of the nodal values."""
    if not ax:
        _, ax = plt.subplots()
    ax.plot(xx.ravel(order='F'), u.ravel(order='F'), '-o', ms=3)


def plotmodal(xx, u, u2=None):
    """High fidelity modal."""
    c = Vinv @ u
    m = 40
    v = np.zeros((m, nel))
    zz = np.zeros((m, nel))
    z = np.linspace(-1, 1, m)
    for k in range(nel):
        zz[:, k] = xx[0, k] + (xx[p, k] - xx[0, k])*(z+1)/2
        for q in range(p+1):
            v[:, k] += c[q, k] * nlegendre(q)(z)

    zz = zz.ravel(order='F')
    plt.plot(zz, v.ravel(order='F'), '-o', ms=3)
    if u2 is not None:
        plt.plot(zz, u2(zz), color='tab:red', lw=1)
    plt.show()


# --------------------------------------------------------------------
# 6
# --------------------------------------------------------------------
# initial condition
u = uinit(xx)
u0 = u.copy()
nstep = int(np.ceil(T/ht))
step = 0

if applylimiter:
    limiter = slopelimiter
else:
    limiter = blanklimiter

plt.ion()
fig, ax = plt.subplots()
uline, = ax.plot(xx.ravel(order='F'), u0.ravel(order='F'), '-', color='tab:red', lw=1)
txt = ax.text(1.0 / 3.0, 0.8, 't=%g, i=%g' % (0.0, 0), fontsize=16)

while step * ht < T:

    # SSP 3 time stepping
    u1 = u + ht * L(u)
    u1 = limiter(u1)
    u2 = (1/4) * (3 * u + u1 + ht * L(u1))
    u2 = limiter(u2)
    u[:] = (1/3) * (u + 2 * u2 + 2 * ht * L(u2))
    u = limiter(u)

    step += 1

    #color = adjust_lightness('tab:blue', amount=amount[i])
    uline.set_ydata(u.ravel(order='F'))
    #uline.set_color('tab:blue')
    txt.set_text('t=%g, i=%g' % ((n + 1) * ht, n))

    fig.canvas.draw()
    fig.canvas.flush_events()
plt.show(block=True)

