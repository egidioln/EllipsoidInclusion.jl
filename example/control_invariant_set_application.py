import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import scipy.spatial as sp
from Ellipsoids import *

np.random.seed(1)

plt.rcParams.update({
    "text.usetex": True,
    'text.latex.preamble': r'\usepackage{amsfonts}',
    "font.family": "Helvetica"
})


n = 2
m = 1
A = np.array([[0.0, 1],
             [0.1, 0.3]])
B = np.array([[0.0],[0.5]])
H = np.array([[-0.3],[0.6]])
Q = np.eye(n)
R = np.eye(m)

P = la.solve_continuous_are(A, B, Q, R)
K = np.linalg.inv(R) @ B.transpose() @ P
Acl = A - B @ K

wmax = 0.5
PHw = P @ H*wmax
S = (-Acl.transpose() @ P - P @ Acl)
Si = np.linalg.inv(S)
cw = Si @ PHw
Sw = S*( 1 /(cw.transpose() @ S @ cw) )
badSet = Ellipsoid(Sw, cw)
lyapSet = Ellipsoid(P, np.zeros([n,1]))

# Define ell
def build_ell(el,el0):
    lb, cb = el.get_lb_cb(el0)
    cb2 = (cb**2).squeeze()
    def eval_ell(b):
        try:
            blb = b*lb 
            blb_m_1  =  b*lb - 1
            blb_m_1_sq = blb_m_1*blb_m_1
            blb_m_1_cb = blb_m_1_sq*blb_m_1
            polPos =  -b - np.sum((blb / blb_m_1)*cb2)
            dpolPos =  -1 + np.sum((lb / blb_m_1_sq)*cb2)
            ddpolPos =  -2*np.sum((lb**2 /blb_m_1_cb)*cb2)
            return (polPos, dpolPos, ddpolPos)
        except ZeroDivisionError:
            return (-np.infty,)*3
    return (lambda b: eval_ell(b), lb, cb)


eval_ell, lb, cb = build_ell(badSet,lyapSet)

# Maximize ell
def maximize(eval_ell, lb, epsilon = 1e-10, upperbound=1e4):
    interval = [1 / min(lb), upperbound]
    b = None
    f = None
    df = None
    while np.diff(interval)[0]>epsilon:
        b = np.mean(interval)
        f, df,_ = eval_ell(b)
        if df > 0:
            interval[0] = b
        else:
            interval[1] = b
    return (b,f,df)

b, f, df = maximize(eval_ell, lb)

lyapSetIn = Ellipsoid(P*(-1/f), np.zeros([n,1]))

# Simulation
dt = 0.01
T = 20
tspan = np.linspace(0, T, int(T/dt))
Ad = la.expm(Acl*dt)
Hd = la.inv(Acl) @ (Ad - np.eye(n)) @ H

x = np.zeros([n, len(tspan)])
x0 = -1 * np.ones(n)
x[:,0] = x0
w = wmax
for k in range(1,len(tspan)):
    if np.random.rand() > 0.99:
        w *= -1
    x[:,k] = Ad @ x[:,k-1] + (Hd * w).squeeze()
    
bad_data  = badSet.plot_data_2d()
all_bad_data = np.concatenate([bad_data, -bad_data],1)
bad_hull = sp.ConvexHull(all_bad_data.transpose())
hull_vertices = np.append(bad_hull.vertices, bad_hull.vertices[0])
plt.figure(figsize=[3.7,2.5])
h_bad1, = plt.plot(*bad_data, label="$\mathcal{B}_i$")
h_bad2, = plt.plot(*(-bad_data), label='_nolegend_')
h_bhull, = plt.fill(*all_bad_data[:,hull_vertices], label='_nolegend_')
h_in, = plt.plot(*lyapSetIn.plot_data_2d(), label="$\mathcal{V}$")
h_traj, =  plt.plot(*x, label="$x(t)$")
h_x0 = plt.scatter(x0[0],x0[1], label='$x_0$')
h_x0.set_color(h_traj.get_color())
h_bad2.set_color(h_bad1.get_color())
h_bhull.set_color(h_bad1.get_color())
h_bhull.set_color(tuple([(x+2)/3 for x in h_bhull.get_facecolor()]))

ax = plt.gca()
ax.set_xlim([-2, 2])
ax.set_ylim([-1.5, 1.5])
ax.set_xlabel("$x_1$")
ax.set_ylabel("$x_2$")
plt.legend()
#plt.legend([",","])
plt.tight_layout()
plt.savefig(f'control_ex.eps', format='eps')
