import numpy as np
import matplotlib.pyplot as plt
from Ellipsoids import *

plt.rcParams.update({
    "text.usetex": True,
    'text.latex.preamble': r'\usepackage{amsfonts}',
    "font.family": "Helvetica"
})


def ell_f (b, lb, cb):
    blb = np.outer(b, lb)
    return -b + np.sum((cb*cb)*blb/(1-blb), axis=1)




ell0 = Ellipsoid(np.array([[1/4, 0],[0, 1/9]]), np.array([-1,1]))

ell1 = Ellipsoid(np.array([[1/3, -1/5],[-1/5, 1/2]]), np.array([-1.5,1.5]))
ell2 = Ellipsoid(np.array([[1/1, 1/3],[1/3, 2/1]]), np.array([-0.6,-0.5]))
ell3 = Ellipsoid(np.array([[1/2, 0],[0, 12]]), np.array([-1, 2.4]))

ells = [ell1, ell2, ell3]
plt.figure(figsize=[5,3])
for (i, ell) in enumerate(ells):
    lb, cb = ell.get_lb_cb(ell0)

    b = np.linspace(1/min(lb), 1.3, 500)
    ellb = ell_f(b, lb, cb)
    ellb[np.abs(ellb)>10] = np.nan


    plt.subplot(2, len(ells), i+1)
    h_ell = plt.plot(*ell.plot_data_2d())
    h_ell0 = plt.plot(*ell0.plot_data_2d())
    ax = plt.gca()
    xlim = [-4.2, 2.2]
    ylim = [-2.2, 4.2]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if i == 0:
        ax.set_ylabel("Ellipsoids")
    if i == 2:
        plt.plot(np.ones(2)*ell0.c[0], np.array(ylim), linestyle=":", linewidth=0.4, color=[0.6, 0.6, 0.6])
        plt.plot(np.array(xlim), np.ones(2)*ell.c[1], linestyle=":", linewidth=0.4, color=[0.6, 0.6, 0.6])

    plt.subplot(2, len(ells), len(ells)+i+1)    

    h_ellb = plt.plot(b, ellb, color=[0.2, 0.2, 0.6])
    h_none = plt.plot(b,b*0-1, color=[0.7, 0.3, 0.3])
    ax = plt.gca()
    ax.set_xlim([1/min(lb),1.3])
    ax.set_ylim([-2.6,-0.6])
    # ax.grid(linestyle="-", linewidth=0.1)
    if i == 0:
        ax.set_ylabel("$\\ell_{\\tilde{\\mathbf{c}},\\tilde{\\mathbf{P}}}(\\beta)$")
    ax.set_xlabel("$\\beta$")

    print(max(ellb[np.isnan(ellb)==False]))

plt.tight_layout()
plt.savefig('ellipsoids.eps', format='eps')

#plt.show()


