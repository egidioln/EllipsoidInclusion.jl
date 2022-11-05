import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    'text.latex.preamble': r'\usepackage{amsfonts}',
    "font.family": "Helvetica"
})


def ell_f (b, lb, cb):
    blb = np.outer(b, lb)
    return -b + np.sum((cb*cb)*blb/(1-blb), axis=1)


class Ellipsoid:
    def __init__(self,P,c):
        self.P = P
        self.c = c
        self.L = np.linalg.cholesky((P+P.transpose())/2)
        
    def get_lb_cb(self, ell0):
        L0 = ell0.L
        L0i = np.linalg.inv(L0)
        P = L0i*self.P*L0i.transpose()
        specDecomp = np.linalg.svd(P)
        lb = specDecomp[1]
        cb = specDecomp[0].transpose() @ L0.transpose() @ (self.c -ell0.c)
        return (lb,cb)

    def plot_data_2d(self, points=100):
        if len(self.c) != 2:
            raise Exception("only for 2D ellipses")
        t = np.linspace(0, 2*np.pi, 100)
        return  np.linalg.inv(self.L.transpose()) @ np.array([np.cos(t), np.sin(t)]) + self.c.reshape([2,1])



ell0 = Ellipsoid(np.array([[1/4, 0],[0, 1/9]]), np.array([-1,1]))

ell1 = Ellipsoid(np.array([[1/3, -1/5],[-1/5, 1/2]]), np.array([-1.5,1.5]))
ell2 = Ellipsoid(np.array([[1/1, 1/3],[1/3, 2/1]]), np.array([-0.6,-0.5]))
ell3 = Ellipsoid(np.array([[1/2.56, 0],[0, 12]]), np.array([-1, 2.4]))

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
    ax.set_xlim([-4.2, 2.2])
    ax.set_ylim([-2.2, 4.2])
    if i == 0:
        ax.set_ylabel("Ellipsoids")

    plt.subplot(2, len(ells), len(ells)+i+1)    

    h_ellb = plt.plot(b, ellb, color=[0.2, 0.2, 0.6])
    h_none = plt.plot(b,b*0-1, color=[0.7, 0.3, 0.3])
    ax = plt.gca()
    ax.set_xlim([1/min(lb),1.3])
    ax.set_ylim([-1.6,-0.6])
    ax.grid(linestyle="-", linewidth=0.1)
    if i == 0:
        ax.set_ylabel("$\\ell_{\\mathbf{c},\\mathbf{P}}(\\beta)$")
    ax.set_xlabel("$\\beta$")

    print(max(ellb[np.isnan(ellb)==False]))

plt.tight_layout()
plt.savefig('ellipsoids.eps', format='eps')

#plt.show()


