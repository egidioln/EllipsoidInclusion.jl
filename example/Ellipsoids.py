import numpy as np

class Ellipsoid:
    def __init__(self,P,c):
        self.P = P
        self.c = c
        self.L = np.linalg.cholesky((P+P.transpose())/2)
        
    def get_lb_cb(self, ell0):
        L0 = ell0.L
        L0i = np.linalg.inv(L0)
        P = L0i @ self.P @ L0i.transpose()
        specDecomp = np.linalg.svd(P)
        lb = specDecomp[1]
        cb = specDecomp[0].transpose() @ L0.transpose() @ (self.c -ell0.c)
        return (lb,cb)

    def plot_data_2d(self, points=100):
        if len(self.c) != 2:
            raise Exception("only for 2D ellipses")
        t = np.linspace(0, 2*np.pi, 100)
        return  np.linalg.inv(self.L.transpose()) @ np.array([np.cos(t), np.sin(t)]) + self.c.reshape([2,1])

    def inside(self, ell0):
        lb, cb = self.get_lb_cb(ell0)
        pass
        