import numpy as np
import matplotlib.pyplot as plt

class RelPerm:
    """ Relative Permeability class, 
        Base class calls linear rel. perm
    """
    def __init__(self, sw_res, rp_type = 'linear'):
        self.rp_type = rp_type
        self.sw_res = sw_res
    def kr_wet(self, sw):
        if sw <= self.sw_res:
            krw = 0.
        elif sw >= 1.:
            krw = 1.
        else:
            krw = (sw - self.sw_res) / (1. - self.sw_res)
        return krw
    def kr_non(self, sw):
        sn = 1. - sw
        if sn <= 0.:
            krn = 0.
        elif sn >= 1.:
            krn = 1.
        else:
            krn = (sn) / (1.)
        return krn
    def plot(self, n_spaces=100., fmt = 'png'):
        sw = np.linspace(self.sw_res, 1., n_spaces)
        kr_n = np.zeros(len(sw))
        kr_w = np.zeros(len(sw))
        for i in range(len(sw)):
            kr_n[i] = self.kr_wet(sw[i])
            kr_w[i] = self.kr_non(sw[i])
        fig = plt.figure(num=None, dpi=480,\
                facecolor='w', edgecolor = 'k')
        fig.suptitle('Relative Permeability Curves')
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wetting Fluid Saturation []')
        ax.set_ylabel('Relative Permeability []')
        p1 = plt.plot(sw, kr_w, label = 'wetting')
        p2 = plt.plot(sw, kr_n, label = 'nonwetting')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor= (1, 0.5))
        fig.savefig('rel_perm_curves_' + str(self.rp_type) + '.' + fmt)
        plt.clf()
        plt.close()

class CapPres:
    def __init__(self, sw_res, cp_type = 'constant'):
        self.sw_res = sw_res
        self.cp_type = cp_type
    def pcap(self, sw):
        pcap = 0.
        return pcap
    def plot(self, n_spaces=100., fmt = 'png'):
        sw = np.linspace(self.sw_res, 1., n_spaces)
        pc = np.zeros(len(sw))
        for i in range(len(sw)):
            pc[i] = self.pcap(sw[i])
        fig = plt.figure(num=None, dpi=480,\
                facecolor='w', edgecolor = 'k')
        fig.suptitle('Capillary Pressure Curve')
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wetting Fluid Saturation []')
        ax.set_ylabel('Capillary Pressure [Pa]')
        p1 = plt.plot(sw, pc, label = 'wetting')
        fig.savefig('capillary_pressure_' + str(self.cp_type) + '.' + fmt)
        plt.clf()
        plt.close()
if __name__ == '__main__':
    sw_res = 0.2
    rp_type = 'linear'
    rp_linear = RelPerm(sw_res, rp_type = rp_type)
    rp_linear.plot()
    cp_type = 'constant'
    cp_constant = CapPres(sw_res, cp_type = cp_type)
    cp_constant.plot()
