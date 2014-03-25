import numpy as np
import matplotlib.pyplot as plt

# ===========================================================================
# BASE CLASSES
# ===========================================================================
class RelPerm:
    """ Relative Permeability class, 
        Base class calls linear rel. perm
    """
    def __init__(self, sw_res, rp_type = 'linear'):
        RelPerm.rp_type = rp_type
        RelPerm.sw_res = sw_res

    def set_rp_type(self, rp_type):
        print "Changing relative permeability type from " + str(self.rp_type)
        RelPerm.rp_type = rp_type
        print " to " + str(self.rp_type)
        return 0
    def get_rp_type(self):
        return RelPerm.rp_type

    def set_sw_res(self, sw_res):
        print "Changing sw_res from " + str(self.sw_res)
        RelPerm.sw_res = sw_res
        print " to " + str(self.sw_res)
        return 0

    def get_sw_res(self):
        return RelPerm.sw_res

    def kr_wet(self, sw):
        sw_res = self.get_sw_res()
        if sw <= sw_res:
            krw = 0.
        elif sw >= 1.:
            krw = 1.
        else:
            krw = (sw - sw_res) / (1. - sw_res)
        return krw
    def d_ds_kr_wet(self, sw):
        return 1./ (1. - self.get_sw_res())
    def d_ds_kr_non(self, sw):
        return -1.
    def kr_non(self, sw):
        sn = 1. - sw
        if sn <= 0.:
            krn = 0.
        elif sn >= 1.:
            krn = 1.
        else:
            krn = (sn) / (1.)
        return krn
    def plot_value(self, n_spaces=100., fmt = 'png'):
        sw_res = self.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
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
        rp_type = self.get_rp_type()
        fig.savefig('rel_perm_curves_' + str(rp_type) + '.' + fmt)
        plt.clf()
        plt.close()

    def plot_derivative(self, n_spaces=100., fmt = 'png'):
        sw_res = self.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
        kr_n = np.zeros(len(sw))
        kr_w = np.zeros(len(sw))
        for i in range(len(sw)):
            kr_n[i] = self.d_ds_kr_wet(sw[i])
            kr_w[i] = self.d_ds_kr_non(sw[i])
        fig = plt.figure(num=None, dpi=480,\
                facecolor='w', edgecolor = 'k')
        fig.suptitle('Relative Permeability Derivatives')
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wetting Fluid Saturation []')
        ax.set_ylabel('Rel Perm Derivative []')
        p1 = plt.plot(sw, kr_w, label = 'wetting')
        p2 = plt.plot(sw, kr_n, label = 'nonwetting')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor= (1, 0.5))
        rp_type = self.get_rp_type()
        fig.savefig('rel_perm_derivatives' + str(rp_type) + '.' + fmt)
        plt.clf()
        plt.close()

class CapPres:
    """ Base class for capillary pressure.
        Base class with no modification calls zero capillary pressure
    """
    def __init__(self, sw_res, cp_type = 'constant'):
        CapPres.sw_res = sw_res
        CapPres.cp_type = cp_type

    def set_cp_type(self, cp_type):
        print "changing capillary pressure type: " + str(CapPres.cp_type)
        CapPres.cp_type = cp_type
        print " to: " + str(CapPres.cp_type)
        return 0
    def get_cp_type(self):
        return CapPres.cp_type

    def set_sw_res(self, sw_res):
        print "Changing sw_res from " + str(self.sw_res)
        RelPerm.sw_res = sw_res
        print " to " + str(self.sw_res)
        return 0
    def get_sw_res(self):
        return CapPres.sw_res

    def pcap(self, sw):
        pcap = 0.
        return pcap

    def d_dsw_pcap(self, sw):
        dpds = 0.
        return dpds

    def plot_value(self, n_spaces=100., fmt = 'png'):
        sw_res = self.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
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
        cp_type = self.get_cp_type()
        fig.savefig('capillary_pressure_' + str(cp_type) + '.' + fmt)
        plt.clf()
        plt.close()
    def plot_derivative(self, n_spaces=100., fmt = 'png'):
        sw_res = self.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
        pc = np.zeros(len(sw))
        for i in range(len(sw)):
            pc[i] = self.d_dsw_pcap(sw[i])
        fig = plt.figure(num=None, dpi=480,\
                facecolor='w', edgecolor = 'k')
        fig.suptitle('d/dsw (pcap) [Pa]/[]')
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wetting Fluid Saturation []')
        ax.set_ylabel('Capillary Pressure Derivative [Pa]/[]')
        p1 = plt.plot(sw, pc, label = 'wetting')
        cp_type = self.get_cp_type()
        fig.savefig('d_ds_pcap' + str(cp_type) + '.' + fmt)
        plt.clf()
        plt.close()
# ===========================================================================
# END ---------------------BASE CLASSES
# ===========================================================================

class RPVanGenuchten(RelPerm):
    def __init__(self, sw_res, lamb, s_lr, s_ls, s_gr, \
            rp_type = 'van_genuchten'):
        self.set_sw_res(sw_res)
        self.set_rp_type(rp_type)
        self.lamb = lamb
        self.s_ls = s_ls
        self.s_lr = s_lr
        self.s_gr = s_gr
    def kr_wet(self, sw):
        s_star = (sw - self.s_lr) / (self.s_ls - self.s_lr)
        s_hat = (sw - self.s_lr) / (1 - self.s_lr - self.s_gr)
        if sw < self.s_ls:
            inner = 1 - pow(1 - pow(s_star, 1./self.lamb), self.lamb)
            krw = np.sqrt(s_star) * pow(inner, 2.)
        elif sw >= self.s_ls:
            krw = 1.
        return krw
    def kr_non(self, sw):
        if self.s_gr == 0:
            krn = 1 - self.kr_wet
        elif self.s_gr > 0.:
            s_hat = (sw - self.s_lr) / (1 - self.s_lr - self.s_gr)
            krn = pow(1 - s_hat, 2.) * ( 1 - pow(s_hat, 2.))
        return krn
    def d_ds_kr_wet(self, sw):
        """ krw = a(ss) * f(g(h(ss)))
            dkrw/ds = 
        a'(ss) * ss' * f(g(h(sstar)) + 
        a * f'(g(h(ss) * g'(h(ss)) * h'(ss) * ss'
            a = sqrt(sstar)
            f = g^2
            g = 1 - h^lamb
            h = 1 - sstar^(1/lamb)
            the prime multiplication is nested here instead 
            of multiplied at the end. 
        """
        #TODO This does not correctly evaluate for the end, when sw = 1.
        ss = (sw - self.s_lr) / (self.s_ls - self.s_lr)
        ssp = 1. / (self.s_ls - self.s_lr)
        h = 1. - pow(ss, 1./self.lamb)
        hp = -1./self.lamb * pow(ss, 1./self.lamb - 1.) * ssp
        g = 1. - pow(h, self.lamb)
        gp = -self.lamb * pow(h, self.lamb -1.) * hp
        f = pow(g, 2.)
        fp = 2. * g * gp
        a = np.sqrt(ss)
        ap = 1. / (a * 2.) * ssp
        dkrw_ds = ap * f + a * fp
        return dkrw_ds

    def d_ds_kr_non(self, sw):
        if self.s_gr == 0:
            dkrn_ds = - dkrw_ds
        elif self.s_gr > 0.:
            sh = (sw - self.s_lr) / (1 - self.s_lr - self.s_gr)
            shp = 1. / (1 - self.s_lr - self.s_gr)
            h = 1. - pow(sh, 2.)
            hp = -2. * sh * shp
            g = 1. - sh
            gp = -shp
            f = pow(g,2.)
            fp = 2 * g * gp
            dkrn_ds = f * hp + h * fp
        return dkrn_ds

class CPVanGenuchten(CapPres):
    def __init__(self, sw_res, lamb, s_lr, p_0, p_max, s_ls, \
            cp_type = 'van_genuchten'):
        self.set_sw_res(sw_res)
        self.set_cp_type(cp_type)
        self.lamb = lamb
        self.s_lr = s_lr
        self.p_0 = p_0
        self.p_max = p_max
        self.s_ls = s_ls
    def pcap(self, sw):
        s_star = (sw - self.s_lr) / (self.s_ls - self.s_lr)
        pcap = self.p_0 * \
                pow(pow(s_star, -1./self.lamb) - 1., 1. - self.lamb)
        if pcap > self.p_max:
            print "hitting the max capillary pressure"
            pcap = self.p_max
        elif pcap < 0.:
            print "less than zero pressure"
            pcap = 0.
        return pcap

    def d_dsw_pcap(self, sw):
        """ f = pO (g) ^ (1- lamb)
            g = (h)^(-1/lamb) -1
            h = s_star
        """
        h = (sw - self.s_lr) / (self.s_ls - self.s_lr)
        g = pow(h, -1. / self.lamb) - 1
        f = self.p_0 * pow(g, 1.- self.lamb)
        hp = 1. / (self.s_ls - self.s_lr)
        gp = -1./self.lamb * pow(h, -1./self.lamb -1)
        fp = self.p_0 * (1 - self.lamb) * pow(g, -self.lamb)
        dpds = fp * gp * hp
        return dpds

if __name__ == '__main__':
    sw_res = 0.2
    rp_linear = RelPerm(sw_res)
    rp_linear.plot_value()
    cp_constant = CapPres(sw_res)
    cp_constant.plot_value()
    rp_lamb = 0.8
    rp_s_lr = 0.2
    rp_s_ls = 1.0
    s_gr = 0.05
    rp_van_genuchten = RPVanGenuchten(rp_s_lr, rp_lamb, rp_s_lr, rp_s_ls, s_gr)
    rp_van_genuchten.plot_value()
    rp_van_genuchten.plot_derivative()
    cp_lamb = 0.4
    cp_s_lr = 0.2
    cp_p_0 = 1. / 1.61e-5
    cp_p_max = 1.e7
    cp_s_ls = 0.999
    cp_van_genuchten = CPVanGenuchten(cp_s_lr, cp_lamb, cp_s_lr, cp_p_0,\
            cp_p_max, cp_s_ls)
    cp_van_genuchten.plot_value()
    cp_van_genuchten.plot_derivative()

