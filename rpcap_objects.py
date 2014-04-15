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
        RelPerm.rp_type = rp_type
        return 0
    def get_rp_type(self):
        return RelPerm.rp_type

    def set_sw_res(self, sw_res):
        RelPerm.sw_res = sw_res
        return 0

    def set_viscosities(self, mu_w, mu_n):
        RelPerm.mu_w = mu_w
        RelPerm.mu_n = mu_n
        return 0

    def get_viscosities(self):
        return RelPerm.mu_w, RelPerm.mu_n

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
    def kr_non(self, sw):
        sn = 1. - sw
        if sn <= 0.:
            krn = 0.
        elif sn >= 1.:
            krn = 1.
        else:
            krn = (sn) / (1.)
        return krn
    def d_ds_kr_non(self, sw):
        return -1.
    def fracflow_wet(self, sw):
        """ MUST ADD VISCOSITIES FIRST
        """
        muw, mun = self.get_viscosities()
        krw = self.kr_wet(sw)
        krn = self.kr_non(sw)
        ff = (krw / muw) / (krw / muw + krn / mun)
        return ff
    def fracflow_prime_wet(self, sw):
        muw, mun = self.get_viscosities()
        krw = self.kr_wet(sw)
        krn = self.kr_non(sw)
        high = (krw / muw) 
        low = (krw / muw + krn / mun)
        d_high = self.d_ds_kr_wet(sw) / muw
        d_low = self.d_ds_kr_wet(sw) / muw +\
                self.d_ds_kr_non(sw) / mun
        ffp = (low * d_high - high * d_low) / pow(low, 2.)
        return ffp


    def fracflow_non(self, sw):
        """ MUST ADD VISCOSITIES FIRST
        """
        muw, mun = self.get_viscosities()
        krw = self.kr_wet(sw)
        krn = self.kr_non(sw)
        ff = (krn / mun) / (krw / muw + krn / mun)
        return ff

    def plot_value(self, n_spaces=100., fmt = 'png'):
        sw_res = self.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
        kr_n = np.zeros(len(sw))
        kr_w = np.zeros(len(sw))
        for i in range(len(sw)):
            kr_w[i] = self.kr_wet(sw[i])
            kr_n[i] = self.kr_non(sw[i])
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

    def plot_fracflow(self, n_spaces=100., fmt = 'png'):
        sw_res = self.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
        ff = np.zeros(len(sw))
        for i in range(len(sw)):
            ff[i] = self.fracflow_wet(sw[i])
        fig = plt.figure(num=None, dpi=480,\
                facecolor='w', edgecolor = 'k')
        fig.suptitle('Fractional Flow Function')
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wetting Fluid Saturation []')
        ax.set_ylabel('f(sw) []')
        p1 = plt.plot(sw, ff)
        rp_type = self.get_rp_type()
        fig.savefig('fracflow_wet'+ '_' + str(rp_type) + '.' + fmt)
        plt.clf()
        plt.close()

    def plot_derivative(self, n_spaces=100., fmt = 'png'):
        sw_res = self.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
        dkr_n = np.zeros(len(sw))
        dkr_w = np.zeros(len(sw))
        for i in range(len(sw)):
            dkr_w[i] = self.d_ds_kr_wet(sw[i])
            dkr_n[i] = self.d_ds_kr_non(sw[i])
        fig = plt.figure(num=None, dpi=480,\
                facecolor='w', edgecolor = 'k')
        fig.suptitle('Relative Permeability Derivatives')
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wetting Fluid Saturation []')
        ax.set_ylabel('Rel Perm Derivative []')
        p1 = plt.plot(sw, dkr_w, label = 'wetting')
        p2 = plt.plot(sw, dkr_n, label = 'nonwetting')
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
        CapPres.cp_type = cp_type
        return 0
    def get_cp_type(self):
        return CapPres.cp_type

    def set_sw_res(self, sw_res):
        RelPerm.sw_res = sw_res
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
class Sens:
    def __init__(self, rel_perm_object, cap_pres_object,\
            permeability = 2.e-12):
        self.rel_perm = rel_perm_object
        self.cap_pres = cap_pres_object
        self.k = permeability
        self.g = 9.81
    def set_density_viscosity(self, mu_w, mu_n, rho_w, rho_n):
        self.rel_perm.set_viscosities(mu_w, mu_n)
        self.rho_w = rho_w
        self.rho_n = rho_n
    def q_total(self, sw):
        rho = 600.
        q_total = 0.1 * pow(10.,9.) / (50 * 50 * 24 * 365.25 * 3600) * 1./ rho
        q_total = 0.001
        print "q_total", q_total
        return q_total
    def q_prime(self, sw):
        q_prime = 0.
        return q_prime
    def advection_term(self, sw):
        fw = self.rel_perm.fracflow_wet(sw)
        qt = self.q_total(sw)
        lamb_n = self.rel_perm.kr_non(sw) / self.rel_perm.mu_w
        del_rho_g = (self.rho_w - self.rho_n) * self.g
        F_w = fw * (qt + self.k * lamb_n * del_rho_g)
        return F_w
    def advection_term_prime(self, sw):
        fw = self.rel_perm.fracflow_wet(sw)
        fp = self.rel_perm.fracflow_prime_wet(sw)
        qt = self.q_total(sw)
        qtp = self.q_prime(sw)
        lamb_n = self.rel_perm.kr_non(sw) / self.rel_perm.mu_n
        lamb_n_p = self.rel_perm.d_ds_kr_non(sw) / self.rel_perm.mu_n
        del_rho_g = (self.rho_w - self.rho_n) * self.g
        term1 = fp * (qt + self.k * lamb_n * del_rho_g)
        term2 = fw * (qtp + self.k * del_rho_g * lamb_n_p)
        print "term1,", term1, "term2", term2
        F_w = term1 + term2
        return F_w
    def diffusion_term(self, sw):
        fw = self.rel_perm.fracflow_wet(sw)
        lamb_n = self.rel_perm.kr_non(sw) / self.rel_perm.mu_n
        dpc_dsw = self.cap_pres.d_dsw_pcap(sw)
        D = - fw * self.k * lamb_n * dpc_dsw
        print "term3"
        print D
        return D
    def pseudo_peclet_term(self, sw, delta_x):
        F_w = self.advection_term_prime(sw)
        D = self.diffusion_term(sw)
        pec = F_w * delta_x / D
        return pec
    def plot_diffusion_term(self, n_spaces=100., fmt = 'png'):
        sw_res = self.rel_perm.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
        diffusion = np.zeros(len(sw))
        for i in range(len(sw)):
            diffusion[i] = self.diffusion_term(sw[i])
        fig = plt.figure(num=None, dpi=480,\
                facecolor='w', edgecolor = 'k')
        fig.suptitle("diffusion term D")
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wetting Fluid Saturation []')
        ax.set_ylabel('D []')
        p1 = plt.plot(sw, diffusion)
        rp_type = self.rel_perm.get_rp_type()
        cp_type = self.cap_pres.get_cp_type()
        fig.savefig('diffusion_term' + \
                '_cp_' + str(cp_type) + '_rp_' + str(rp_type) + '.' + fmt)
        plt.clf()
        plt.close()

    def plot_peclet_term(self, delta_x, n_spaces=100., fmt = 'png'):
        sw_res = self.rel_perm.get_sw_res()
        sw = np.linspace(sw_res, 1., n_spaces)
        peclet = np.zeros(len(sw))
        for i in range(len(sw)):
            peclet[i] = self.pseudo_peclet_term(sw[i], delta_x)
        fig = plt.figure(num=None, dpi=480,\
                facecolor='w', edgecolor = 'k')
        fig.suptitle("Peclet term F_w' * dx / D, dx = " + str(delta_x))
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wetting Fluid Saturation []')
        ax.set_ylabel('peclet(sw) []')
        p1 = plt.plot(sw, peclet)
        rp_type = self.rel_perm.get_rp_type()
        cp_type = self.cap_pres.get_cp_type()
        fig.savefig('peclet_term'+ '_dx_' + str(int(delta_x)) + \
                '_cp_' + str(cp_type) + '_rp_' + str(rp_type) + '.' + fmt)
        plt.clf()
        plt.close()

# ===========================================================================
# END ---------------------BASE CLASSES
# ===========================================================================

class RPVanGenuchten(RelPerm):
    def __init__(self, sw_res, lamb, s_lr, s_ls, s_gr):
        self.set_sw_res(sw_res)
        self.set_rp_type('van_genuchten')
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
        if sw == 1.:
            dkrw_ds = self.d_ds_kr_wet(sw-0.001)
        return dkrw_ds

    def d_ds_kr_non(self, sw):
        if self.s_gr == 0:
            dkrn_ds = - dkrw_ds
        elif self.s_gr > 0.:
            sh = (sw - self.s_lr) / (1 - self.s_lr - self.s_gr)
            shp = 1. / (1 - self.s_lr - self.s_gr)
            a = pow(1. - sh, 2.)
            b = 1 - pow(sh, 2.)
            ap = -2. * (1. - sh) * shp
            bp = -2. * sh * shp
            dkrn_ds = a * bp + b * ap
        return dkrn_ds

class RPCubic(RelPerm):
    def __init__(self):
        self.set_sw_res(0.)
        self.set_rp_type('cubic')
    def kr_wet(self, sw):
        krw = pow(sw, 3.)
        return krw
    def d_ds_kr_wet(self, sw):
        dkrw_ds = 3. * pow(sw, 2.)
        return dkrw_ds
    def kr_non(self, sw):
        krn = pow(1. - sw, 3.)
        return krn
    def d_ds_kr_non(self, sw):
        drkn_ds = 3 * pow(1. - sw, 2.) * -1.
        return drkn_ds

class RPGrant(RelPerm):
    def __init__(self, slr, sgr):
        self.set_sw_res(0.)
        self.s_lr = slr
        self.s_gr = sgr
        self.set_rp_type('grant')
    def kr_wet(self, sw):
        sh = (sw - self.s_lr) / (1. - self.s_lr - self.s_gr)
        krw = pow(sh, 4.)
        return krw
    def d_ds_kr_wet(self, sw):
        sh = (sw - self.s_lr) / (1. - self.s_lr - self.s_gr)
        shp = 1. / (1. - self.s_lr - self.s_gr)
        dkrw_ds = 4 * pow(sh, 3.) * shp
        return dkrw_ds
    def kr_non(self, sw):
        krn = 1. - self.kr_wet(sw)
        return krn
    def d_ds_kr_non(self, sw):
        drkn_ds = -1. * self.d_ds_kr_wet(sw)
        return drkn_ds

class CPLeverett(CapPres):
    def __init__(self, s_lr, p0, sigma):
        self.s_lr = s_lr
        self.p0 = p0
        self.sigma = sigma
        self.set_cp_type('leverett')
    def pcap(self, sw):
        ss = (sw - 1.) / (1. - self.s_lr)
        f = 1.417 * (1 - ss) - 2.120 * pow(1 - ss, 2.) + \
                1.263 * pow(1 - ss, 3.)
        pcap = self.p0 * self.sigma * f
        return pcap
    def d_dsw_pcap(self, sw):
        ss = (1. - sw) / (1. - self.s_lr)
        return dpcapdsw
class CPVanGenuchten(CapPres):
    def __init__(self, sw_res, lamb, s_lr, p_0, p_max, s_ls):
        self.set_sw_res(sw_res)
        self.set_cp_type('van_genuchten')
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

    mu_w = 6.9e-4 # Pa*s
    mu_n = 5.45e-5
    rho_w = 1020
    rho_n = 688

    perm = 2.e-12

    rp_lamb = 0.8
    rp_s_lr = 0.2
    rp_s_ls = 1.0
    s_gr = 0.05
    rp_linear = RelPerm(rp_s_lr)

    #rp_van_genuchten = RPVanGenuchten(rp_s_lr, rp_lamb, rp_s_lr, rp_s_ls, s_gr)
    #rp_van_genuchten.plot_value()
    #rp_van_genuchten.plot_derivative()

    #rp_grant = RPGrant(rp_s_lr, s_gr)
    #rp_grant.plot_value()
    #rp_grant.plot_derivative()

    rp_cubic = RPCubic()
    rp_cubic.plot_value()
    rp_cubic.plot_derivative()
    
    cp_lamb = 0.4
    cp_s_lr = 0.0
    cp_p_0 = 1. / 1.61e-3
    cp_p_max = 1.e5
    cp_s_ls = 1.

    cp_pres_constant = CapPres(rp_s_lr)

    #cp_leverett = CPLeverett(0.2, 1.e4, 0.027)
    #cp_leverett.plot_value()

    cp_van_genuchten = CPVanGenuchten(cp_s_lr, cp_lamb, cp_s_lr, cp_p_0,\
            cp_p_max, cp_s_ls)
    cp_van_genuchten.plot_value()
    cp_van_genuchten.plot_derivative()

    s = Sens(rp_cubic, cp_van_genuchten, permeability = perm)
    s.set_density_viscosity(mu_w, mu_n, rho_w, rho_n)
    delta_x = 1.
    s.plot_peclet_term(delta_x)
    s.plot_diffusion_term()
    s.rel_perm.plot_fracflow()
