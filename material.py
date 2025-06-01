#  ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <olivier.dazel@univ-lemans.fr> and <alan.geslain@univ-lemans.fr> wrote this file.  
#  As long as you retain this notice you can do whatever you want with this stuff. 
#  If we meet some day, and you think this stuff is worth it, you can buy us a beer in return.       
# Olivier DAZEL and Alan GESLAIN
#  ----------------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt

from numpy import pi, sqrt


possible_fluid_models = ["Delany Bazley", "Miki", "Johnson-Champoux-Allard", "Johnson-Lafarge", "Pride-Lafarge"]
possible_solid_models = ["Rigid frame", "Biot", "Limp"]


class Air():

    P = 1.01325e5  # atmospheric Pressure [Pa]
    gamma = 1.400  # polytropic coefficient []
    mu = 0.1839e-4  # dynamic viscosity [kg.m^-1.s^-1]
    Pr = 0.710  # Prandtl's number []
    rho = 1.213  # density [kg.m^-3]
    K = gamma*P  # adiabatic bulk modulus
    c = np.sqrt(K/rho)  # adiabatic sound speed
    Z = rho*c  # characteristic impedance
    nu = mu/rho  # kinematic viscosity [m.s^-2]
    nu_prime = nu/Pr  # viscothermal losses
    C_p = 1006  # (mass) specific heat capacity as constant pressure [J.K^-1]
    lambda_ = 0.0262  # thermal conductivity [W.m^-1.K^-1]
    

class Material():
    def __init__(self, params):
        self.update_from_params(params)

    def update_from_params(self, params):
        self.fluid_model = params['fluid_model']
        self.solid_model = params['solid_model']
        self.phi = params['Porosity'].value
        self.sigma = params['Resistivity'].value
        self.alpha_inf = params['Tortuosity'].value
        self.Lambda = params['Viscous characteristic length'].value / 1e6
        self.Lambda_prime = params['Thermal characteristic length'].value / 1e6
        self.k0prime= params['Static thermal permeability'].value
        self.rho_1 = params['Frame density'].value
        self.E = params['Young\'s modulus'].value
        self.nu = params['Poisson\'s ratio'].value
        self.eta_s = params['Loss factor'].value
        self.alpha_0 = params['Static tortuosity'].value

    def equivalent_fluid_parameters(self, f):
        if self.fluid_model == "Delany Bazley":
            self.delany_bazley(f)
        elif self.fluid_model == "Miki":
            self.miki(f)
        elif self.fluid_model == "Johnson-Champoux-Allard":
            self.johnson_koplik_dashen(f)
            self.champoux_allard(f)
        elif self.fluid_model == "Johnson-Lafarge":
            self.johnson_koplik_dashen(f)
            self.lafarge(f)
        elif self.fluid_model == "Pride-Lafarge":
            self.pride(f)
            self.lafarge(f)
        else:
            raise ValueError("Unknown fluid model: {}".format(self.fluid_model))

    def delany_bazley(self, f):
        omega = 2*pi*f
        X = f/self.sigma
        Z = Air.rho*Air.c*( 1 + 9.08*(X*1000)**(-0.75) - 1j*11.9*(X*1000)**(-0.73) ); 
        k = omega/Air.c * (-1j) * ( 10.3*(X*1000)** (-0.59) +1j* ( 1 + 10.8*(X*1000)**(-0.70) ) );
        self.K_eq_til = Z*omega/k
        self.rho_eq_til = k*Z/omega

    def miki(self, f):
        omega = 2*pi*f
        X = f/self.sigma
        Z = Air.rho*Air.c*( 1 + 5.50*(X*1000)**(-0.632)- 1j*8.43*(X*1000)**(-0.632) ); 
        k = omega/Air.c * (-1j) * ( 11.41*(X*1000)**(-0.618)+ 1j* (1 + 7.81*(X*1000)**(-0.618) ) )
        self.K_eq_til = Z*omega/k
        self.rho_eq_til = k*Z/omega

    def johnson_koplik_dashen(self, f):
        # JKD model for rho_eq_til
        omega = 2*pi*f
        omega_0 = self.sigma*self.phi/(Air.rho*self.alpha_inf)
        omega_inf= (self.sigma*self.phi*self.Lambda)**2/(4*Air.mu*Air.rho*self.alpha_inf**2)
        self.F = sqrt(1+1j*omega/omega_inf)
        self.rho_eq_til = (Air.rho*self.alpha_inf/self.phi)*(1+(omega_0/(1j*omega))*self.F)

    def pride(self, f):
        # Pride model for rho_eq_til
        omega = 2*pi*f
        omega_0 = self.sigma*self.phi/(Air.rho*self.alpha_inf)
        omega_inf= (self.sigma*self.phi*self.Lambda)**2/(4*Air.mu*Air.rho*self.alpha_inf**2)
        beta = self.alpha_0/self.alpha_inf-1
        p = omega_0/(2*beta*omega_inf)
        F = 1-p+p*sqrt(1+1j*omega/(omega_inf*p**2))
        self.rho_eq_til = (Air.rho*self.alpha_inf/self.phi)*(1+(omega_0/(1j*omega))*F)

    def champoux_allard(self, f):
        #  Champoux-Allard model for K_eq_til
        omega = 2*pi*f
        self.omega_prime_infty = (16*Air.nu_prime)/(self.Lambda_prime**2)
        self.F_prime_CA = sqrt(1+1j*omega/self.omega_prime_infty)
        self.alpha_prime_til = 1+self.omega_prime_infty*self.F_prime_CA/(2*1j*omega)
        self.K_eq_til = (Air.gamma*Air.P/self.phi)/(Air.gamma-(Air.gamma-1)/self.alpha_prime_til)

    def lafarge(self, f):
        #  Lafarge model for K_eq_til
        omega = 2*pi*f
        self.F_prime_CAL = sqrt(1+4j*omega*self.k0prime**2*Air.C_p*Air.rho/(self.Lambda_prime**2*self.phi**2*Air.lambda_))
        self.alpha_prime_til = 1 - 1j*self.phi*Air.lambda_/(self.k0prime*Air.C_p*Air.rho*omega)*self.F_prime_CAL
        self.K_eq_til = (Air.gamma*Air.P/self.phi)/(Air.gamma-(Air.gamma-1)/self.alpha_prime_til)

    def limp_parameters(self,f):
        self.rho_12 = -self.phi*Air.rho*(self.alpha_inf-1)
        self.rho_11 = self.rho_1-self.rho_12
        self.rho_2 = self.phi*Air.rho
        self.rho_22 = self.rho_2-self.rho_12

        self.rho_22_til = self.phi**2*self.rho_eq_til
        self.rho_12_til = self.rho_2-self.rho_22_til
        self.rho_11_til = self.rho_1-self.rho_12_til
        self.rho_til = self.rho_11_til-((self.rho_12_til**2)/self.rho_22_til)

        self.gamma_til = self.phi*(self.rho_12_til/self.rho_22_til-(1-self.phi)/self.phi)
        self.rho_s_til = self.rho_til+self.gamma_til**2*self.rho_eq_til


        self.rho_eq_til *= self.rho_til/self.rho_s_til


    def biot_parameters(self, f):
        omega = 2*pi*f

        self.lambda_ = self.E*self.nu/((1+self.nu)*(1-2*self.nu))
        self.rho_12 = -self.phi*Air.rho*(self.alpha_inf-1)
        self.rho_11 = self.rho_1-self.rho_12
        self.rho_2 = self.phi*Air.rho
        self.rho_22 = self.rho_2-self.rho_12

        self.rho_22_til = self.phi**2*self.rho_eq_til
        self.rho_12_til = self.rho_2-self.rho_22_til
        self.rho_11_til = self.rho_1-self.rho_12_til
        self.rho_til = self.rho_11_til-((self.rho_12_til**2)/self.rho_22_til)

        self.gamma_til = self.phi*(self.rho_12_til/self.rho_22_til-(1-self.phi)/self.phi)
        self.rho_s_til = self.rho_til+self.gamma_til**2*self.rho_eq_til
        self.structural_loss = 1+1j*self.eta_s

        self.N = self.E/(2*(1+self.nu))*self.structural_loss
        self.A_hat = (self.E*self.nu)/((1+self.nu)*(1-2*self.nu))*self.structural_loss
        self.P_hat = self.A_hat+2*self.N

        # Biot 1956 elastic coefficients
        self.R_til = self.K_eq_til*self.phi**2
        self.Q_til = ((1-self.phi)/self.phi)*self.R_til
        self.P_til = self.P_hat+self.Q_til**2/self.R_til

        delta_eq = omega*sqrt(self.rho_eq_til/self.K_eq_til)
        delta_s_1 = omega*sqrt(self.rho_til/self.P_hat)
        delta_s_2 = omega*sqrt(self.rho_s_til/self.P_hat)

        Psi = ((delta_s_2**2+delta_eq**2)**2-4*delta_eq**2*delta_s_1**2)
        sdelta_total = sqrt(Psi)

        delta_1 = sqrt(0.5*(delta_s_2**2+delta_eq**2+sdelta_total))
        delta_2 = sqrt(0.5*(delta_s_2**2+delta_eq**2-sdelta_total))
        delta_3 = omega*sqrt(self.rho_til/self.N)

        mu_1 = self.gamma_til*delta_eq**2/(delta_1**2-delta_eq**2)
        mu_2 = self.gamma_til*delta_eq**2/(delta_2**2-delta_eq**2)
        mu_3 = -self.gamma_til

        self.delta_s_1 = delta_s_1
        self.delta_s_2 = delta_s_2
        self.delta_1 = delta_1
        self.delta_2 = delta_2
        self.delta_3 = delta_3
        self.delta_eq = delta_eq
        self.mu_1 = mu_1
        self.mu_2 = mu_2
        self.mu_3 = mu_3




