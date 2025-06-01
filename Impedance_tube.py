#  ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <olivier.dazel@univ-lemans.fr> and <alan.geslain@univ-lemans.fr> wrote this file.  
#  As long as you retain this notice you can do whatever you want with this stuff. 
#  If we meet some day, and you think this stuff is worth it, you can buy us a beer in return.       
# Olivier DAZEL and Alan GESLAIN
#  ----------------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt

from numpy import pi, sqrt, cos, sin
from material import Air, Material

possible_output_types = ["Absorption", "Surface Impedance", "R"]

class Impedance_tube(Material):
    def __init__(self, params):
        super().__init__(params)
        self.termination_condition = params['termination_condition']

    def update_from_params(self, params):
        self.d = params['Thickness'].value/ 1e2
        self.termination_condition = "rigid"


    def compute_surface_impedance_rigid_backing(self, h, f):
        omega = 2*pi*f
        if self.solid_model in ["Rigid frame", "Limp"]:
            k = omega*sqrt(self.rho_eq_til/self.K_eq_til)
            Z = sqrt(self.rho_eq_til*self.K_eq_til)
            Z_s = -1j*Z/np.tan(k*h)
        else:
            N = self.K_eq_til*self.delta_1*self.delta_2*(self.mu_2-self.mu_1)
            D = 1j*omega*(self.delta_1*self.mu_2*np.tan(self.delta_2*h)-self.delta_2*self.mu_1*np.tan(self.delta_1*h))
            Z_s= N/D
        return Z_s/(Air.Z)

    def compute_indicators(self , h, f):
        if self.termination_condition == "rigid":
            Z_s= self.compute_surface_impedance_rigid_backing(h, f)
            R = (Z_s-1)/(Z_s+1)
            T = 0*R
        else:
            jkljkl
        alpha = 1 - np.abs(R) **2
        d ={'frequency': f, 'Absorption': alpha, "R": R, "Surface Impedance": Z_s}
        return d

    def run(self,f, params):
        Material.update_from_params(self, params)
        self.update_from_params(params)
        Material.equivalent_fluid_parameters(self, f)
        if self.solid_model == "Limp":
            self.limp_parameters(f)
        elif self.solid_model == "Biot":
            self.biot_parameters(f) 



        return self.compute_indicators(self.d, f)


            

