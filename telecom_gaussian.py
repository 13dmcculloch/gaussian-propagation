# -*- coding: utf-8 -*-
"""
Gaussian beam propagator

Created on Thu May 30 16:22:33 2024

@author: drm1g20
"""

import math
import numpy as np
import matplotlib.pyplot as plt

import json

class Beam:
    def __init__(self, filename='', 
                 wavelength=800e-9, waist_init=4.6876E-6/2, R_init=10e10):
        if filename:
            self.load_parameters(filename)
        else:
            self.wavelength = wavelength
            self.w_init = waist_init
            self.R_init = R_init
            
            self.z_max = 1
            self.sim_res = int(1e5)
            
        self.q_init = self.q_scalar()
        self.q = False  # haven't evaluated q yet
        
        self.z_R = math.pi * self.w_init**2 / self.wavelength
        
        self.z, self.z_inc = np.linspace(0, self.z_max, self.sim_res, 
                                         retstep=True)
        
    
    def print_details(self):
        print("Beam parameters")
        print("\tWavelength:", self.wavelength * 1e9, "nm")
        print("\tInitial half waist:", self.w_init * 1e6, "um")
        print("\tInitial bend radius:", self.R_init, "m\n")
        
        print("Simulation parameters")
        print("\tLength:", self.z_max, "m")
        print("\tSteps:", self.sim_res / 1e6, "M")
        print("\tResolution:", self.z_max * 1e6 / self.sim_res, "um\n")
        
        print("Lenses")
        for alias in self.lenses.keys():
            lens = self.lenses[alias]
            print("\tLens", alias)
            print("\t\tFocal length:", lens['focal_length'] * 1e3, "mm")
            print("\t\tz position:", lens['z_position'], "m\n")
        
        
    def load_parameters(self, filename):
        """
        Load beam initial conditions from json file, 1 deep, and assign
        """
        with open(filename, 'r') as f:
            self.parameters = json.load(f)
        
        self.wavelength = self.parameters.get('wavelength', 1550e-9)
        self.w_init = self.parameters.get('initial_waist', 5e-6)
        self.R_init = self.parameters.get('initial_curvature', 10e10)
        
        self.z_max = self.parameters.get('simulation_length', 1)
        self.sim_res = int(self.parameters.get('simulation_resolution', 
                                                  1e5))
        
        
    def load_lenses(self, filename):
        """
        Load lenses from json file with format:
            "alias" {
                "z_position": [z (m)],
                "focal_length": [f (m)]
            }
        """
        self.lenses = {}
        
        with open(filename, 'r') as f:
            self.lenses = json.load(f) 
            
        # check lens position is in z
        for alias in self.lenses.keys():
            lens = self.lenses[alias]
            if (lens.get('z_position', 0) > np.max(self.z) or
                lens.get('z_position', 0) < np.min(self.z)):
                print("Lens", alias, "is outside of simulation region!")
               
                
    """
    def q_fast(self):
        q = np.zeros(len(self.z), dtype=np.complex_)
        q[0] = 1/self.q_init
        
        float_tolerance = self.z_inc  # I need to tweak
        done_lenses = []
        
        z_fs = 0
        z_last = 0
        
        print("Fitting lenses...", end='')
                
        for k in self.lenses.keys():
            z_pos_lens = float(self.lenses[k]['z_position'])
            
            for idx in range(1, len(self.z)):
                z_pos = self.z[idx]
                
                if (z_pos > z_pos_lens - float_tolerance and
                    z_pos < z_pos_lens + float_tolerance):
                    
                    if k in done_lenses:
                        break
                    
                    # apply free space propagation to q between lenses
                    z_fs = z_pos - z_last
                    s = Space(z_fs, z_pos)
                   
                    f = self.lenses[k]['focal_length']
                    
                    z_fs = 0
                    z_last = z_pos
        
        pass
    """
      
    
    def q_evaluate(self):
        """
        Returns 1/q ndarray and also stores it in object as self.q
        """
        q = [1/self.q_init]  # z == 0
        
        float_tolerance = self.z_inc  # I need to tweak
        done_lenses = []
        done_lens_info = []
        
        stop_idx = len(self.z)
        start_idx = 1
        print("Evaluating...", end='')
        for idx in range(start_idx, stop_idx):
            z_pos = self.z[idx]

            q_p = q[idx - 1]
            
            at_lens = False
            for k in self.lenses.keys():
                z_pos_lens = float(self.lenses[k]['z_position'])
                # evaluate q_n at the lens if deemed to be at one
                if (z_pos > z_pos_lens - float_tolerance and 
                           z_pos < z_pos_lens + float_tolerance):
                    
                    if k in done_lenses:
                        break
                    
                    at_lens = True
                    
                    f = self.lenses[k]['focal_length']
                    
                    lens = Lens(f, z_pos_lens)
                    q_lens = self.Mq(q_p, lens.matrix)
                    
                    done_lenses.append(k)
                    done_lens_info.append((k, z_pos, f))
                    
                    break
                    
            
            s = Space(self.z_inc, z_pos)
            if at_lens:
                q_n = self.Mq(q_lens, s.matrix)
            else:
                q_n = self.Mq(q_p, s.matrix)
                
            q.append(q_n)
            
            
            if idx % (len(self.z) / 20) == 0:
                pass
            
        print("Done")
        
        for alias, z_pos, f in done_lens_info:
            print("Fitted lens", alias, 
                  "at", z_pos, 
                  "m with f =", f * 1e3,
                  "mm")
        
        print("")
        
        self.q = 1/np.array(q)
        
        return self.q
            
    
    def focal_points(self):
        if not self.q.any():
            print("Need to evaluate q first!\nobj.q_evaluate()")
            return
        
        w = np.sqrt((self.wavelength / math.pi) * -1/np.imag(self.q))
        w_p = np.gradient(w)
        
        print("Finding focal points...", end='')
        focal_points = []
        z_prev = self.z[0]
        for idx in range(1, len(w_p)):
            z_pos = self.z[idx]
            if w_p[idx] > 0 and w_p[idx - 1] < 0:
                w_loc = z_pos
                w_val = w[idx]
                
                focal_points.append((w_loc, w_val))
           
        print("Done")
        
        for w_loc, w_val in focal_points:
            print("At z=", w_loc, "m, w=", w_val, "um")        
    
    
    def plot_both_norm(self):
        self.fig, self.ax = plot_setup("Radius of curvature and beam waist",
                                       "z/z_R (m)",
                                       "Norm (m)")
        
        R = 1/np.real(self.q)
        self.ax.plot(self.z/self.z_R, R/self.z_R, label='R(z)/z_R')
        
        
        w = np.sqrt((self.wavelength / math.pi) * -1/np.imag(self.q))
        self.ax.plot(self.z/self.z_R, w/self.w_init, label='w(z)/w(0)')
        
        self.ax.set_ylim([0, 10])
        
        self.ax.legend()


    def plot_R(self):
        self.fig, self.ax = plot_setup("Gaussian Beam 2",
                                       "z (m)",
                                       "Radius of curvature (m)")
        R = np.real(self.q)
        
        self.ax.plot(self.z, R, label='R(z)')
        
    
    def plot_q(self):
        self.fig, self.ax = plot_setup("Gaussian Beam 2",
                                       "z (m)",
                                       "Parameter (m)")
        self.ax.plot(self.z, np.real(self.q), label='Re(q(z))')
        self.ax.plot(self.z, np.imag(self.q), label='Im(q(z))')
        
        self.ax.legend()
        
    
    def plot_w(self):
        w = np.sqrt((self.wavelength / math.pi) * -1/np.imag(self.q))
        
        self.fig, self.ax = plot_setup("Gaussian Beam 2",
                                       "z (m)",
                                       "Beam waist (mm)")
        
        self.ax.plot(self.z, w * 1e3, label='w(z)')
        
        self.ax.legend()
        
    
    def q_scalar(self):
        """
        Get first 1/q value from intials
        """
        return (self.R_init**-1 - 1j * (self.wavelength / 
                                    (np.pi * self.w_init**2)))
    
    
    def add_lenses(self, *components):
        pass
    
        
    def Mq(self, q_in, M):
        """
        Returns q_out from M.q_in
        q_in - current q
        M - matrix representation of optical component
        """
        num = M[0][0] * q_in + M[0][1]
        den = M[1][0] * q_in + M[1][1]
        
        return num / den
    

def plot_generic(x, y):
    fig, ax = plot_setup("Gaussian beam derivatives",
                        "z",
                        "beam")
    
    ax.plot(x,y)
    
    ax.set_xlim([0.899, 0.901])
    ax.set_ylim([-0.8e-6, 0.8e-6])
    
        
class opticalComponent():
    """
    Contains matrix setters for ABCD representation of component
    """
    def __init__(self, a, b, c, d, z):
        """
        Constructor returns matrix
        """
        self.matrix = np.array([[a, b], [c, d]])
        self.z_pos = z
    
    
class Lens(opticalComponent):
    """
    Constructor
    """
    def __init__(self, f, z_pos):
        super().__init__(1, 0, -1/f, 1, z_pos)
    
    
class Space(opticalComponent):
    """
    Constructor 
    """
    def __init__(self, d, z_pos):
        super().__init__(1, d, 0, 1, z_pos)

    
def plot_setup(title, xlabel, ylabel):
    plt.rcParams.update({'font.size': 13})
    fig, ax = plt.subplots(dpi=600)
        
    ax.minorticks_on()
    ax.tick_params(direction='in', which='both')
    ax.grid(visible=True, axis='both', which='major', linestyle='--', 
                alpha=0.5)
    
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    return fig, ax


def loading_bar(norm, length):
    """
    norm = pos/end: float
    """
    # clear line
    print('\r', end='')
    for i in range(length + len("Loading" + ' xx% []')):
        print(' ', end='')

    n = math.floor(norm * length)
    if norm > 0.98:
        n = length
        
    print("\rLoading", end=' ')
    print(str(n * 100 / length) + '% [', end ='')
	
    for i in range(n):
        print('=', end='')
    
    print(']', end='')
    
    
beam = Beam('beam.json')
beam.load_lenses("lenses.json")
beam.print_details()

beam.q_evaluate()

beam.plot_w()
beam.focal_points()