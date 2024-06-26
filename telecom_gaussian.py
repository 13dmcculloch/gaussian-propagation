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
            
        self.sim_resolution = self.z_max / self.sim_res
        if self.sim_resolution > 1e-6:
            print("[WARNING]: Simulation resolution exceeds 1 um"
                  " and will lead to significant linear diffraction effects.")
            
        self.q_init = self.q_scalar()
        self.q = False  # haven't evaluated q yet
        self.w = False
        self.R = False
        
        self.z_R = math.pi * self.w_init**2 / self.wavelength
        
        self.z, self.z_inc = np.linspace(0, self.z_max, self.sim_res, 
                                         retstep=True)
        
        self.focals = [(np.min(self.z), self.w_init)]
        
    
    def print_details(self):
        print("Beam parameters")
        print("\tWavelength:", self.wavelength * 1e9, "nm")
        print("\tInitial half waist:", self.w_init * 1e6, "um")
        print("\tInitial bend radius:", self.R_init, "m\n")
        
        print("Simulation parameters")
        print("\tLength:", self.z_max, "m")
        print("\tSteps:", self.sim_res / 1e6, "M")
        print("\tResolution:", self.sim_resolution * 1e6, "um\n")
        
        print("Lenses")
        for alias in self.lenses.keys():
            lens = self.lenses[alias]
            print("\tLens", alias)
            print("\t\tFocal length:", lens['focal_length'] * 1e3, "mm")
            print("\t\tz position:", lens['z_position'], "m\n")
        
    
    def print_focals(self):
        i = 0
        for w_loc, w_val in self.focals:
            print(i, "At z=", w_loc, "m, w=", w_val * 1e6, "um")
            i += 1
            
        print("")
        
        w_prev = self.focals[0][1]
        for idx in range(1, len(self.focals)):
            w_cur = self.focals[idx][1]
            
            print("Ratio between", idx, "and", idx - 1, "=", 
                  w_cur / w_prev)
            
            w_prev = w_cur
            
    
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
                    
                    done_lenses.append(k)
        
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
                    radius = self.lenses[k]['radius']
                    
                    lens = Lens(f, z_pos_lens, radius)
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
        
        self.done_lens_info = done_lens_info
        self.q = 1/np.array(q)
        
        print("Calculating beam waist...", end='')
        self.w = np.sqrt((self.wavelength / math.pi) * -1/np.imag(self.q))
        print("Done")
        
        print("Calculating bend radius...", end='')
        self.R = np.real(self.q)
        print("Done\n")
        
        return self.q
            
    
    def print_q_missing(self):
        print("[ERROR]: Need to evaluate q first!")
        
    
    def focal_points(self):
        if type(self.q) == bool or type(self.w) == bool:
            self.print_q_missing()
            return
        
        w_p = np.gradient(self.w)
        
        print("Finding focal points...", end='')
        for idx in range(1, len(w_p)):
            z_pos = self.z[idx]
            if w_p[idx] > 0 and w_p[idx - 1] < 0:
                w_loc = z_pos
                w_val = self.w[idx]
                
                self.focals.append((w_loc, w_val))
           
        print("Done")
        
        self.print_focals()    
    
    
    def waist_radius_check(self):
        if type(self.w) == bool:
            self.print_q_missing()
            return 
        
        for alias, z_pos, f in self.done_lens_info:
            radius = self.lenses[alias]['radius']
            
            idx = int(z_pos / self.z_inc)
            waist = self.w[idx]
            
            if waist > radius:
                print("Beam waist is outside OD of lens", alias)
    
    
    def plot_lens_radius(self):
        c = 0.0
        for k in self.lenses.keys():
            lens = self.lenses[k]
            self.ax.vlines(x=lens['z_position'], 
                            ymin=0, 
                            ymax=lens['radius'] * 1e3,
                            color=str(c),
                            label='lens ' + k)
            
            c += 0.2
            
        self.ax.set_ylim([0, lens['radius'] * 1e3 + 3])
    
    
    def plot_both_norm(self):
        self.fig, self.ax = plot_setup("Radius of curvature and beam waist",
                                       "z/z_R (m)",
                                       "Norm (m)")
        
        self.ax.plot(self.z/self.z_R, self.R/self.z_R, label='R(z)/z_R')
        
        self.ax.plot(self.z/self.z_R, self.w/self.w_init, label='w(z)/w(0)')
        
        self.ax.set_ylim([0, 10])
        
        self.ax.legend()


    def plot_R(self):
        self.fig, self.ax = plot_setup("Gaussian Beam 2",
                                       "z (m)",
                                       "Radius of curvature (m)")
        
        self.ax.plot(self.z, self.R, label='R(z)')
        
    
    def plot_q(self):
        self.fig, self.ax = plot_setup("Gaussian Beam 2",
                                       "z (m)",
                                       "Parameter (m)")
        self.ax.plot(self.z, np.real(self.q), label='Re(q(z))')
        self.ax.plot(self.z, np.imag(self.q), label='Im(q(z))')
        
        self.ax.legend()
        
    
    def plot_w(self):
        self.fig, self.ax = plot_setup("Gaussian Beam 2",
                                       "z (m)",
                                       "Beam waist (mm)")
        
        self.plot_lens_radius()
        
        self.ax.plot(self.z, self.w * 1e3, label='w(z)')
        
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
    def __init__(self, f, z_pos, radius):
        super().__init__(1, 0, -1/f, 1, z_pos)
        self.radius = radius
        
    
    def check_radius(self, radius):
        """
        Returns false if argument is less than lens radius
        """
        return radius < self.radius
    
    
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
beam.plot_lens_radius()

beam.focal_points()
beam.waist_radius_check()
