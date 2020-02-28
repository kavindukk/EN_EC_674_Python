"""
Class to determine wind velocity at any given moment,
calculates a steady wind speed and uses a stochastic
process to represent wind gusts. (Follows section 4.4 in uav book)
"""
import sys
sys.path.append('..')
from tools.transfer_function import transfer_function
from tools.tools import parameter_calc as pc
import numpy as np

import math as m
root = m.sqrt

class wind_simulation:
    def __init__(self, Ts):
        # steady state wind defined in the inertial frame

        # a1, b1, a2, a3, b2, a4, a5, b3 = pc()
        Lu, Lv, Lw =200, 200, 50
        sigma_u, sigma_v  = 1.06, 1.06
        sigma_w = .7
        V_a = 20

        c_u = sigma_u*root(2*V_a/Lu/np.pi)
        c_v = sigma_v*root(3*V_a/Lv/np.pi)
        c_w = sigma_w*root(3*V_a/Lw/np.pi)

        a1 = c_u
        b1 = V_a/Lu
        a2 = c_v
        a3 = c_v*V_a/Lv/root(3)
        b2 = V_a/Lv
        a4 = c_w
        a5 = c_w*V_a/Lw/root(3)
        b3 = V_a/Lw
        self._steady_state = np. array([ [0,0,0]])

        self.u_w = transfer_function(num=np.array([[a1]]),
                                     den=np.array([[1, b1]]),
                                     Ts=Ts)
        self.v_w = transfer_function(num=np.array([[a2, a3]]),
                                     den=np.array([[1, 2*b2, b2**2.0]]),
                                     Ts=Ts)
        self.w_w = transfer_function(num=np.array([[a4, a5]]),
                                     den=np.array([[1, 2*b3, b3**2.0]]),
                                     Ts=Ts)
        self._Ts = Ts

    def update(self):
        # returns a six vector.
        #   The first three elements are the steady state wind in the inertial frame
        #   The second three elements are the gust in the body frame
        gust = np.array([[self.u_w.update(np.random.randn())],
                         [self.v_w.update(np.random.randn())],
                         [self.w_w.update(np.random.randn())]]).T
        # gust = np.array([[0.],[0.],[0.]]).T
        return np.concatenate(( self._steady_state, gust ))

