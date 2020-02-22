"""
pid_control
    - Beard & McLain, PUP, 2012
    - Last Update:
        2/6/2019 - RWB
"""
import sys
import numpy as np
sys.path.append('..')

class pid_control:
    def __init__(self, kp=0.0, ki=0.0, kd=0.0, Ts=0.01, sigma=0.05, limit=1.0):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.Ts = Ts
        self.limit = limit
        self.integrator = 0.0
        self.error_delay_1 = 0.0
        self.error_dot_delay_1 = 0.0
        # gains for differentiator
        self.a1 = (2.0 * sigma - Ts) / (2.0 * sigma + Ts)
        self.a2 = 2.0 / (2.0 * sigma + Ts)

    def update(self, y_ref, y, reset_flag=False):
        if reset_flag == True:
            self.integrator = 0.0
            self.error_delay_1 = 0.0
            self.error_dot_delay_1 =0.0
        #compute the error
        error = y_ref - y
        #update the integrator with trapazoidal rule
        self.integrator = self.integrator + (self.Ts/2)*(error + self.error_delay_1)
        #update the differentiator
        error_dot = self.a1*self.error_dot_delay_1 + self.a2*(error - self.error_delay_1)
        #PID control
        u = self.kp*error + self.ki*self.integrator + self.kd*error_dot
        #saturate PID control at limit
        u_sat = self._saturate(u)
        #integral anti-windup
        #adjust integrator to keep u out of saturation
        if np.abs(self.ki)>0.0001:
            self.integrator = self.integrator + (self.Ts/self.ki)*(u_sat - u)
        #update the delayed variables
        self.error_delay_1 = error
        self.error_dot_delay_1 = error_dot
        return u_sat

    def update_with_rate(self, y_ref, y, ydot, reset_flag=False):
        if reset_flag == True:
            self.integrator = 0.0
            self.error_delay_1 = 0.0
        #compute the error
        error = y_ref - y
        #update the integrator with trapazoidal rule
        self.integrator = self.integrator + (self.Ts/2)*(error + self.error_delay_1)
        #PID control
        u = self.kp*error + self.ki*self.integrator - self.kd*ydot
        #saturate PID control at limit
        u_sat = self._saturate(u)
        #integral anti-windup
        #adjust integrator to keep u out of saturation
        if np.abs(self.ki)>0.0001:
            self.integrator = self.integrator + (self.Ts/self.ki)*(u_sat - u)
        #update the delayed variables
        self.error_delay_1 = error
        return u_sat

    def _saturate(self, u):
        # saturate u at +- self.limit
        if u >= self.limit:
            u_sat = self.limit
        elif u <= -self.limit:
            u_sat = -self.limit
        else:
            u_sat = u
        return u_sat

class pi_control:
    def __init__(self, kp=0.0, ki=0.0, Ts=0.01, limit=1.0):
        self.kp = kp
        self.ki = ki
        self.Ts = Ts
        self.limit = limit
        self.integrator = 0.0
        self.error_delay_1 = 0.0

    def update(self, y_ref, y):
        #compute the error
        error = y_ref - y
        #update the integrator with trapazoidal rule
        self.integrator = self.integrator + (self.Ts/2)*(error + self.error_delay_1)
        #PI control
        u = self.kp*error + self.ki*self.integrator
        #saturate PID control at limit
        u_sat = self._saturate(u)
        #integral anti-windup
        #adjust integrator to keep u out of saturation
        if np.abs(self.ki)>0.0001:
            self.integrator = self.integrator + (self.Ts/self.ki)*(u_sat - u)
            #update the delayed variables
        self.error_delay_1 = error
        # self.error_dot_delay_1 = error_dot
        return u_sat

    def _saturate(self, u):
        # saturate u at +- self.limit
        if u >= self.limit:
            u_sat = self.limit
        elif u <= -self.limit:
            u_sat = -self.limit
        else:
            u_sat = u
        return u_sat

class pd_control_with_rate:
    # PD control with rate information
    # u = kp*(yref-y) - kd*ydot
    def __init__(self, kp=0.0, kd=0.0, limit=1.0):
        self.kp = kp
        self.kd = kd
        self.limit = limit

    def update(self, y_ref, y, ydot):
        #compute the error
        error = y_ref - y
        #PD control
        u = self.kp*error - self.kd*ydot
        #saturate PD control at limit
        u_sat = self._saturate(u)
        return u_sat

    def _saturate(self, u):
        # saturate u at +- self.limit
        if u >= self.limit:
            u_sat = self.limit
        elif u <= -self.limit:
            u_sat = -self.limit
        else:
            u_sat = u
        return u_sat