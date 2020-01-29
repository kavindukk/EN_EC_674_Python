"""
mav_dynamics
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state

chap_edit_branch
"""
import sys
sys.path.append('..')
import numpy as np
import math as m

# load message types
from message_types.msg_state import msg_state

from chap3.equations import equations

import parameters.aerosonde_parameters as MAV
from tools.tools import Quaternion2Euler, Inertial2Body, Body2Inertia

# Quaternion2Rotation, 

class mav_dynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.pn0],  # (0)
                               [MAV.pe0],   # (1)
                               [MAV.pd0],   # (2)
                               [MAV.u0],    # (3)
                               [MAV.v0],    # (4)
                               [MAV.w0],    # (5)
                               [MAV.e0],    # (6)
                               [MAV.e1],    # (7)
                               [MAV.e2],    # (8)
                               [MAV.e3],    # (9)
                               [MAV.p0],    # (10)
                               [MAV.q0],    # (11)
                               [MAV.r0]])   # (12)
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        self._update_velocity_data()
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0
        self._beta = 0
        # initialize true_state message
        self.msg_true_state = msg_state()

    ###################################
    # public functions
    def update_state(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)

        # Integrate ODE using Runge-Kutta RK4 algorithm
        time_step = self._ts_simulation
        k1 = self._derivatives(self._state, forces_moments)
        k2 = self._derivatives(self._state + time_step/2.*k1, forces_moments)
        k3 = self._derivatives(self._state + time_step/2.*k2, forces_moments)
        k4 = self._derivatives(self._state + time_step*k3, forces_moments)
        self._state += time_step/6 * (k1 + 2*k2 + 2*k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = np.sqrt(e0**2+e1**2+e2**2+e3**2)
        self._state[6][0] = self._state.item(6)/normE
        self._state[7][0] = self._state.item(7)/normE
        self._state[8][0] = self._state.item(8)/normE
        self._state[9][0] = self._state.item(9)/normE

        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)

        # update the message class for the true state
        self._update_msg_true_state()

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):   

        pn = self._state[0][0]
        pe = self._state[1][0]
        pd = self._state[2][0]
        u = self._state[3][0]
        v = self._state[4][0]
        w = self._state[5][0]
        e0 = self._state[6][0]
        e1 = self._state[7][0]
        e2 = self._state[8][0]
        e3 = self._state[9][0]
        p = self._state[10][0]
        q = self._state[11][0]
        r = self._state[12][0]

        # extract forces and moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        l = forces_moments.item(3)
        m = forces_moments.item(4)
        n = forces_moments.item(5)

        eq = equations(u, v, w, MAV.Jx, MAV.Jy, MAV.Jz, MAV.Jxz, p, q, r, l, m, n)

        # position kinematics
        pn_dot = (e1**2+e0**2-e2**2-e3**2)*u +2*(e1*e2-e3*e0)*v +2*(e1*e3+e2*e0)*w
        pe_dot = 2*(e1*e2 +e3*e0)*u + (e2**2+ e0**2 -e1**2-e3**2)*v +2*(e2*e3 - e1*e0)*w        
        pd_dot = -2*(e1*e3 -e2*e0)*u -2*(e2*e3 +e1*e0)*v - (e3**2+e0**2-e1**2-e2**2)

        # position dynamics
        u_dot = r*v -q*w + fx/MAV.mass
        v_dot = p*w - r*u + fy/MAV.mass
        w_dot = q*u - p*v + fz/MAV.mass
  
        e0_dot = 0.5*(-p*e1 - q*e2 - r*e3)
        e1_dot = 0.5*(p*e0 + r*e2 - q*e3)
        e2_dot = 0.5*(q*e0 - r*e1 + +p*e3)
        e3_dot = 0.5*(r*e0 + q*e1 - p*e2 )
        
        # rotatonal dynamics
        p_dot = eq.rot_dyn()[0]
        q_dot = eq.rot_dyn()[1]
        r_dot = eq.rot_dyn()[2]

        # collect the derivative of the states
        x_dot = np.array([[pn_dot, pe_dot, pd_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T      
         
       
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6,1))):
        # compute airspeed
        phi, theta, psi = Quaternion2Euler(self._state[6:10])

        Wnb, Web, Wdb = Inertial2Body(phi, theta, psi, wind.item(0), wind.item(1), wind.item(2))

        u, v, w = self._state[3], self._state[4], self._state[5]

        u_w, v_w, w_w = Wnb+wind.item(3), Web+wind.item(4), Wdb+wind.item(5) 

        u_r, v_r, w_r = u - u_w, v - v_w, w - w_w

        # print(u)
        # print(v)
        # print(w)

        self._Va = np.sqrt(u_r**2 + v_r**2 + w_r**2)
        # compute angle of attack
        self._alpha = m.atan(w_r/u_r)
        # compute sideslip angle
        self._beta = m.asin(v_r/self._Va)

    def _forces_moments(self, delta):
        """        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        phi, theta, psi = Quaternion2Euler([self._state[6], self._state[7], self._state[8], self._state[9]])

        C_X_alpha = -MAV.C_D_alpha*m.cos(self._alpha) + MAV.C_L_alpha*m.sin(self._alpha)
        C_X_q_alpha = -MAV.C_D_q*m.cos(self._alpha) + MAV.C_L_q*m.sin(self._alpha)
        C_X_delta_e_alpha = -MAV.C_D_delta_e*m.cos(self._alpha) + MAV.C_L_delta_e*m.sin(self._alpha)
        C_Z_alpha = -MAV.C_D_alpha*m.sin(self._alpha) - MAV.C_L_alpha*m.cos(self._alpha)
        C_Z_q_alpha = -MAV.C_D_q*m.sin(self._alpha) - MAV.C_L_q*m.cos(self._alpha)
        C_Z_delta_e_alpha = -MAV.C_D_delta_e*m.sin(self._alpha) - MAV.C_L_delta_e*m.cos(self._alpha)

        rho, Va, S, q, mass, g, c = MAV.rho, self._Va, MAV.S_wing, self._state[11], MAV.mass, MAV.gravity, MAV.c

        fx = -mass*g*m.sin(theta) + .5*rho*Va**2*S*(C_X_alpha + C_X_q_alpha*c*q/2/Va  + C_X_delta_e_alpha*delta[1] )

        fy = mass*g*m.cos(theta)*m.sin(phi) + .5*rho*Va**2*S*(MAV.C_Y_0 + MAV.C_Y_beta*self._beta + MAV.C_Y_p*MAV.b*self._state[12]/2/Va + MAV.C_Y_delta_a*delta[0] + MAV.C_Y_delta_r*delta[2] )

        fz = mass*g*m.cos(theta)*m.cos(phi) + .5*rho*Va**2*S*(C_Z_alpha + C_Z_q_alpha*c*self._state[11]/2/Va + C_Z_delta_e_alpha*delta[1])

        C_Y_0 = MAV.C_Y_0
        C_Y_beta = MAV.C_Y_beta
        C_Y_p = MAV.C_Y_p
        C_Y_r = MAV.C_Y_r
        C_Y_delta_a = MAV.C_Y_delta_a
        C_Y_delta_r = MAV.C_Y_delta_r
        C_ell_0 = MAV.C_ell_0
        C_ell_beta = MAV.C_ell_beta
        C_ell_p = MAV.C_ell_p
        C_ell_r = MAV.C_ell_r
        C_ell_delta_a = MAV.C_ell_delta_a
        C_ell_delta_r = MAV.C_ell_delta_r
        C_n_0 = MAV.C_n_0
        C_n_beta = MAV.C_n_beta
        C_n_p = MAV.C_n_p
        C_n_r = MAV.C_n_r
        C_n_delta_a = MAV.C_n_delta_a
        C_n_delta_r = -MAV.C_n_delta_r
        C_m_0 = MAV.C_m_0
        C_m_alpha = MAV.C_m_alpha
        C_m_q = MAV.C_m_q
        C_m_delta_e = MAV.C_m_delta_e
        b= MAV.b
        c= MAV.c
        p = self._state[10]
        q= self._state[11]
        r = self._state[12]

        Mx = .5*rho*Va**2*S*b*(C_ell_0 + C_ell_beta*self._beta + C_ell_p*b*p/2/Va + C_ell_r*b*r/2/Va + C_ell_delta_a*delta[0] + C_n_delta_r*delta[2])

        My = .5*rho*Va**2*S*c*(C_m_0 + C_m_alpha*self._alpha + C_m_q*c*q/2/Va +C_m_delta_e*delta[1])

        Mz = .5*rho*Va**2*S*b*(C_n_0 + C_n_beta*self._beta + C_n_p*b*p/2/Va + C_n_r*b*r/2/Va + C_n_delta_a*delta[0] + C_n_delta_r*delta[2] )

        V_in = MAV.V_max*delta[3]
        a = MAV.C_Q0*rho*np.power(MAV.D_prop, 5)
        b = (MAV.C_Q1*rho*np.power(MAV.D_prop, 4)/(2*np.pi))*Va + MAV.KQ**2/MAV.R_motor
        c = MAV.C_Q2*rho*np.power(MAV.D_prop, 3)*Va**2 - (MAV.KQ/MAV.R_motor)*V_in + MAV.KQ*MAV.i0
        Omega_op = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
        J_op = 2*np.pi*self._Va/(Omega_op*MAV.D_prop)
        C_T = MAV.C_T2*J_op**2 + MAV.C_T1*J_op + MAV.C_T0
        C_Q = MAV.C_Q2*J_op**2 + MAV.C_Q1*J_op + MAV.C_Q0
        n = Omega_op/(2*np.pi)
        fx += rho*n**2*np.power(MAV.D_prop, 4)*C_T
        # Mx += -rho*n**2*np.power(MAV.D_prop, 5)*C_Q

        self._forces[0] = fx
        self._forces[1] = fy
        self._forces[2] = fz
        return np.array([[fx, fy, fz, Mx, My, Mz]]).T

    def _update_msg_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        self.msg_true_state.pn = self._state.item(0)
        self.msg_true_state.pe = self._state.item(1)
        self.msg_true_state.h = -self._state.item(2)
        self.msg_true_state.Va = self._Va
        self.msg_true_state.alpha = self._alpha
        self.msg_true_state.beta = self._beta
        self.msg_true_state.phi = phi
        self.msg_true_state.theta = theta
        self.msg_true_state.psi = psi
        self.msg_true_state.Vg = self._Va + self._wind 
        vg_i, vg_j, vg_k =  Body2Inertia(phi, theta, psi, self.msg_true_state.Vg )

        self.msg_true_state.gamma = m.asin( -vg_k/m.sqrt(vg_i**2+vg_j**2 + vg_k**2))
        self.msg_true_state.chi = m.atan(vg_j/vg_i)
        self.msg_true_state.p = self._state.item(10)
        self.msg_true_state.q = self._state.item(11)
        self.msg_true_state.r = self._state.item(12)
        self.msg_true_state.wn = self._wind.item(0)
        self.msg_true_state.we = self._wind.item(1)
