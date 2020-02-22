"""
mav_dynamics
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state

chap_edit_branch
"""
import sys
sys.path.append('..')
import numpy as np
import math 

# load message types
from message_types.msg_state import msg_state
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

        pn = self._state.item(0)
        pe = self._state.item(1)
        pd = self._state.item(2)
        u = self._state.item(3)
        v = self._state.item(4)
        w = self._state.item(5)
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)

        # extract forces and moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        l = forces_moments.item(3)
        m = forces_moments.item(4)
        n = forces_moments.item(5)

        # position kinematics
        pn_dot = (e1**2+e0**2-e2**2-e3**2)*u +2*(e1*e2-e3*e0)*v +2*(e1*e3+e2*e0)*w
        pe_dot = 2*(e1*e2 +e3*e0)*u + (e2**2+ e0**2 -e1**2-e3**2)*v +2*(e2*e3 - e1*e0)*w        
        pd_dot = 2*(e1*e3 -e2*e0)*u +2*(e2*e3 +e1*e0)*v +(e3**2+e0**2-e1**2-e2**2)*w

        # position dynamics
        u_dot = r*v -q*w + fx/MAV.mass
        v_dot = p*w - r*u + fy/MAV.mass
        w_dot = q*u - p*v + fz/MAV.mass
  
        e0_dot = 0.5*(-p*e1 - q*e2 - r*e3)
        e1_dot = 0.5*(p*e0 + r*e2 - q*e3)
        e2_dot = 0.5*(q*e0 - r*e1 + +p*e3)
        e3_dot = 0.5*(r*e0 + q*e1 - p*e2 )
        
        # rotatonal dynamics   
        p_dot = MAV.gamma1*p*q - MAV.gamma2*q*r + l
        q_dot = MAV.gamma5*p*r - MAV.gamma6*(p**2 - r**2) + m
        r_dot = MAV.gamma7*p*q - MAV.gamma1*q*r + n

        # p_dot = MAV.gamma1*p*q -MAV.gamma2*q*r + MAV.gamma3*l + MAV.gamma4*n
        # q_dot = MAV.gamma5*p*r -MAV.gamma6*(p**2 - r**2) + m/MAV.Jy
        # r_dot = MAV.gamma7*p*q -MAV.gamma1*q*r + MAV.gamma4*l + MAV.gamma8*n

        # collect the derivative of the states
        x_dot = np.array([[pn_dot, pe_dot, pd_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T      
         
        
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6,1))):

        # compute airspeed        
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
       
        #Wnb = Inertial2Body(phi, theta, psi, wind.item(0), wind.item(1), wind.item(2))[0]
        #Web = Inertial2Body(phi, theta, psi, wind.item(0), wind.item(1), wind.item(2))[1]
        #Wdb = Inertial2Body(phi, theta, psi, wind.item(0), wind.item(1), wind.item(2))[2]

        u = self._state[3]
        v = self._state[4] 
        w = self._state[5]

        #u_w = Wnb+wind.item(3)
        #v_w = Web+wind.item(4)
        #w_w = Wdb+wind.item(5) 

        u_r = u #- u_w
        v_r = v #- v_w
        w_r = w #- w_w

        self._Va = math.sqrt(u_r**2 + v_r**2 + w_r**2)
        
        # compute angle of attack      
        self._alpha = np.arctan2(w_r,u_r)        
        # compute sideslip angle
        self._beta = math.asin(v_r/self._Va)


    def _forces_moments(self, delta):
        """        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """  
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)

        # Gravity Force
        mass, g = MAV.mass, MAV.gravity
        Fgn = mass*g*2*(e1*e3-e2*e0)
        Fge = mass*g*2*(e2*e3+ e1*e0)
        Fgd = mass*g*(e3**2 + e0**2 - e1**2 -e2**2)
        Fg = [ Fgn, Fge, Fgd ]

        # Longitudinal Aerodynamics
        rho, Va, S, alpha, beta = MAV.rho, self._Va, MAV.S_wing, self._alpha, self._beta
        #-----Linear_model-----------
        # CL_alpha_ = MAV.C_L_0 + MAV.C_L_alpha*self._alpha
        # CD_alpha_ = MAV.C_D_0 + MAV.C_D_alpha*self._alpha
        #-----Non_Linear_model-------
        M, alpha_0 = MAV.M, MAV.alpha0
        blending_num_ = 1 + math.exp(-M*(alpha-alpha_0)) + math.exp(M*(alpha+alpha_0))
        blending_den_ = (1 + math.exp(-M*(alpha-alpha_0)))*(1 + math.exp(M*(alpha+alpha_0)))
        sigma_alpha_ = blending_num_/blending_den_
        CL_alpha_ = (1-sigma_alpha_)*(MAV.C_L_0 + MAV.C_L_alpha*self._alpha) + sigma_alpha_*(2*np.sign(alpha)*(np.sin(alpha))**2*np.cos(alpha))
        AR = MAV.b**2/S
        CL_alpha = np.pi*AR/(1 + np.sqrt(1 + (AR/2)**2))
        CD_alpha_ = MAV.C_D_p + (MAV.C_L_0 + CL_alpha*alpha)**2/(np.pi*MAV.e*AR)
        #-----------------------------
        c, s = np.cos, np.sin
        CX_alpha_ = -CD_alpha_*c(alpha) + CL_alpha_*s(alpha)
        CLq, CLde, CDq, CDde = MAV.C_L_q, MAV.C_L_delta_e, MAV.C_D_q, MAV.C_D_delta_e
        CXq_alpha_ = -CDq*c(alpha) + CLq*s(alpha)
        CXde_alpha_ = -CDde*c(alpha) + CLde*s(alpha)
        CZ_alpha_ = -CD_alpha_*s(alpha) - CL_alpha_*c(alpha)
        CZq_alpha_ = -CDq*s(alpha) - CLq*c(alpha)
        CZde_alpha_ = -CDde*s(alpha) - CLde*c(alpha)

        Cm0, Cma, Cmq, Cmde = MAV.C_m_0, MAV.C_m_alpha, MAV.C_m_q, MAV.C_m_delta_e
        cc, q, bb, p, r = MAV.c, self._state.item(11), MAV.b, self._state.item(10), self._state.item(12)
        d_e, d_t, d_a, d_r = delta[0], delta[1], delta[2], delta[3]

        Fx = 0.5*rho*(Va**2)*S*(CX_alpha_ + CXq_alpha_*cc*q/2/Va + CXde_alpha_)
        Fz = 0.5*rho*(Va**2)*S*(CZ_alpha_ + CZq_alpha_*cc*q/2/Va + CZde_alpha_*d_e)
        m = 0.5*rho*(Va**2)*S*cc*(Cm0 + Cma*alpha + Cmq*cc*q/Va/2 + Cmde*d_e)/MAV.Jy

        # Fx = 0.5*rho*(Va**2)*S*( (-CD_alpha_*c(alpha)+CL_alpha_*s(alpha)) + (-CDq*c(alpha)+CLq*s(alpha))*cc*q/2/Va + (-CDde*c(alpha)+CLde*s(alpha)*d_e) )
        # Fz = 0.5*rho*(Va**2)*S*( (-CD_alpha_*s(alpha)-CL_alpha_*c(alpha)) + (-CDq*s(alpha)-CLq*c(alpha))*cc*q/2/Va + (-CDde*s(alpha)-CLde*c(alpha)*d_e) )
        # m  = 0.5*rho*(Va**2)*S*cc*(Cm0 + Cma*alpha + Cmq*cc*q/Va/2 + Cmde*d_e)

        # Lateral Aerodynamics
        CY0, CYb, CYp, CYr, CYda, CYdr = MAV.C_Y_0, MAV.C_Y_beta, MAV.C_Y_p, MAV.C_Y_r, MAV.C_Y_delta_a, MAV.C_Y_delta_r
        Cl0, Clb, Clp, Clr, Clda, Cldr = MAV.C_ell_0, MAV.C_ell_beta, MAV.C_ell_p, MAV.C_ell_r, MAV.C_ell_delta_a, MAV.C_ell_delta_r
        Cn0, Cnb, Cnp, Cnr, Cnda, Cndr = MAV.C_n_0, MAV.C_n_beta, MAV.C_n_p, MAV.C_n_r, MAV.C_n_delta_a, MAV.C_n_delta_r
        Cp0, Cpb, Cpp, Cpr, Cpda, Cpdr = MAV.C_p_0, MAV.C_p_beta, MAV.C_p_p, MAV.C_p_r, MAV.C_p_delta_a, MAV.C_p_delta_r
        Cr0, Crb, Crp, Crr, Crda, Crdr = MAV.C_r_0, MAV.C_r_beta, MAV.C_r_p, MAV.C_r_r, MAV.C_r_delta_a, MAV.C_r_delta_r

        Fy = 0.5*rho*(Va**2)*S*(CY0 + CYb*beta + CYp*bb*p/Va /
                                2 + CYr*bb*r/Va/2 + CYda*d_a + CYdr*d_r)
        l = 0.5*rho*(Va**2)*S*bb*(Cp0 + Cpb*beta + Cpp*bb*p /
                                  Va/2 + Cpr*bb*r/Va/2 + Cpda*d_a + Cpdr*d_r)
        n = 0.5*rho*(Va**2)*S*bb*(Cr0 + Crb*beta + Crp*bb*p /
                                  Va/2 + Crr*bb*r/Va/2 + Crda*d_a + Crdr*d_r)

        # Fy = 0.5*rho*(Va**2)*S*( CY0 + CYb*beta + CYp*bb*p/Va/2 + CYr*bb*r/Va/2 + CYda*d_a +CYdr*d_r)
        # l  = 0.5*rho*(Va**2)*S*bb*(Cl0 + Clb*beta + Clp*bb*p/Va/2 + Clr*bb*r/Va/2 + Clda*d_a + Cldr*d_r)
        # n  = 0.5*rho*(Va**2)*S*bb*(Cn0 + Cnb*beta + Cnp*bb*p/Va/2 + Cnr*bb*r/Va/2 + Cnda*d_a + Cndr*d_r)



        D = MAV.D_prop
        C_Q2 = MAV.C_Q2  
        C_Q1 = MAV.C_Q1  
        C_Q0 = MAV.C_Q0  
        C_T2 = MAV.C_T2  
        C_T1 = MAV.C_T1  
        C_T0 = MAV.C_T0  
        K_V  = MAV.K_V                  # from datasheet RPM/V
        KQ   = MAV.KQ                   # KQ in N-m/A, V-s/rad
        R    = MAV.R_motor              # ohms
        V_in = MAV.V_max*d_t
        KQ_i0 = MAV.i0

        #Propeller Model in slides-----------

        a = rho*np.power(D,5)*C_Q0/(2*np.pi)**2
        b = rho*np.power(D,4)*C_Q1*self._Va/(2*np.pi) + KQ*K_V/R
        c = rho*np.power(D,3)*C_Q2*(self._Va)**2 - KQ*V_in/R + KQ_i0

        omega_op = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
        J_op = 2*np.pi*self._Va/omega_op/D 

        CT = C_T0 + C_T1*J_op + C_T2*J_op**2
        CQ = C_Q0 + C_Q1*J_op + C_Q2*J_op**2

        # n = omega_op/(2*np.pi)

        # Tp = rho*n**2*np.power(D,4)*CT
        # Qp = rho*n**2*np.power(D,5)*CQ

        Tp = rho*np.power(D,4)*omega_op**2*CT/(2*np.pi)**2
        Qp = rho*np.power(D,5)*omega_op**2*CQ/(2*np.pi)**2

        Ftx = Fx + Fgn + Tp
        Fty = Fy + Fge
        Ftz = Fz + Fgd
        lt = l + Qp

        self._forces[0] = Ftx
        self._forces[1] = Fty
        self._forces[2] = Ftz
        return np.array([[Ftx, Fty, Ftz, lt, m, n]]).T

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
        self.msg_true_state.Vg = math.sqrt((self._state[3][0])**2+(self._state[4][0])**2+(self._state[5][0])**2)

        e0 = self._state[6][0]
        e1 = self._state[7][0]
        e2 = self._state[8][0]
        e3 = self._state[9][0]
        u = self._state[3][0]
        v = self._state[4][0]
        w = self._state[5][0]

        pn_dot = (e1**2+e0**2-e2**2-e3**2)*u +2*(e1*e2-e3*e0)*v +2*(e1*e3+e2*e0)*w
        pe_dot = 2*(e1*e2 +e3*e0)*u + (e2**2+ e0**2 -e1**2-e3**2)*v +2*(e2*e3 - e1*e0)*w        
        pd_dot = 2*(e1*e3 -e2*e0)*u +2*(e2*e3 +e1*e0)*v +(e3**2+e0**2-e1**2-e2**2)*w

        # pn_dot, pe_dot, vg_k = Body2Inertia(phi, theta, psi, [self.state[3][0], self.state[4][0], self.state[5][0]] )

        self.msg_true_state.gamma = math.asin(-pe_dot/math.sqrt(pn_dot**2+pe_dot**2 + pd_dot**2))
        self.msg_true_state.chi = math.atan(pe_dot/pn_dot)
        self.msg_true_state.p = self._state.item(10)
        self.msg_true_state.q = self._state.item(11)
        self.msg_true_state.r = self._state.item(12)
        self.msg_true_state.wn = self._wind.item(0)
        self.msg_true_state.we = self._wind.item(1)
