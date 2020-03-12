"""
observer
    - Beard & McLain, PUP, 2012
    - Last Update:
        3/2/2019 - RWB
"""
import sys
import numpy as np
sys.path.append('..')
import parameters.control_parameters as CTRL
import parameters.simulation_parameters as SIM
import parameters.sensor_parameters as SENSOR
import parameters.aerosonde_parameters as MAV
from tools.rotations import Euler2Rotation
from tools.wrap import wrap

from message_types.msg_state import msg_state

class observer:
    def __init__(self, ts_control):
        # initialized estimated state message
        self.estimated_state = msg_state()
        # use alpha filters to low pass filter gyros and accels
        self.lpf_gyro_x = alpha_filter(alpha=0.5)
        self.lpf_gyro_y = alpha_filter(alpha=0.5)
        self.lpf_gyro_z = alpha_filter(alpha=0.5)
        self.lpf_accel_x = alpha_filter(alpha=0.5)
        self.lpf_accel_y = alpha_filter(alpha=0.5)
        self.lpf_accel_z = alpha_filter(alpha=0.5)
        # use alpha filters to low pass filter static and differential pressure
        self.lpf_static = alpha_filter(alpha=0.9)
        self.lpf_diff = alpha_filter(alpha=0.5)
        # ekf for phi and theta
        self.attitude_ekf = ekf_attitude()
        # ekf for pn, pe, Vg, chi, wn, we, psi
        self.position_ekf = ekf_position()

    def update(self, measurements):

        # estimates for p, q, r are low pass filter of gyro minus bias estimate
        self.estimated_state.p = self.lpf_gyro_x.update(measurements.gyro_x)
        self.estimated_state.q = self.lpf_gyro_y.update(measurements.gyro_y)
        self.estimated_state.r = self.lpf_gyro_z.update(measurements.gyro_z)

        # invert sensor model to get altitude and airspeed
        self.estimated_state.h = self.lpf_static.update(measurements.static_pressure)/MAV.rho/MAV.gravity
        self.estimated_state.Va = np.sqrt(2*self.lpf_diff.update(measurements.diff_pressure)/MAV.rho)

        # estimate phi and theta with simple ekf
        self.attitude_ekf.update(self.estimated_state, measurements)

        # estimate pn, pe, Vg, chi, wn, we, psi
        self.position_ekf.update(self.estimated_state, measurements)

        # not estimating these
        self.estimated_state.alpha = self.estimated_state.theta
        self.estimated_state.beta = 0.0
        self.estimated_state.bx = 0.0
        self.estimated_state.by = 0.0
        self.estimated_state.bz = 0.0
        return self.estimated_state

class alpha_filter:
    # alpha filter implements a simple low pass filter
    # y[k] = alpha * y[k-1] + (1-alpha) * u[k]
    def __init__(self, alpha=0.5, y0=0.0):
        self.alpha = alpha  # filter parameter
        self.y = y0  # initial condition

    def update(self, u):
        self.y = self.alpha*self.y + (1-self.alpha)*u
        return self.y

class ekf_attitude:
    # implement continous-discrete EKF to estimate roll and pitch angles
    def __init__(self):
        self.Q = 1e-9*np.diag(1,1)
        self.Q_gyro = SENSOR.gyro_sigma**2*np.diag(1,1,1)
        self.R_accel = SENSOR.accel_sigma**2*np.diag(1,1,1)
        self.N = 2  # number of prediction step per sample
        self.xhat =  np.array([0,0]).T # initial state: phi, theta
        self.P = np.diag(10,10)
        self.Ts = SIM.ts_control/self.N

    def update(self, state, measurement):
        self.propagate_model(state)
        self.measurement_update(state, measurement)
        state.phi = self.xhat.item(0)
        state.theta = self.xhat.item(1)

    def f(self, x, state):  
        # system dynamics for propagation model: xdot = f(x, u)
        p,q,r,phi,theta = state.p,state.q, state.r, x.phi, x.theta

        _f = np.array([ p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta),
                        q*np.cos(phi)-r*np.sin(phi) ]).T
        return _f

    def h(self, x, state):
        # measurement model y
        g = MAV.gravity
        p, q, r, phi, theta, Va= state.p, state.q, state.r, x.phi, x.theta, state.Va
        _h = np.array([
                    q*Va*np.sin(theta) + MAV.gravity*np.sin(theta),
                    r*Va*np.cos(theta) - p*Va*np.sin(theta) - g*np.cos(theta)*np.sin(phi),
                    -q*Va*np.cos(theta) - g*np.cos(theta)*np.cos(phi)
        ]).T 
        return _h

    def propagate_model(self, state):
        # model propagation
        for i in range(0, self.N):
             # propagate model
            self.xhat = self.xhat + self.Ts*self.f(self.xhat, state)
            # compute Jacobian
            A = jacobian(self.f, self.xhat, state)
            # compute G matrix for gyro noise
            theta,phi = state.theta,state.phi
            G = np.array([ [1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta), 0],
                            [0, np.cos(phi), -np.sin(phi), 0]  ])
            # update P with continuous time model
            # self.P = self.P + self.Ts * (A @ self.P + self.P @ A.T + self.Q + G @ self.Q_gyro @ G.T)
            # convert to discrete time models
            A_d = np.eye(2) + A*self.Ts + (self.Ts**2)/2*(A@A)
            # G_d = 
            # update P with discrete time model
            self.P = A_d @ self.P @ A_d.T + (self.Q + G @ self.Q_gyro @ G.T)*self.Ts**2 

    def measurement_update(self, state, measurement):
        # measurement updates
        threshold = 2.0
        h = self.h(self.xhat, state)
        C = jacobian(self.h, self.xhat, state)
        y = np.array([measurement.accel_x, measurement.accel_y, measurement.accel_z]).T
        S_inv = np.linalg.inv(self.R_accel+ C @ self.P @ C.T)
        # for i in range(0, 3):
        if np.abs(y[0]-h[0, 0]) or np.abs(y[1]-h[1, 0]) or np.abs(y[2]-h[2, 0]) < threshold:
            # Ci = 
            L = self.P @ C.T @ S_inv
            I_LC = np.eye(2)-L @ C
            self.P = I_LC @ self.P @ I_LC.T + L @ self.R_accel @ L.T
            self.xhat = self.xhat + L @ (y-h)

class ekf_position:
    # implement continous-discrete EKF to estimate pn, pe, chi, Vg
    def __init__(self):
        self.Q = 1e-9*np.diag(1,1,1,1,1,1,1)
        self.R_gps = np.diag(SENSOR.gps_n_sigma**2, SENSOR.gps_e_sigma**2, SENSOR.gps_Vg_sigma**2, SENSOR.gps_course_sigma**2)
        self.R_psudo = np.diag(0.01**2, 0.01**2)
        self.N = 2   # number of prediction step per sample
        self.Ts = (SIM.ts_control / self.N)
        self.xhat = np.array([0,0,0,0,0,0,0]).T
        self.P = np.diag(10,10,10,10,10,10,10)
        self.gps_n_old = 9999
        self.gps_e_old = 9999
        self.gps_Vg_old = 9999
        self.gps_course_old = 9999


    def update(self, state, measurement):
        self.propagate_model(state)
        self.measurement_update(state, measurement)
        state.pn = self.xhat.item(0)
        state.pe = self.xhat.item(1)
        state.Vg = self.xhat.item(2)
        state.chi = self.xhat.item(3)
        state.wn = self.xhat.item(4)
        state.we = self.xhat.item(5)
        state.psi = self.xhat.item(6)

    def f(self, x, state):
        # system dynamics for propagation model: xdot = f(x, u)
        Vg, chi, wn, we, psi = x.item(2), x.item(3), x.item(4), x.item(5), x.item(6)
        g = MAV.gravity
        phi, q, r, theta, Va = state.phi, state.q, state.r, state.theta, state.Va
        pn_dot = Vg*np.cos(chi)
        pe_dot = Vg*np.sin(chi)
        chi_dot = (g/Vg)*np.tan(phi)*np.cos(chi - psi)
        wn_dot = 0
        we_dot = 0
        psi_dot = q*(np.sin(phi)/np.cos(theta)) + r*(np.cos(phi)/np.cos(theta))
        Vg_dot = ((Va*np.cos(psi)+wn)*(-Va*psi_dot*np.sin(psi)) + (Va*np.sin(psi)+we)*Va*psi_dot*np.cos(psi))/Vg
        _f = np.array([
            pn_dot, pe_dot, Vg_dot, chi_dot, wn_dot, we_dot, psi_dot  
        ]).T 

        return _f

    def h_gps(self, x, state):
        # measurement model for gps measurements
        pn, pe, Vg, chi = x.item(0), x.item(1), x.item(2), x.item(3)    
        _h =np.array([
            pn,
            pe,
            Vg,
            chi,
        ]).T 
        return _h


    def h_pseudo(self, x, state):
        # measurement model for wind triangale pseudo measurement
        Vg, chi, wn, we, psi = x.item(2), x.item(3), x.item(4), x.item(5), x.item(6)        
        Va = state.Va
        _h = np.array([            
            Va*np.cos(psi)+wn-Vg*np.cos(chi),
            Va*np.cos(psi) + we - Vg*np.sin(chi)
        ]).T
        return _h

    def propagate_model(self, state):
        # model propagation
        for i in range(0, self.N):
            # propagate model
            self.xhat = self.xhat + self.Ts*self.f(self.xhat,state)
            # compute Jacobian
            A = jacobian(self.f, self.xhat, state)
            # update P with continuous time model
            # self.P = self.P + self.Ts * (A @ self.P + self.P @ A.T + self.Q + G @ self.Q_gyro @ G.T)
            # convert to discrete time models
            A_d = np.eye(7) + A @ self.Ts + A @ A * self.Ts**2
            # update P with discrete time model
            self.P = A_d @ P @ A_d.T + self.Ts**2*self.Q 

    def measurement_update(self, state, measurement):
        # always update based on wind triangle pseudu measurement
        h = self.h_pseudo(self.xhat, state)
        C = jacobian(self.h_pseudo, self.xhat, state)
        y = np.array([0, 0]).T 
        # for i in range(0, 2):
            # Ci = 
        S_inv = np.linalg.inv(self.R_psudo + C @ self.P @ C.T)
        L = self.P @ C.T @ S_inv
        I_LC = np.eye(7) - L @ C
        self.P = I_LC @ self.P @ I_LC.T + L @ self.R @ L.T
        self.xhat = self.xhat + L @ (y - h)

        # only update GPS when one of the signals changes
        if (measurement.gps_n != self.gps_n_old) \
            or (measurement.gps_e != self.gps_e_old) \
            or (measurement.gps_Vg != self.gps_Vg_old) \
            or (measurement.gps_course != self.gps_course_old):

            h = self.h_gps(self.xhat, state)
            C = jacobian(self.h_gps, self.xhat, state)
            y = np.array([measurement.gps_n, measurement.gps_e, measurement.gps_Vg, measurement.gps_course]).T
            # for i in range(0, 4):
                # Ci = 
            L = self.P @ C.T @ S_inv
            I_LC = np.eye(7) - L @ C
            self.P = self.P = I_LC @ self.P @ I_LC.T + L @ self.R_gps @ L.T
            self.xhat = self.xhat = self.xhat + L @ (y - h)
            # update stored GPS signals
            self.gps_n_old = measurement.gps_n
            self.gps_e_old = measurement.gps_e
            self.gps_Vg_old = measurement.gps_Vg
            self.gps_course_old = measurement.gps_course

def jacobian(fun, x, state):
    # compute jacobian of fun with respect to x
    f = fun(x, state)
    m = f.shape[0]
    n = x.shape[0]
    eps = 0.01  # deviation
    J = np.zeros((m, n))
    for i in range(0, n):
        x_eps = np.copy(x)
        x_eps[i][0] += eps
        f_eps = fun(x_eps, state)
        df = (f_eps - f) / eps
        J[:, i] = df[:, 0]
    return J
