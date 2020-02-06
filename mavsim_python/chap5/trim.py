"""
compute_trim 
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        2/5/2019 - RWB
"""
from tools.tools import Euler2Quaternion, Quaternion2Euler
import sys
sys.path.append('..')
import numpy as np
from scipy.optimize import minimize
from parameters import aerosonde_parameters as a


def compute_trim(mav, Va, gamma):
    # define initial state and input

    #   Initial conditions for MAV
    pn0 = 0.  # initial north position
    pe0 = 0.  # initial east position
    pd0 = -20.0  # initial down position
    u0 = 25.  # initial velocity along body x-axis
    v0 = 0.  # initial velocity along body y-axis
    w0 = 0.  # initial velocity along body z-axis
    phi0 = 0.  # initial roll angle
    theta0 = 0.  # initial pitch angle
    psi0 = 0.0  # initial yaw angle
    p0 = 0  # initial roll rate
    q0 = 0  # initial pitch rate
    r0 = 0  # initial yaw rate
    Va0 = np.sqrt(u0 ** 2 + v0 ** 2 + w0 ** 2)
    #   Quaternion State
    e = Euler2Quaternion(phi0, theta0, psi0)
    e0 = e[0]
    e1 = e[1]
    e2 = e[2]
    e3 = e[3]
    state0 = np.array([a.pn0, a.pe0, a.pd0, a.u0, a.v0, a.w0, a.e0, a.e1, a.e2, a.e3, a.p0, a.q0, a.r0]).T
    delta0 = np.array([0., 0., 0., 0.]).T
    x0 = np.concatenate((state0, delta0), axis=0)
    # define equality constraints
    cons = ({'type': 'eq',
             'fun': lambda x: np.array([
                                x[3]**2 + x[4]**2 + x[5]**2 - Va**2,  # magnitude of velocity vector is Va
                                x[4],  # v=0, force side velocity to be zero
                                x[6]**2 + x[7]**2 + x[8]**2 + x[9]**2 - 1.,  # force quaternion to be unit length
                                x[7], # e1=0  - forcing e1=e3=0 ensures zero roll and zero yaw in trim
                                x[9], # e3=0
                                x[10], # p=0  - angular rates should all be zero
                                x[11], # q=0
                                x[12], # r=0
                                ]),
             'jac': lambda x: np.array([
                                [0., 0., 0., 2*x[3], 2*x[4], 2*x[5], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 2*x[6], 2*x[7], 2*x[8], 2*x[9], 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
                                ])
             })
    # solve the minimization problem to find the trim states and inputs
    res = minimize(trim_objective, x0, method='SLSQP', args = (mav, Va, gamma),
                   constraints=cons, options={'ftol': 1e-10, 'disp': True})
    # extract trim state and input and return
    trim_state = np.array([res.x[0:13]]).T
    trim_input = np.array([res.x[13:17]]).T
    return trim_state, trim_input

# objective function to be minimized
def trim_objective(x, mav, Va, gamma):
    R = 150;

    state = x[0:13]
    delta = x[13:17]
    x_dot = np.array([0., 0., Va*np.sin(gamma), 0., 0., 0., 0., 0., Va*np.sin(gamma)/R, 0., 0., 0.])
    mav._state = state
    mav._update_velocity_data()
    forces_moments = mav._forces_moments(delta)
    f = mav._derivatives(state, forces_moments)
    f_changed = [f[0], f[1], -f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12]]
    # tef = x_dot - f
    tef = x_dot - f_changed
    J = np.linalg.norm(tef)

    return J

