import numpy as np
import math as m
pi = m.pi
s = m.sin
c = m.cos
root = m.sqrt

def Quaternion2Euler(e):
	e0 = e[0]
	e1 = e[1]
	e2 = e[2]
	e3 = e[3]
	phi = np.arctan2(2*(e0*e1+e2*e3), (e0*e0 + e3*e3 - e1*e1 - e2*e2))
	theta = np.arcsin(2*(e0*e2 -e1*e3))
	psi = np.arctan2(2*(e0*e3 + e1*e2), (e0*e0 + e1*e1 - e2*e2 - e3*e3))
	return [phi, theta, psi]

def Inertial2Body (phi, theta, psi, w_n, w_e, w_d):
	w_nb = c(theta)*c(psi)*w_n	+ c(theta)*s(psi)*w_e - s(theta)*w_d
	w_eb = (s(phi)*s(theta)*c(psi) - c(phi)*s(psi))*w_n + (s(phi)*s(theta)*s(psi)+c(phi)*c(psi))*w_e + s(phi)*c(theta)*w_d
	w_db = (c(phi)*s(theta)*c(psi)+s(phi)*s(psi))*w_n + (c(phi)*s(theta)*s(psi)-s(phi)*c(psi))*w_e + c(theta)*c(phi)*w_d

	wind_b = [w_nb, w_eb, w_db]
	return wind_b

def Quaternion2Rotation(quat):
    e0 = quat[0]
    e1 = quat[1]
    e2 = quat[2]
    e3 = quat[3]

    R = np.array([[e1**2 + e0**2 - e2**2 -e3**2, 2*(e1*e2 - e3*e0), 2*(e1*e3 + e2*e0)],
                  [2*(e1*e2 + e3*e0), e2**2 + e0**2 - e1**2 - e3**2, 2*(e2*e3 - e1*e0)],
                  [2*(e1*e3 - e2*e0), 2*(e2*e3 + e1*e0), e3**2 + e0**2 - e1**2 - e2**2]])

    return R

def Euler2Rotation(phi, theta, psi):
    """
    Converts euler angles to rotation matrix (R_b^i, i.e., body to inertial)
    """
    # only call sin and cos once for each angle to speed up rendering
    c_phi = np.cos(phi)
    s_phi = np.sin(phi)
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)

    R_roll = np.array([[1, 0, 0],
                       [0, c_phi, s_phi],
                       [0, -s_phi, c_phi]])
    R_pitch = np.array([[c_theta, 0, -s_theta],
                        [0, 1, 0],
                        [s_theta, 0, c_theta]])
    R_yaw = np.array([[c_psi, s_psi, 0],
                      [-s_psi, c_psi, 0],
                      [0, 0, 1]])

    R = R_roll @ R_pitch @ R_yaw  # inertial to body (Equation 2.4 in book)
    return R.T  # transpose to return body to inertial