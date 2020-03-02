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