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


def Euler2Quaternion(phi, theta, psi):		
	e0 = c(psi/2)*c(theta/2)*c(phi/2) + s(psi/2)*s(theta/2)*s(phi/2)
	e1 = c(psi/2)*c(theta/2)*s(phi/2) - s(psi/2)*s(theta/2)*c(phi/2)
	e2 = c(psi/2)*s(theta/2)*c(phi/2) + s(psi/2)*c(theta/2)*s(phi/2)	
	e3 = s(psi/2)*c(theta/2)*c(phi/2) - c(psi/2)*s(theta/2)*s(psi/2)
	return [e0, e1, e2, e3]

#Dryden Gust Model
def parameter_calc():
	Lu, Lv, Lw =200, 200, 50
	sigma_u, sigma_v  = 1.06, 1.06
	sigma_w = .7
	V_a = 20

	c_u = sigma_u*root(2*V_a/Lu)
	c_v = sigma_v*root(3*V_a/Lv)
	c_w = sigma_w*root(3*V_a/Lw)

	a_1 = c_u
	b_1 = V_a/Lu
	a_2 = c_v
	a_3 = c_v*V_a/Lv/root(3)
	b_2 = V_a/Lv
	a_4 = c_w
	a_5 = c_w*V_a/Lw/root(3)
	b_3 = V_a/Lw

	params = [a_1, b_1, a_2, a_3, b_2, a_4, a_5, b_3 ]
	return params


def Inertial2Body (phi, theta, psi, w_n, w_e, w_d):
	w_nb = c(theta)*c(psi)*w_n	+ c(theta)*s(psi)*w_e - s(theta)*w_d
	w_eb = (s(phi)*s(theta)*c(psi) - c(phi)*s(psi))*w_n + (s(phi)*s(theta)*s(psi)+c(phi)*c(psi))*w_e + s(phi)*c(theta)*w_d
	w_db = (c(phi)*s(theta)*c(psi)+s(phi)*s(psi))*w_n + (c(phi)*s(theta)*s(psi)-s(phi)*c(psi))*w_e + c(theta)*c(phi)*w_d

	wind_b = [w_nb, w_eb, w_db]
	return wind_b
		