import numpy as np
import math as m
pi = m.pi
sin = m.sin
cos = m.cos


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
	e0 = cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2)
	e1 = cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2)
	e2 = cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2)	
	e3 = sin(psi/2)*cos(theta/2)*cos(phi/2) - cos(psi/2)*sin(theta/2)*sin(psi/2)
	return [e0, e1, e2, e3]








			
		