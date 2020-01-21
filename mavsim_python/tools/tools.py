import numpy as np


class Quaternion2Euler:	
	def __init__(self, e0, e1, e2, e3 ):
		self.e0 = e0
		self.e1 = e1
		self.e2 = e2
		self.e3 = e3


	def calc_phi(self):
		phi = np.arctan2(
				2*(self.e0*self.e1+self.e2*self.e3), (self.e0*self.e0 + self.e3*self.e3 - self.e1*self.e1 - self.e2*self.e2)
			)
		return phi

	def calc_theta(self):
		theta = np.arcsin(
					2*(self.e0*self.e2 -self.e1*self.e3)
			)
		return theta

	def  calc_psi(self):
		psi = np.arctan2(
				2*(self.e0*self.e3 + self.e1*self.e2), (self.e0*self.eo + self.e1*self.e1 - self.e2*self.e2 - self.e3*self.e3)
			)
		return psi

	def calc_angles():
		angles = [calc_phi(), calc_theta(), calc_psi()]
		return angles
		