import math as m
pi = m.pi
sin = m.sin
cos = m.cos

class equations:
	def __init__(self, psi, theta, phi, u, v, w, jx, jy, jz, jxz, p, q, r, l, m, n):
		self.psi = psi
		self.theta = theta
		self. phi = phi
		self.u = u
		self.v = v
		self.w = w
		self.jx = jx
		self.jy = jy
		self.jz = jz
		self.jxz = jxz
		self.p = p
		self. q = q
		self.r = r
		self.l = l
		self.m = m
		self.n = n


	def  pn_dot(self):
		t1 = cos(self.theta)*cos(self.psi)*self.u
		t2 = (sin(self.phi)*sin(self.theta)*cos(self.psi)-cos(self.phi)*sin(self.psi))*self.v
		t3 = (cos(self.psi)*sin(self.theta)*cos(self.psi)+sin(self.phi)*sin(self.psi))*self.w
		return t1 + t2 + t3

	def pe_dot(self):
		t1 = cos(self.theta)*sin(self.psi)*self.u
		t2 = (sin(self.phi)*sin(self.theta)*sin(self.psi)+cos(self.phi)*cos(self.psi))*self.v
		t3 = (cos(self.phi)*sin(self.theta)*sin(self.psi)-sin(self.phi)*cos(self.psi))*self.w
		return t1 + t2 + t3

	def pd_dot(self):
		t1 = -sin(self.theta)*self.u
		t2 = sin(self.phi)*cos(self.theta)*self.v
		t3 = cos(self.phi)*cos(self.theta)*self.w
		return t1 + t2 + t3

	def rot_dyn(self):
		T = self.jx*self.jz-self.jxz*self.jxz
		T1 = (self.jx-self.jy+self.jz)*self.jxz/T
		T2 = ((self.jz-self.jy)*self.jz+self.jxz*self.jxz)/T
		T3 = self.jz/T
		T4 = self.jxz/T
		T5 = (self.jz-self.jx)/self.jy
		T6 = self.jxz/self.jy
		T7 = ((self.jx - self.jy)*self.jx +self.jxz*self.jxz)/T
		T8 = self.jx/T

		p_dot = T1*self.p*self.q - T2*self.q*self.r +T3*self.l +T4*self.n
		q_dot = T5*self.p*self.r -T6*(self.p*self.p - self.r*self.r) + self.m/self.jy
		r_dot = T7*self.p*self.q -T1*self.q*self.r + T4*self.l + T8*self.n

		dot_values = [p_dot, q_dot, r_dot ]
		return dot_values

