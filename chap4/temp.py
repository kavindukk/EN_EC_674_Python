phi, theta, psi = Quaternion2Euler(
    [self._state[6], self._state[7], self._state[8], self._state[9]])

# C_X_alpha = -MAV.C_D_alpha*m.cos(self._alpha) + MAV.C_L_alpha*m.sin(self._alpha)
# C_X_q_alpha = -MAV.C_D_q*m.cos(self._alpha) + MAV.C_L_q*m.sin(self._alpha)
# C_X_delta_e_alpha = -MAV.C_D_delta_e*m.cos(self._alpha) + MAV.C_L_delta_e*m.sin(self._alpha)
# C_Z_alpha = -MAV.C_D_alpha*m.sin(self._alpha) - MAV.C_L_alpha*m.cos(self._alpha)
# C_Z_q_alpha = -MAV.C_D_q*m.sin(self._alpha) - MAV.C_L_q*m.cos(self._alpha)
# C_Z_delta_e_alpha = -MAV.C_D_delta_e*m.sin(self._alpha) - MAV.C_L_delta_e*m.cos(self._alpha)

# fx = -mass*g*m.sin(theta) + .5*rho*Va**2*S*(C_X_alpha + C_X_q_alpha*c*q/2/Va  + C_X_delta_e_alpha*delta[0] )
# fy = mass*g*m.cos(theta)*m.sin(phi) + .5*rho*Va**2*S*(MAV.C_Y_0 + MAV.C_Y_beta*self._beta + MAV.C_Y_p*MAV.b*self._state[12]/2/Va + MAV.C_Y_delta_a*delta[2] + MAV.C_Y_delta_r*delta[3] )
# fz = mass*g*m.cos(theta)*m.cos(phi) + .5*rho*Va**2*S*(C_Z_alpha + C_Z_q_alpha*c*self._state[11]/2/Va + C_Z_delta_e_alpha*delta[0])

# fx = 2*mass*g*(e1*e3-e2*e0) + .5*rho*Va**2*S*(C_X_alpha + C_X_q_alpha*c*q/2/Va  + C_X_delta_e_alpha*delta[0] )
# fy = 2*mass*g*(e2*e3+e1*e0) + .5*rho*Va**2*S*(MAV.C_Y_0 + MAV.C_Y_beta*self._beta + MAV.C_Y_p*MAV.b*self._state[12]/2/Va + MAV.C_Y_delta_a*delta[2] + MAV.C_Y_delta_r*delta[3] )
# fz = mass*g*(e3**2+e0**2-e1**2-e2**2) + .5*rho*Va**2*S*(C_Z_alpha + C_Z_q_alpha*c*self._state[11]/2/Va + C_Z_delta_e_alpha*delta[0])

# C_L_0 = 0.28
# C_L_alpha = 3.45
# C_L_q = 0.0
# C_L_delta_e = -0.36
# C_D_0 = 0.03
# C_D_alpha = 0.3
# C_D_p = 0.0437
# C_D_q = 0.0
# C_D_delta_e = 0.0
