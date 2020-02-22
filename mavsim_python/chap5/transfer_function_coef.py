import numpy as np
import parameters.aerosonde_parameters as MAV
import sys
sys.path.append('..')

CY0, CYb, CYp, CYr, CYda, CYdr = MAV.C_Y_0, MAV.C_Y_beta, MAV.C_Y_p, MAV.C_Y_r, MAV.C_Y_delta_a, MAV.C_Y_delta_r
Cp0, Cpb, Cpp, Cpr, Cpda, Cpdr = MAV.C_p_0, MAV.C_p_beta, MAV.C_p_p, MAV.C_p_r, MAV.C_p_delta_a, MAV.C_p_delta_r
rho, S, mass, bb = MAV.rho, MAV.S_wing, MAV.mass, MAV.b
Va = 17

#Roll_Autopilot
a_phi_1 = - (0.25*rho*(Va)*S*(bb**2)*Cpp)
a_phi_2 = 0.5*rho*(Va**2)*S*bb*Cpda

#Slideslip_Autopilot
a_beta_1 = rho*Va*S*CYb/2/mass
a_beta_2 = rho*Va*S*CYdr/2/mass

#Pitch_Autopilot
Cm0, Cma, Cmq, Cmde = MAV.C_m_0, MAV.C_m_alpha, MAV.C_m_q, MAV.C_m_delta_e
c = MAV.c
a_theta_1 = -rho*Va**2*c*S*Cmq*c/2/MAV.Jy/2/Va
a_theta_2 = -rho*Va**2*c*S*Cma/2/MAV.Jy
a_theta_3 = rho*Va**2*c*S*Cmde/2/MAV.Jy

#Airspeed_Coefficients
a_v1 = rho*Va*S*(MAV.C_D_0 + MAV.C_D_alpha*4.4*np.pi/180 +
                 MAV.C_D_delta_e*.95) + rho*S*MAV.C_prop*Va/mass
a_v2 = rho*S*MAV.C_prop*((MAV.k_motor/100)**2)*0.95/mass

print(a_v1)
print(a_v2)
