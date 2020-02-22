import parameters.aerosonde_parameters as MAV
import sys
sys.path.append('..')
import numpy as np
import chap5.transfer_function_coef as TF

gravity = MAV.gravity
sigma = 0.05
Va0 = 15

#----------roll loop-------------
roll_kp = 45/15
roll_kd = .4

#----------course loop-------------
course_kp = 2.8
course_ki = -0.02029

#----------sideslip loop-------------
sideslip_ki = 1
sideslip_kp = 1

#----------yaw damper-------------
yaw_damper_tau_r = 1
yaw_damper_kp = 1

#----------pitch loop-------------
pitch_kp = -45/10
pitch_kd = -1.16364
K_theta_DC = pitch_kp*TF.a_theta_3/(TF.a_theta_2 + pitch_kp*TF.a_theta_3)

#----------altitude loop-------------
altitude_kp = 0.9723
altitude_ki = 0.06876
altitude_zone = 10

#---------airspeed hold using throttle---------------
airspeed_throttle_kp = (TF.a_v1 - 2*.707*1)/K_theta_DC/MAV.gravity
airspeed_throttle_ki = 1/K_theta_DC/MAV.gravity
