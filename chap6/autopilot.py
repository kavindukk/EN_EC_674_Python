"""
autopilot block for mavsim_python
    - Beard & McLain, PUP, 2012
    - Last Update:
        2/6/2019 - RWB
"""
import sys
import numpy as np
sys.path.append('..')
import parameters.control_parameters as AP
from tools.transfer_function import transfer_function
from tools.wrap import wrap
from chap6.pid_control import pid_control#, pi_control, pd_control_with_rate
from message_types.msg_state import msg_state


class autopilot:
    def __init__(self, ts_control):
        # instantiate lateral controllers
        self.roll_from_aileron = pid_control(
                        kp=AP.roll_kp,
                        kd=AP.roll_kd,
                        limit=np.radians(45))
        self.course_from_roll = pid_control(
                        kp=AP.course_kp,
                        ki=AP.course_ki,
                        Ts=ts_control,
                        limit=np.radians(30))
        # self.sideslip_from_rudder = pid_control(
        #                 kp=AP.sideslip_kp,
        #                 ki=AP.sideslip_ki,
        #                 Ts=ts_control,
        #                 limit=np.radians(45))
        self.yaw_damper = transfer_function(
                        num=np.array([[AP.yaw_damper_kp, 0]]),
                        den=np.array([[1, 1/AP.yaw_damper_tau_r]]),
                        Ts=ts_control)

        # instantiate lateral controllers
        self.pitch_from_elevator = pid_control(
                        kp=AP.pitch_kp,
                        kd=AP.pitch_kd,
                        limit=np.radians(45))
        self.altitude_from_pitch = pid_control(
                        kp=AP.altitude_kp,
                        ki=AP.altitude_ki,
                        Ts=ts_control,
                        limit=np.radians(30))
        self.airspeed_from_throttle = pid_control(
                        kp=AP.airspeed_throttle_kp,
                        ki=AP.airspeed_throttle_ki,
                        Ts=ts_control,
                        limit=1.5,
                        throttle_flag=True)
        self.commanded_state = msg_state()

    def update(self, cmd, state):
        # lateral autopilot
        # cmd.course_command= wrap(cmd.course_command, state.chi)        
        # phi_c = self.course_from_roll.update(cmd.course_command, state.chi, reset_flag=True)
        phi_c = self.course_from_roll.update(cmd.course_command, state.chi)
        delta_a = self.roll_from_aileron.update_with_rate(phi_c, state.phi, state.p)
        delta_a = np.asscalar(delta_a)
        # delta_r = self.sideslip_from_rudder.update(0, state.beta)
        delta_r = self.yaw_damper.update(state.r)

        # longitudinal autopilot
        h_c = self.saturate(cmd.altitude_command, state.h-AP.altitude_zone, state.h + AP.altitude_zone)
        theta_c = self.altitude_from_pitch.update(h_c, state.h)
        delta_e = self.pitch_from_elevator.update_with_rate(theta_c, state.theta, state.q)
        # delta_e = np.asscalar(delta_e)
        # delta_e = delta_e.item(0)
        delta_t = self.airspeed_from_throttle.update(cmd.airspeed_command, state.Va)
        delta_t = self.saturate(delta_t, 0.0, 1)

        # construct output and commanded states
        # delta = np.array([[delta_e], [delta_a], [delta_r], [delta_t]]) 
        delta = np.array([[delta_e], [delta_t], [delta_a], [delta_r]])
        self.commanded_state.h = cmd.altitude_command
        self.commanded_state.Va = cmd.airspeed_command
        self.commanded_state.phi = phi_c
        self.commanded_state.theta = theta_c
        self.commanded_state.chi = cmd.course_command
        
        return delta, self.commanded_state

    def saturate(self, input, low_limit, up_limit):
        if input <= low_limit:
            output = low_limit
        elif input >= up_limit:
            output = up_limit
        else:
            output = input
        return output
