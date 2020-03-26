import numpy as np
from math import sin, cos, atan, atan2
import sys

sys.path.append('..')
from message_types.msg_autopilot import msg_autopilot

class path_follower:
    def __init__(self):
        self.chi_inf = np.radians(80)  # approach angle for large distance from straight-line path
        self.k_path = 0.02  # proportional gain for straight-line path following
        self.k_orbit = 2.5  # proportional gain for orbit following
        self.gravity = 9.8
        self.autopilot_commands = msg_autopilot()  # message sent to autopilot

    def update(self, path, state):
        if path.flag=='line':
            self._follow_straight_line(path, state)
        elif path.flag=='orbit':
            self._follow_orbit(path, state)
        return self.autopilot_commands

    def _follow_straight_line(self, path, state):
        #Course command
        p = np.array([[state.pn], [state.pe], [-state.h]])
        r = path.line_origin
        q = path.line_direction
        Xq = np.arctan2(q[1][0],q[0][0])
        cXq, sXq = np.cos(Xq), np.sin(Xq)
        R = np.array([
            [ cXq, sXq, 0],
            [-sXq, cXq, 0],
            [0, 0, 0]
        ])
        e_p = R @ (p - r)
        #Altitude command        
        ep_i = np.array([state.pn, state.pe, -state.h]) - np.array([r[0][1],r[1][0], r[2][0]])
        cross = np.cross([q[0][0],q[1][0],q[2][0]], [0,0,1])
        n = cross/np.linalg.norm(cross)
        s = ep_i - np.dot(ep_i,n)*n


        self.autopilot_commands.airspeed_command = path.airspeed
        self.autopilot_commands.course_command = Xq-self.chi_inf*(2/np.pi)*np.arctan(self.k_path * e_p[1][0])   
        self.autopilot_commands.altitude_command = -r[2][0] - np.sqrt(s[0]**2+s[1]**2)*q[2][0]/np.sqrt(q[0][0]**2+q[1][0]**2)
        self.autopilot_commands.phi_feedforward = 0

    def _follow_orbit(self, path, state):
        c = path.orbit_center
        d = np.sqrt((state.pn-c[0][0])**2 + (state.pe-c[1][0])**2)
        rho = path.orbit_radius
        if path.orbit_direction =='CW':
            l= 1
        else:
            l=-1
        psi = np.arctan2(state.pe-c[1][0], state.pn-c[0][0])
        psi = self._wrap(psi, state.chi)

        Vg = state.Vg
        chi = state.chi
        psi = state.psi
        phi_feedfw = atan(Vg**2 / ( self.gravity * rho * cos( chi - psi ) ))

        chi_c = psi + l*(np.pi/2  +  np.arctan(self.k_orbit)*(d-rho)/rho)
        self.autopilot_commands.airspeed_command = path.airspeed
        self.autopilot_commands.course_command = chi_c
        self.autopilot_commands.altitude_command = -path.orbit_center[2][0]
        self.autopilot_commands.phi_feedforward = phi_feedfw

    def _wrap(self, chi_c, chi):
        while chi_c-chi > np.pi:
            chi_c = chi_c - 2.0 * np.pi
        while chi_c-chi < -np.pi:
            chi_c = chi_c + 2.0 * np.pi
        return chi_c

