import numpy as np
import sys
sys.path.append('..')
from chap11.dubins_parameters import dubins_parameters
from message_types.msg_path import msg_path
import tools.wrap

class path_manager:
    def __init__(self):
        # message sent to path follower
        self.path = msg_path()
        # pointers to previous, current, and next waypoints
        self.ptr_previous = 0
        self.ptr_current = 1
        self.ptr_next = 2        
        # flag that request new waypoints from path planner
        self.flag_need_new_waypoints = True
        self.halfspace_n = np.inf * np.ones((3,1))
        self.halfspace_r = np.inf * np.ones((3,1))
        # state of the manager state machine
        self.manager_state = 1
        # dubins path parameters
        self.dubins_path = dubins_parameters()
        self.new_waypoint_path = True 
        self.flag_path_changed = False 
        self.state_changed = True

    def update(self, waypoints, radius, state):
        # this flag is set for one time step to signal a redraw in the viewer
        if self.path.flag_path_changed == True:
            self.path.flag_path_changed = False
        if waypoints.num_waypoints == 0:
            waypoints.flag_manager_requests_waypoints = True
        else:
            if waypoints.type == 'straight_line':
                self.line_manager(waypoints, state)
            elif waypoints.type == 'fillet':
                self.fillet_manager(waypoints, radius, state)
            elif waypoints.type == 'dubins':
                self.dubins_manager(waypoints, radius, state)
            else:
                print('Error in Path Manager: Undefined waypoint type.')

        return self.path

    def line_manager(self, waypoints, state):

        p = np.array([[state.pn], [state.pe], [-state.h]])
        w = waypoints.ned.T              

        if self.new_waypoint_path:
            self.initialize_pointers()
   
        self.halfspace_r = np.array([w[self.ptr_current]]).T
        qi_p = (w[self.ptr_current]-w[self.ptr_previous])/np.linalg.norm(w[self.ptr_current]-w[self.ptr_previous])
        qi = (w[self.ptr_next]-w[self.ptr_current])/np.linalg.norm(w[self.ptr_next]-w[self.ptr_current])
        self.halfspace_n = np.array([(qi_p+qi)/np.linalg.norm(qi_p+qi)]).T        
        self.path.line_origin = np.array([w[self.ptr_current]]).T
        self.path.line_direction = np.array([qi_p]).T        

        if self.inHalfSpace(p):  
            self.increment_pointers(waypoints.num_waypoints)          
            self.path.flag_path_changed = True
            self.path.line_origin = waypoints.ned[self.ptr_previous]
            self.path.line_direction = qi_p

    
    def fillet_manager(self, waypoints, radius, state):        
        p = np.array([[state.pn], [state.pe], [-state.h]])
        w = waypoints.ned.T     
       
        if self.new_waypoint_path:
            self.initialize_pointers()
            self.manager_state = 1

        qi_p = (w[self.ptr_current]-w[self.ptr_previous])/np.linalg.norm(w[self.ptr_current]-w[self.ptr_previous])
        qi = (w[self.ptr_next]-w[self.ptr_current])/np.linalg.norm(w[self.ptr_next]-w[self.ptr_current])
        ro_var = np.arccos(-qi_p.T@qi)

        #straight line
        if self.manager_state == 1:            
            self.path.flag = 'line'
            r = w[self.ptr_current]            
            z = w[self.ptr_current] - radius/np.tan(ro_var/2.0)*qi_p
            self.halfspace_r = np.array([z]).T
            self.halfspace_n = np.array([qi_p]).T
            #set the path for path follower
            self.path.line_origin = np.array([r]).T
            self.path.line_direction = np.array([qi_p]).T

            switch = self.inHalfSpace(p)
            if switch:
                self.manager_state = 2
                self.state_changed = True

        elif self.manager_state == 2:
            self.path.flag_path_changed = self.state_changed
            self.state_changed = False
            qn = qi_p - qi
            qn /= np.linalg.norm(qn)
            c = w[self.ptr_current] - (radius/np.sin(ro_var/2.0))*qn
            lam = np.sign(qi_p.item(0)*qi.item(1) - qi_p.item(1)*qi.item(0))

            self.path.flag = 'orbit'
            self.path.orbit_center = np.array([c]).T
            self.path.orbit_radius = radius            
            if lam > 0:
                self.path.orbit_direction = 'CCW'
            else:
                self.path.orbit_direction = 'CW'

            z = w[self.ptr_current] + (radius/np.tan(ro_var/2))*qi
            self.halfspace_r = np.array([z]).T
            self.halfspace_n = np.array([qi]).T

            if self.inHalfSpace(p):
                self.increment_pointers(waypoints.num_waypoints)
                self.manager_state = 1
                self.state_changed = True


    def dubins_manager(self, waypoints, radius, state):
        p = np.array([[state.pn], [state.pe], [-state.h]])
        w = waypoints.ned.T
        chi = waypoints.course    

        if self.new_waypoint_path:
            self.initialize_pointers()
            self.manager_state = 1
        #dubins parameters
        self.dubins_path.update(w[self.ptr_previous], chi[self.ptr_previous][0], w[self.ptr_current], chi[self.ptr_current][0], radius)

        #Incoming Orbit 
        if self.manager_state == 1:          
            self.path.flag = 'orbit'           
            self.path.orbit_center = self.dubins_path.center_s
            self.path.orbit_direction = self.dubins_path.dir_s
            self.path.orbit_radius = radius
            self.halfspace_n = -self.dubins_path.n1
            self.halfspace_r = self.dubins_path.r1
            if self.flag_path_changed:
                self.path.flag_path_changed = True
                self.flag_path_changed = False
            if self.inHalfSpace(p):
                self.manager_state = 2
                self.flag_path_changed = True  

        #first orbit
        elif self.manager_state == 2:
            self.halfspace_n = self.dubins_path.n1
            self.halfspace_r = self.dubins_path.r1
            if self.flag_path_changed:
                self.path.flag_path_changed = True
                self.flag_path_changed = False
            if self.inHalfSpace(p):
                self.manager_state = 3
                self.flag_path_changed = True  

        #straight line 
        elif self.manager_state == 3:
            self.path.flag = 'line'                      
            self.path.line_origin = self.dubins_path.r2
            self.path.line_direction = self.dubins_path.n1 
            self.halfspace_r = self.dubins_path.r2
            self.halfspace_n = self.dubins_path.n1 
            if self.flag_path_changed:
                self.path.flag_path_changed = True
                self.flag_path_changed = False
            if self.inHalfSpace(p):
                self.manager_state = 4
                self.flag_path_changed = True  

        # Start second orbit 
        elif self.manager_state == 4:            
            self.path.flag = 'orbit'           
            self.path.orbit_center = self.dubins_path.center_e
            self.path.orbit_direction = self.dubins_path.dir_e
            self.path.orbit_radius = radius
            self.halfspace_r = self.dubins_path.r3
            self.halfspace_n = self.dubins_path.n3
            if self.flag_path_changed:
                self.path.flag_path_changed = True
                self.flag_path_changed = False
            if self.inHalfSpace(p):
                self.manager_state = 5
                self.flag_path_changed = True  

        # end second orbit 
        elif self.manager_state == 5:           
            self.path.flag = 'line'
            self.halfspace_r = self.dubins_path.r3
            self.halfspace_n = self.dubins_path.n3
            if self.flag_path_changed:
                self.path.flag_path_changed = True
                self.flag_path_changed = False
            if self.inHalfSpace(p):
                self.manager_state = 1
                self.flag_path_changed = True  
                self.increment_pointers(waypoints.num_waypoints)
                self.dubins_path.update(w[self.ptr_previous], chi[self.ptr_previous][0], w[self.ptr_current], chi[self.ptr_current][0], radius)
            self.path.line_origin = self.dubins_path.r3
            self.path.line_direction = self.dubins_path.n3   


    def initialize_pointers(self):
        self.ptr_previous = 0
        self.ptr_current = 1
        self.ptr_next = 2
        self.new_waypoint_path = False

    def increment_pointers(self, num_waypoints):
        self.ptr_previous = self.ptr_current
        self.ptr_current = self.ptr_next
        self.ptr_next = self.ptr_next + 1

        if self.ptr_next == num_waypoints:
            self.ptr_next = 0

    def inHalfSpace(self, pos):
        self.halfspace_r = self.halfspace_r.reshape((3,1))
        if (pos-self.halfspace_r).T @ self.halfspace_n >= 0:
            return True
        else:
            return False