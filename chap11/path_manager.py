import numpy as np
import sys
sys.path.append('..')
from chap11.dubins_parameters import dubins_parameters
from message_types.msg_path import msg_path

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
        self.num_waypoints = 0
        self.halfspace_n = np.inf * np.ones((3,1))
        self.halfspace_r = np.inf * np.ones((3,1))
        # state of the manager state machine
        self.manager_state = 1
        # dubins path parameters
        self.dubins_path = dubins_parameters()

    def update(self, waypoints, radius, state):
        # this flag is set for one time step to signal a redraw in the viewer
        if self.flag_need_new_waypoints == True:
            self.flag_need_new_waypoints = False
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

        p = np.array([[state.pn, state.pe, state.h]]).T

        wim1 = np.array([[waypoints.ned[0][self.ptr_previous],waypoints.ned[1][self.ptr_previous],waypoints.ned[2][self.ptr_previous]]]).T
        wi = np.array([[waypoints.ned[0][self.ptr_current],waypoints.ned[1][self.ptr_current],waypoints.ned[2][self.ptr_current]]]).T
        wip1 = np.array([[waypoints.ned[0][self.ptr_next],waypoints.ned[1][self.ptr_next],waypoints.ned[2][self.ptr_next]]]).T

        # return(r)
        self.path.line_origin = wim1

        qim1 = (wi - wim1) / np.linalg.norm(wi - wim1)
        qi = (wip1 - wi) / np.linalg.norm(wip1 - wi)

        # check the Halfspace
        # send ri and ni into function inHalfSpace() check location
        # H(self.halfspace_r, self.halfspace_n)
        self.halfspace_r = wi
        self.halfspace_n = (qim1 + qi)/ np.linalg.norm(qim1 + qi)
        if self.inHalfSpace(p):
            self.increment_pointers()

        # return(q)
        self.path.line_direction = qim1
        ### Return
        # self.path.airspeed = waypoints.airspeed.item(self.ptr_previous) # why???

    def fillet_manager(self, waypoints, radius, state):

        p = np.array([[state.pn, state.pe, state.h]]).T

        R = radius

        wim1 = np.array([[waypoints.ned[0][self.ptr_previous], waypoints.ned[1][self.ptr_previous],
                          waypoints.ned[2][self.ptr_previous]]]).T
        wi = np.array([[waypoints.ned[0][self.ptr_current], waypoints.ned[1][self.ptr_current],
                        waypoints.ned[2][self.ptr_current]]]).T
        wip1 = np.array(
            [[waypoints.ned[0][self.ptr_next], waypoints.ned[1][self.ptr_next], waypoints.ned[2][self.ptr_next]]]).T

        qim1 = (wi - wim1) / np.linalg.norm(wi - wim1)
        qi = (wip1 - wi) / np.linalg.norm(wip1 - wi)

        varrho = np.arccos(-qim1.T @ qi)

        if self.manager_state == 1:
            # return flag, r, q
            self.path.flag = 'line'
            self.path.line_origin = wim1
            self.path.line_direction = qim1

            # check the Halfspace
            # send ri and ni into function inHalfSpace() check location
            # H(self.halfspace_r, self.halfspace_n)
            z = wi - (R/np.tan(varrho/2))*qim1
            self.halfspace_r = z
            self.halfspace_n = qim1
            if self.inHalfSpace(p):
                self.manager_state = 2

        elif self.manager_state == 2:
            # return flag, c, rho, lambda
            self.path.flag = 'orbit'
            self.path.orbit_center = wi + (R/np.sin(varrho/2)) * ((qim1-qi)/np.linalg.norm(qim1-qi))
            self.path.orbit_radius = R
            self.path.orbit_direction = np.sign(qim1.item(0)*qi.item(1) - qim1.item(1)*qi.item(0))

            # check the Halfspace
            # send ri and ni into function inHalfSpace() check location
            # H(self.halfspace_r, self.halfspace_n)
            z = wi + (R / np.tan(varrho / 2)) * qi
            self.halfspace_r = z
            self.halfspace_n = qi
            if self.inHalfSpace(p):
                self.increment_pointers()
                self.manager_state = 1

        ### Return
        # self.path.airspeed = waypoints.airspeed.item(self.ptr_previous) # why???


    def dubins_manager(self, waypoints, radius, state):
        p = np.array([[state.pn, state.pe, state.h]]).T

        R = radius

        # w = [n,e,d]
        wim1 = np.array([[waypoints.ned[0][self.ptr_previous], waypoints.ned[1][self.ptr_previous],
                          waypoints.ned[2][self.ptr_previous]]]).T
        wi = np.array([[waypoints.ned[0][self.ptr_current], waypoints.ned[1][self.ptr_current],
                        waypoints.ned[2][self.ptr_current]]]).T
        chim1 = waypoints.course[0][self.ptr_previous]
        chii = waypoints.course[0][self.ptr_current]

        self.dubins_path.update(wim1, chim1, wi, chii, R)
        # self.p_s
        # self.chi_s
        # self.p_e
        # self.chi_e
        # self.radius
        # self.length
        # self.center_s
        # self.dir_s
        # self.center_e
        # self.dir_e
        # self.r1 (z1)
        # self.r2 (z2)
        # self.r3 (z3)
        # self.n1 (q1)
        # self.n3 (q3)


        if self.manager_state == 1:
            # return
            self.path.flag = 'orbit'
            self.path.orbit_center = self.dubins_path.center_s
            self.path.orbit_radius = R
            self.path.orbit_direction = self.dubins_path.dir_s

            # check the Halfspace
            # send ri and ni into function inHalfSpace() check location
            # H(self.halfspace_r, self.halfspace_n)
            self.halfspace_r = self.dubins_path.r1
            self.halfspace_n = -self.dubins_path.n1
            if self.inHalfSpace(p):
                self.manager_state = 2

        elif self.manager_state == 2:
            # check the Halfspace
            # send ri and ni into function inHalfSpace() check location
            # H(self.halfspace_r, self.halfspace_n)
            self.halfspace_r = self.dubins_path.r1
            self.halfspace_n = self.dubins_path.n1
            if self.inHalfSpace(p):
                self.manager_state = 3

        elif self.manager_state == 3:
            # return
            self.path.flag = 'line'
            self.path.line_origin = self.dubins_path.r1
            self.path.line_direction = self.dubins_path.n1
            # check the Halfspace
            # send r and n into function inHalfSpace() check location
            # H(self.halfspace_r, self.halfspace_n)
            self.halfspace_r = self.dubins_path.r2
            self.halfspace_n = self.dubins_path.n1
            if self.inHalfSpace(p):
                self.manager_state = 4

        elif self.manager_state == 4:
            # return
            self.path.flag = 'orbit'
            self.path.orbit_center = self.dubins_path.center_e
            self.path.orbit_radius = self.dubins_path.radius
            self.path.orbit_direction = self.dubins_path.dir_e
            # check the Halfspace
            # send r and n into function inHalfSpace() check location
            # H(self.halfspace_r, self.halfspace_n)
            self.halfspace_r = self.dubins_path.r3
            self.halfspace_n = -self.dubins_path.n3
            if self.inHalfSpace(p):
                self.manager_state = 5

        elif self.manager_state == 5:
            # check the Halfspace
            # send r and n into function inHalfSpace() check location
            # H(self.halfspace_r, self.halfspace_n)
            self.halfspace_r = self.dubins_path.r3
            self.halfspace_n = self.dubins_path.n3
            if self.inHalfSpace(p):
                self.manager_state = 1
                self.increment_pointers()
                self.dubins_path.update(wim1, chim1, wi, chii, R)

    # def initialize_pointers(self):

    def increment_pointers(self):
        self.ptr_previous = self.ptr_previous + 1
        self.ptr_current = self.ptr_current + 1
        self.ptr_next = self.ptr_next + 1

    def inHalfSpace(self, pos):
        if (pos-self.halfspace_r).T @ self.halfspace_n >= 0:
            return True
        else:
            return False
