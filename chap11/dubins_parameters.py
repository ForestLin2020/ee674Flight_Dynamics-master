# dubins_parameters
#   - Dubins parameters that define path between two configurations
#
# mavsim_matlab 
#     - Beard & McLain, PUP, 2012
#     - Update history:  
#         3/26/2019 - RWB

import numpy as np
import sys
sys.path.append('..')

cos = np.cos
sin = np.sin

class dubins_parameters:
    def __init__(self):
        self.p_s = np.inf*np.ones((3,1))  # the start position in re^3
        self.chi_s = np.inf  # the start course angle
        self.p_e = np.inf*np.ones((3,1))  # the end position in re^3
        self.chi_e = np.inf  # the end course angle
        self.radius = np.inf  # turn radius
        self.length = np.inf  # length of the Dubins path
        self.center_s = np.inf*np.ones((3,1))  # center of the start circle
        self.dir_s = np.inf  # direction of the start circle
        self.center_e = np.inf*np.ones((3,1))  # center of the end circle
        self.dir_e = np.inf  # direction of the end circle
        self.r1 = np.inf*np.ones((3,1))  # vector in re^3 defining half plane H1
        self.r2 = np.inf*np.ones((3,1))  # vector in re^3 defining position of half plane H2
        self.r3 = np.inf*np.ones((3,1))  # vector in re^3 defining position of half plane H3
        self.n1 = np.inf*np.ones((3,1))  # unit vector in re^3 along straight line path
        self.n3 = np.inf*np.ones((3,1))  # unit vector defining direction of half plane H3

    def update(self, ps, chis, pe, chie, R):
        ell = np.linalg.norm(ps - pe)
        if ell < 2 * R:
            print('Error in Dubins Parameters: The distance between nodes must be larger than 2R.')
        else:
            self.p_s = ps
            self.chi_s = chis
            self.p_e = pe
            self.chi_e = chie
            self.radius = R
            e1 = np.array([[1,0,0]]).T


            #(n,e,d)
            crs = ps + R * rotz(np.pi/2) @ np.array([[cos(chis),sin(chis),0]]).T        # np.array([]) and np.array([[]])  >> [] and [[]]
            cls = ps + R * rotz(-np.pi/2) @ np.array([[cos(chis),sin(chis),0]]).T
            cre = pe + R * rotz(np.pi/2) @ np.array([[cos(chie),sin(chie),0]]).T
            cle = pe + R * rotz(-np.pi/2) @ np.array([[cos(chie),sin(chie),0]]).T


            vartheta = np.arctan2(cre.item(1)-crs.item(1),cre.item(0)-crs.item(0))
            L1_1 = np.linalg.norm(crs-cre)
            L1_2 = R * mod(2*np.pi + mod(vartheta-0.5*np.pi) - mod(chis-0.5*np.pi))
            L1_3 = R * mod(2*np.pi + mod(chie-0.5*np.pi) - mod(vartheta-0.5*np.pi))
            L1 = L1_1 + L1_2 + L1_3

            ell = np.linalg.norm(cle-crs)
            vartheta = np.arctan2(cle.item(1)-crs.item(1),cle.item(0)-crs.item(0))
            vartheta2 = vartheta - 0.5*np.pi + np.arcsin(2*R/ell)
            L2_1 = np.sqrt(ell**2 - 4*R**2)
            L2_2 = R * mod(2*np.pi + mod(vartheta2) - mod(chis-0.5*np.pi))
            L2_3 = R * mod(2*np.pi + mod(vartheta2+np.pi) + mod(chie+0.5*np.pi))
            L2 = L2_1 + L2_2 + L2_3

            ell = np.linalg.norm(cre-cls)
            vartheta = np.arctan2(cre.item(1)-cls.item(1),cre.item(0)-cls.item(0))
            vartheta2 = np.arccos(2*R/ell)
            L3_1 =np.sqrt(ell**2 - 4*R**2)
            L3_2 = R *mod(2*np.pi + mod(chis+0.5*np.pi) - mod(vartheta+vartheta2))
            L3_3 = R *mod(2*np.pi + mod(chie-0.5*np.pi) - mod(vartheta+vartheta2-np.pi))
            L3 = L3_1 + L3_2 + L3_3

            vartheta = np.arctan2(cle.item(1)-cls.item(1),cle.item(0)-cls.item(0))
            L4_1 = np.linalg.norm(cls-cle)
            L4_2 = R * mod(2*np.pi + mod(chis+0.5*np.pi) - mod(vartheta+0.5*np.pi))
            L4_3 = R * mod(2*np.pi + mod(vartheta+0.5*np.pi) - mod(chie+0.5*np.pi))
            L4 = L4_1 + L4_2 + L4_3

            # Length = (L1, L2, L3, L4)
            L = np.array([L1, L2, L3, L4])
            self.length = np.amin(L)    # what is np.amin() and np.argmin()??????

            if np.argmin(L) == 0:
                self.center_s = crs
                self.dir_s = 1
                self.center_e = cre
                self.dir_e = 1
                # q1
                self.n1 = (self.center_e - self.center_s) / (np.linalg.norm(self.center_e - self.center_s))
                # z1
                self.r1 = self.center_s + R * rotz(-0.5*np.pi) @ self.n1
                # z2
                self.r2 = self.center_e + R * rotz(-0.5*np.pi) @ self.n1

            elif np.argmin(L) == 1:
                self.center_s = crs
                self.dir_s = 1
                self.center_e = cle
                self.dir_e = -1
                ell = np.linalg.norm(self.center_e - self.center_s)
                vartheta = np.arctan2(cle.item(1) - crs.item(1), cle.item(0) - crs.item(0))
                vartheta2 = vartheta - 0.5*np.pi + np.arcsin(2*R/ell)
                # q1
                self.n1 = rotz(vartheta2 + 0.5*np.pi) @ e1
                # z1
                self.r1 = self.center_s + R * rotz(vartheta2) @ e1
                # z2
                self.r2 = self.center_e + R * rotz(vartheta2 + np.pi) @ e1

            elif np.argmin(L) == 2:
                self.center_s = cls
                self.dir_s = -1
                self.center_e = cre
                self.dir_e = 1
                ell = np.linalg.norm(self.center_e - self.center_s)
                vartheta = np.arctan2(cre.item(1) - cls.item(1), cre.item(0) - cls.item(0))
                vartheta2 = np.arccos(2 * R / ell)
                # q1
                self.n1 = rotz(vartheta + vartheta2 - 0.5*np.pi) @ e1
                # z1
                self.r1 = self.center_s + R * rotz(vartheta + vartheta2) @ e1
                # z2
                self.r2 = self.center_e + R * rotz(vartheta + vartheta2 - np.pi) @ e1


                vartheta = np.arctan2((cre.item(1) - cls.item(1)),
                                      (cre.item(0) - cls.item(0)))


            elif np.argmin(L) == 3:
                self.center_s = cls
                self.dir_s = -1
                self.center_e = cle
                self.dir_e = -1
                # q1
                self.n1 = (self.center_e - self.center_s) / (np.linalg.norm(self.center_e - self.center_s))
                # z1
                self.r1 = self.center_s + R * rotz(0.5 * np.pi) @ self.n1
                # z2
                self.r2 = self.center_e + R * rotz(0.5 * np.pi) @ self.n1

            # z3
            self.r3 = pe
            # q3
            self.n3 =rotz(chie) @ e1


def rotz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])


def mod(x):
    # make x between 0 and 2*pi
    while x < 0:
        x += 2*np.pi
    while x > 2*np.pi:
        x -= 2*np.pi
    return x


