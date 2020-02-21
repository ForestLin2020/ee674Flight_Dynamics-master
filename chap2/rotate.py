import numpy as np 
import math as m

# c_phi = np.cos(phi)
# s_phi = np.sin(phi)
# c_theta = np.cos(theta)
# s_theta = np.sin(theta)
psi = m.pi 
c_psi = np.cos(psi)
s_psi = np.sin(psi)

R_yaw = np.array([[c_psi, s_psi, 0],
                  [-s_psi, c_psi, 0],
                  [0, 0, 1]])


points = np.array([
                    [-4,  0, 0], # point 0
                    [-2,  2, 2], # point 1
                    [-2, -2, 2], # point 2
                    [-2, -2,-2], # point 3
                    [-2,  2,-2], # point 4
                    [10,  0, 0], # point 5
                    [ 0,  7, 0], # point 6
                    [ 4,  7, 0], # point 7
                    [ 4, -7, 0], # point 8
                    [ 0, -7, 0], # point 9
                    [ 7,  3, 0], # point 10
                    [10,  3, 0], # point 11
                    [10, -3, 0], # point 12
                    [ 7, -3, 0], # point 13
                    [ 7,  0, 0], # point 14
                    [10,  0,-3], # point 15
                  ])




print(R_yaw.shape)
print(points. shape)

new_points = points @ R_yaw

print(new_points)
