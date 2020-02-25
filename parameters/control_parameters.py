import sys
sys.path.append('..')
import numpy as np
import chap5.transfer_function_coef as TF
import parameters.aerosonde_parameters as MAV

gravity = MAV.gravity
sigma = 0.05 # Which value should I put in ??
Va0 = MAV.Va0

#----------roll loop-------------
roll_kp = 3.0
roll_kd = 0.2

#----------course loop-------------
course_kp = 1000
course_ki = 0.9955453022851157

#----------sideslip loop-------------
sideslip_ki = 0     # temperate value
sideslip_kp = 1     # temperate value

#----------yaw damper-------------
# yaw_damper_tau_r =
# yaw_damper_kp =

#----------pitch loop-------------
pitch_kp = -4.5
pitch_kd = -1
K_theta_DC = 0.8555133079847909



#----------altitude loop-------------
altitude_kp = 2.5
altitude_ki = 0.04485491131643955
altitude_zone = 10  # equals Wh?

#---------airspeed hold using throttle---------------


airspeed_throttle_kp = 70
airspeed_throttle_ki = 25.592629794342578
