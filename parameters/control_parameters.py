import sys
sys.path.append('..')
import numpy as np
import chap5.transfer_function_coef as TF
import parameters.aerosonde_parameters as MAV

gravity = MAV.gravity
sigma = 0.05 # why??
Va0 = MAV.Va0

#----------roll loop-------------
roll_kp = 3.0
roll_kd = 0.2

#----------course loop-------------
course_kp = 1000
course_ki = 0.9955453022851157

#----------sideslip loop-------------
# sideslip_ki =
# sideslip_kp =

#----------yaw damper-------------
# yaw_damper_tau_r =
# yaw_damper_kp =

#----------pitch loop-------------
pitch_kp = -4.5
pitch_kd = -1
K_theta_DC = 0.8555133079847909



#----------altitude loop-------------
altitude_kp = 40
altitude_ki = 0.011213727829109888
altitude_zone = 10 # just pick a number?

#---------airspeed hold using throttle---------------
airspeed_throttle_kp =
airspeed_throttle_ki =
