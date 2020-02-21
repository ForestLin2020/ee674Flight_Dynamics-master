
import sys
sys.path.append('..')
import numpy as np
import parameters.aerosonde_parameters as MAV

mass = MAV.mass
gravity = MAV.gravity
rho = MAV.rho
Va = MAV.Va0
Jx = MAV.Jx
Jy = MAV.Jy
Jz = MAV.Jz
Jxz = MAV.Jxz
S = MAV.S_wing
b = MAV.b
c = MAV.c
C_p_0 = MAV.C_p_0
C_p_beta = MAV.C_p_beta
C_p_p = MAV.C_p_p
C_p_r = MAV.C_p_r
C_p_delta_a = MAV.C_p_delta_a
C_p_delta_r = MAV.C_p_delta_r
C_r_0 = MAV.C_r_0
C_r_beta = MAV.C_r_beta
C_r_p = MAV.C_r_p
C_r_r = MAV.C_r_r
C_r_delta_a = MAV.C_r_delta_a
C_r_delta_r = MAV.C_r_delta_r
D = MAV.D_prop
C_Y_0 = MAV.C_Y_0
C_Y_beta = MAV.C_Y_beta
C_Y_p = MAV.C_Y_p
C_Y_r = MAV.C_Y_r
C_Y_delta_a = MAV.C_Y_delta_a
C_Y_delta_r = MAV.C_Y_delta_r
C_ell_0 = MAV.C_ell_0
C_ell_beta = MAV.C_ell_beta
C_ell_p = MAV.C_ell_p
C_ell_r = MAV.C_ell_r
C_ell_delta_a = MAV.C_ell_delta_a
C_ell_delta_r = MAV.C_ell_delta_r
C_n_0 = MAV.C_n_0
C_n_beta = MAV.C_n_beta
C_n_p = MAV.C_n_p
C_n_r = MAV.C_n_r
C_n_delta_a = MAV.C_n_delta_a
C_n_delta_r = MAV.C_n_delta_r
C_L_0 = MAV.C_L_0
C_L_alpha = MAV.C_L_alpha
C_L_q = MAV.C_L_q
C_L_delta_e = MAV.C_L_delta_e
C_D_0 = MAV.C_D_0
C_D_alpha = MAV.C_D_alpha
C_D_p = MAV.C_D_p
C_D_q = MAV.C_D_q
C_D_delta_e = MAV.C_D_delta_e
C_m_0 = MAV.C_m_0
C_m_alpha = MAV.C_m_alpha
C_m_q = MAV.C_m_q
C_m_delta_e = MAV.C_m_delta_e
C_prop = MAV.C_prop
M = MAV.M
alpha0 = MAV.alpha0
epsilon = MAV.epsilon


######################################################################################
                #   Lateral Transfer Function - Roll + Course
######################################################################################
a_phi_1 = -(1/2) * (rho*Va**2*S*b*C_p_p*b)/(2*Va)
a_phi_2 =  (1/2) * (rho*Va**2*S*b*C_p_delta_a)
W_chi = 20      #  outer loop bigger than inner loop more than twice
Vg = 20         # groundspeed ?? depend on me?? is there have any range?
Wn_theta = np.sqrt(3.0*a_phi_2)
W_n_chi = Wn_theta / W_chi
K_i_chi = (W_n_chi**2 * Vg) / gravity # 1.constant >> change Kp ??? 2.how to set up value in simulink

print('a_phi_1 = ', a_phi_1)
print('a_phi_2 = ', a_phi_2)
print('K_i_chi = ', K_i_chi)



######################################################################################
#   Longitudinal Transfer Function - Pitch
######################################################################################
a_theta_1 = -(rho * Va**2 * c * S * C_m_q * c) / (2 * Jy * 2 * Va)
a_theta_2 = -(rho * Va**2 * c * S * C_m_alpha) / (2 * Jy)
a_theta_3 = (rho * Va**2 * c * S * C_m_delta_e) / (2 * Jy)

K_theta_DC = ( (-4.5) * a_theta_3 ) / (a_theta_2 + (-4.5)*a_theta_3)
Wn_theta = np.sqrt(a_theta_2 + (-4.5)*a_theta_3)
Wh = 20     #  1. altitude_zone = 20 ??
Wn_h = Wn_theta / Wh
ki_h = Wn_h**2/(K_theta_DC*25)

print('a_theta_1 = ', a_theta_1)
print('a_theta_2 = ', a_theta_2)
print('a_theta_3 = ', a_theta_3)
print('K_theta_DC = ', K_theta_DC)
print('Wn_theta = ', Wn_theta)
print('Wn_h = ', Wn_h)
print('ki_h = ', ki_h) # 1.constants value >> change kp???   2.how to set up value in simulink