import sys
sys.path.append('..')
import numpy as np
import parameters.aerosonde_parameters as MAV


######################################################################################
                #   MAV
######################################################################################

mass = MAV.mass
gravity = MAV.gravity
rho = MAV.rho
Va = MAV.Va0  # chap4 self._Va = MAV.u0 ???
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
C_L_0 =  MAV.C_L_0
C_L_alpha =  MAV.C_L_alpha
C_L_q =  MAV.C_L_q
C_L_delta_e =  MAV.C_L_delta_e
C_D_0 =  MAV.C_D_0
C_D_alpha =  MAV.C_D_alpha
C_D_p =  MAV.C_D_p
C_D_q =  MAV.C_D_q
C_D_delta_e =  MAV.C_D_delta_e
C_m_0 =  MAV.C_m_0
C_m_alpha =  MAV.C_m_alpha
C_m_q =  MAV.C_m_q
C_m_delta_e =  MAV.C_m_delta_e
C_prop =  MAV.C_prop
M =  MAV.M
alpha0 =  MAV.alpha0
epsilon =  MAV.epsilon

######################################################################################
                #   Lateral Transfer Function - Roll + Course
######################################################################################
a_phi_1 = -(1/2) * (rho*Va**2*S*b*C_p_p*b)/(2*Va)
a_phi_2 =  (1/2) * (rho*Va**2*S*b*C_p_delta_a)

######################################################################################
                #   Lateral Transfer Function - Sideslip
######################################################################################
a_beta_1 = -(rho*Va*S*C_Y_beta) / (2*mass)
a_beta_2 =  (rho*Va*S*C_Y_delta_r) / (2*mass)

######################################################################################
                #   Longitudinal Transfer Function - Pitch
######################################################################################
a_theta_1 = -(rho*Va**2*c*S)*(C_m_q)*(c) / (2*Jy) / (2*Va)
a_theta_2 = -(rho*Va**2*c*S)*(C_m_alpha) / (2*Jy)
a_theta_3 =  (rho*Va**2*c*S)*(C_m_delta_e) / (2*Jy)

######################################################################################
                #   Longitudinal Transfer Function - Airspeed
######################################################################################
# a_V1 =
# a_V2 =
# a_V3 =