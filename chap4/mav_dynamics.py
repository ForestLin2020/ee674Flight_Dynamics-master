"""
mav_dynamics
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
"""
import sys
sys.path.append('..')
import numpy as np
import math

# load message types
from message_types.msg_state import msg_state
from chap4.wind_simulation import wind_simulation
import parameters.aerosonde_parameters as MAV
import parameters.simulation_parameters as SIM
from tools.tools import Quaternion2Euler, TransformationFormInertialToBody




class mav_dynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.pn0],  # (0)
                               [MAV.pe0],   # (1)
                               [MAV.pd0],   # (2)
                               [MAV.u0],    # (3)
                               [MAV.v0],    # (4)
                               [MAV.w0],    # (5)
                               [MAV.e0],    # (6)
                               [MAV.e1],    # (7)
                               [MAV.e2],    # (8)
                               [MAV.e3],    # (9)
                               [MAV.p0],    # (10)
                               [MAV.q0],    # (11)
                               [MAV.r0]])   # (12)
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        self._update_velocity_data()
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.Va0
        self._alpha = 0
        self._beta = 0
        # initialize true_state message
        self.msg_true_state = msg_state()

    ###################################
    # public functions
    def update_state(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)

        # Integrate ODE using Runge-Kutta RK4 algorithm
        time_step = self._ts_simulation
        k1 = self._derivatives(self._state, forces_moments)
        k2 = self._derivatives(self._state + time_step/2.*k1, forces_moments)
        k3 = self._derivatives(self._state + time_step/2.*k2, forces_moments)
        k4 = self._derivatives(self._state + time_step*k3, forces_moments)
        self._state += time_step/6 * (k1 + 2*k2 + 2*k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = np.sqrt(e0**2+e1**2+e2**2+e3**2)

        self._state[6][0] = self._state.item(6)/normE
        self._state[7][0] = self._state.item(7)/normE
        self._state[8][0] = self._state.item(8)/normE
        self._state[9][0] = self._state.item(9)/normE

        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)

        # update the message class for the true state
        self._update_msg_true_state()

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """
        # extract the states
        pn = state.item(0)
        pe = state.item(1)
        pd = state.item(2)
        u = state.item(3)
        v = state.item(4)
        w = state.item(5)
        e0 = state.item(6)
        e1 = state.item(7)
        e2 = state.item(8)
        e3 = state.item(9)
        p = state.item(10)
        q = state.item(11)
        r = state.item(12)

        #   extract forces/moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)

        l = forces_moments.item(3)
        m = forces_moments.item(4)
        n = forces_moments.item(5)

        # position kinematics
        pn_dot = (e1**2 + e0**2 - e2**2 - e3**2) * u + 2*(e1*e2 - e3*e0) * v + 2 * (e1*e3 + e2*e0) * w
        pe_dot = 2 * (e1*e2 + e3*e0) * u + (e2**2 + e0**2 - e1**2 - e3**2) * v + 2 * (e2*e3 - e1*e0) * w
        pd_dot = 2 * (e1*e3 - e2*e0) * u + 2 * (e2*e3 + e1*e0) * v + (e3**2 + e0**2 - e1**2 - e2**2) * w

        # position dynamics
        u_dot = (r*v - q*w) + (fx / MAV.mass)
        v_dot = (p*w - r*u) + (fy / MAV.mass)
        w_dot = (q*u - p*v) + (fz / MAV.mass)

        # rotational kinematics
        e0_dot = (-p * e1 - q * e2 - r * e3) *0.5
        e1_dot = ( p * e0 + r * e2 - q * e3) *0.5
        e2_dot = ( q * e0 - r * e1 + p * e3) *0.5
        e3_dot = ( r * e0 + q * e1 - p * e2) *0.5

        L  = MAV.Jx * MAV.Jz - MAV.Jxz**2
        L1 = MAV.Jxz * (MAV.Jx - MAV.Jy + MAV.Jz) / L
        L2 = (MAV.Jz * (MAV.Jz - MAV.Jy) + MAV.Jxz**2) / L
        L3 = MAV.Jz / L
        L4 = MAV.Jxz / L
        L5 = (MAV.Jz - MAV.Jx) / MAV.Jy
        L6 = MAV.Jxz / MAV.Jy
        L7 = ((MAV.Jx - MAV.Jy) * MAV.Jx + MAV.Jxz**2) / L
        L8 = MAV.Jx / L

        # rotatonal dynamics
        p_dot = (L1 * p * q - L2 * q * r) + (L3 * l + L4 * n)
        q_dot = (L5 * p * r - L6 * (p ** 2 - r ** 2)) + (m / MAV.Jy)
        r_dot = (L7 * p * q - L1 * q * r) + (L4 * l + L8 * n)

        # collect the derivative of the states
        x_dot = np.array([[pn_dot, pe_dot, pd_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6,1))):
        # self.wind = wind_simulation(SIM.ts_simulation) # how can I know there are six vector come here?
        wn_s = wind.item(0)
        we_s = wind.item(1)
        wd_s = wind.item(2)

        wind_ss = [wn_s, we_s, wd_s]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        angles = phi, theta, psi
        wind_ss_b = TransformationFormInertialToBody(wind_ss, angles)

        wn_b = wind_ss_b[0]
        we_b = wind_ss_b[1]
        wd_b = wind_ss_b[2]

        u_w = wn_b + wind.item(3)
        v_w = we_b + wind.item(4)
        w_w = wd_b + wind.item(5)

        u = self._state.item(3)
        v = self._state.item(4)
        w = self._state.item(5)
        
        u_r = u - u_w
        v_r = v - v_w
        w_r = w - w_w

        # compute airspeed
        self._Va = np.sqrt(u_r**2 + v_r**2 + w_r**2)
        # compute angle of attack
        self._alpha = np.arctan(w_r/u_r)
        # compute sideslip angle
        self._beta = np.arcsin(v_r/self._Va)

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_e, delta_t, delta_a, delta_r)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """

        # F(fx,fy,fz) = Fg(x,y,z) + Fa(x,y,z) + Fp(x,y,z)

        phi, theta, psi = Quaternion2Euler(self._state[6:10])

        sin = np.sin
        cos = np.cos
        mass = MAV.mass
        gravity = MAV.gravity
        rho = MAV.rho
        Va = self._Va
        S = MAV.S_wing
        alpha = self._alpha
        beta = self._beta
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)
        D = MAV.D_prop

        Vin = MAV.V_max * delta[1]
        a = (rho * np.power(D,5) * MAV.C_Q0) / (4*np.pi**2)
        b = (rho * np.power(D,4) * MAV.C_Q1 * Va)/(2*np.pi) + (MAV.KQ**2)/(MAV.R_motor)
        c = (rho * np.power(D,3) * MAV.C_Q2 * Va**2) - (MAV.KQ*Vin/MAV.R_motor) + (MAV.KQ*MAV.i0)

        Omaga_p = (-b + math.sqrt(b**2 - (4*a*c))) / (2*a)
        J = (2*np.pi*Va) / (Omaga_p*D)
        C_T = MAV.C_T2*J**2 + MAV.C_T1*J + MAV.C_T0
        C_Q = MAV.C_Q2*J**2 + MAV.C_Q1*J + MAV.C_Q0



        Tp = (rho * np.power(D,4) * Omaga_p**2 * C_T) / (4*np.pi**2)
        Qp = (rho * np.power(D,5) * Omaga_p**2 * C_Q) / (4*np.pi**2)
        # Qp = (rho*np.power(D,5)*MAV.C_Q0/(4*np.pi**2))*Omaga_p**2 + (rho*np.power(D,4)*MAV.C_Q1*Va/(2*np.pi))*Omaga_p + rho*np.power(D,3)*MAV.C_Q2*Va**2

        C_Y_0, C_Y_beta, C_Y_p, C_Y_r, C_Y_delta_a, C_Y_delta_r = MAV.C_Y_0 ,MAV.C_Y_beta, MAV.C_Y_p, MAV.C_Y_r, MAV.C_Y_delta_a, MAV.C_Y_delta_r

        C_D_alpha = MAV.C_D_0 + MAV.C_D_alpha*alpha
        C_L_alpha = MAV.C_L_0 + MAV.C_L_alpha*alpha

        C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha)
        C_X_q_alpha = -MAV.C_D_q*cos(alpha) + MAV.C_L_q*sin(alpha)
        C_X_delta_e_alpha = -MAV.C_D_delta_e*cos(alpha) + MAV.C_L_delta_e*sin(alpha)
        C_Z_alpha = -C_D_alpha*sin(alpha) - C_L_alpha*cos(alpha)
        C_Z_q_alpha = -MAV.C_D_q*sin(alpha) - MAV.C_L_q*cos(alpha)
        C_Z_delta_e_alpha = -MAV.C_D_delta_e*sin(alpha) + MAV.C_L_delta_e*cos(alpha)


        fx = (-mass*gravity*sin(theta)) + Tp + (1/2)*rho*Va**2*S*(C_X_alpha + (C_X_q_alpha*MAV.c*q)/(2*Va)) + (1/2)*rho*Va**2*S*(C_X_delta_e_alpha*delta[0])
        fy = (mass*gravity*cos(theta)*sin(phi)) + (1/2)*rho*Va**2*S*(C_Y_0 + C_Y_beta*beta + (C_Y_p*MAV.b*p)/(2*Va) + (C_Y_r*MAV.b*r)/(2*Va)) + (1/2)*rho*Va**2*S*(C_Y_delta_a*delta[2] + C_Y_delta_r*delta[3])
        fz = (mass*gravity*cos(theta)*cos(phi)) + (1/2)*rho*Va**2*S*(C_Z_alpha + (C_Z_q_alpha*MAV.c*q)/(2*Va)) + (1/2)*rho*Va**2*S*(C_Z_delta_e_alpha*delta[0])

        C_ell_0, C_ell_beta, C_ell_p, C_ell_r, C_ell_delta_a, C_ell_delta_r = MAV.C_ell_0, MAV.C_ell_beta, MAV.C_ell_p, MAV.C_ell_r, MAV.C_ell_delta_a, MAV.C_ell_delta_r
        C_n_0, C_n_beta, C_n_p, C_n_r, C_n_delta_a, C_n_delta_r = MAV.C_n_0, MAV.C_n_beta, MAV.C_n_p, MAV.C_n_r, MAV.C_n_delta_a, MAV.C_n_delta_r

        l = (1/2)*rho*Va**2*S*(MAV.b*( C_ell_0 + C_ell_beta*beta + ((C_ell_p*MAV.b*p)/(2*Va)) + ((C_ell_r*MAV.b*r)/(2*Va)))) + (1/2)*rho*Va**2*S*(MAV.b*(C_ell_delta_a*delta[2] + C_ell_delta_r*delta[3]))
        m = (1/2)*rho*Va**2*S*(MAV.c*( MAV.C_m_0 + MAV.C_m_alpha*alpha + (MAV.C_m_q*MAV.c*q)/(2*Va))) + (1/2)*rho*Va**2*S*(MAV.c*(MAV.C_m_delta_e*delta[0]))
        n = (1/2)*rho*Va**2*S*(MAV.b*( C_n_0 + C_n_beta*beta + (C_n_p*MAV.b*p)/(2*Va) + (C_n_r*MAV.b*r)/(2*Va))) + (1/2)*rho*Va**2*S*(MAV.b*( C_n_delta_a*delta[2] + C_n_delta_r*delta[3]))

        Mx = l + Qp
        My = m
        Mz = n

        # fx = np.array([0])
        # fy = np.array([0])
        # fz = np.array([0])
        # Mx = np.array([0])
        # My = np.array([0])
        # Mz = np.array([0])


        self._forces[0] = fx
        self._forces[1] = fy
        self._forces[2] = fz

        return np.array([[fx, fy, fz, Mx, My, Mz]]).T



    def _update_msg_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        self.msg_true_state.pn = self._state.item(0)
        self.msg_true_state.pe = self._state.item(1)
        self.msg_true_state.h = -self._state.item(2)
        self.msg_true_state.Va = self._Va
        self.msg_true_state.alpha = self._alpha
        self.msg_true_state.beta = self._beta
        self.msg_true_state.phi = phi
        self.msg_true_state.theta = theta
        self.msg_true_state.psi = psi

        u = self._state.item(3)
        v = self._state.item(4)
        w = self._state.item(5)
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)


        # position kinematics
        pn_dot = (e1 ** 2 + e0 ** 2 - e2 ** 2 - e3 ** 2) * u + 2 * (e1 * e2 - e3 * e0) * v + 2 * (e1 * e3 + e2 * e0) * w
        pe_dot = 2 * (e1 * e2 + e3 * e0) * u + (e2 ** 2 + e0 ** 2 - e1 ** 2 - e3 ** 2) * v + 2 * (e2 * e3 - e1 * e0) * w
        pd_dot = 2 * (e1 * e3 - e2 * e0) * u + 2 * (e2 * e3 + e1 * e0) * v + (e3 ** 2 + e0 ** 2 - e1 ** 2 - e2 ** 2) * w


        self.msg_true_state.Vg = np.sqrt((self._state[3][0])**2+(self._state[4][0])**2+(self._state[5][0])**2)
        self.msg_true_state.gamma = np.arcsin(np.sqrt(pn_dot**2+pe_dot**2)/np.sqrt(pn_dot**2+pe_dot**2 + pd_dot**2))
        self.msg_true_state.chi = np.arctan(pe_dot/pn_dot)

        self.msg_true_state.p = self._state.item(10)
        self.msg_true_state.q = self._state.item(11)
        self.msg_true_state.r = self._state.item(12)
        self.msg_true_state.wn = self._wind.item(0)
        self.msg_true_state.we = self._wind.item(1)
