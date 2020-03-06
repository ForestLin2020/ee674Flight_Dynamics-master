"""
observer
    - Beard & McLain, PUP, 2012
    - Last Update:
        3/2/2019 - RWB
"""
import sys
import numpy as np
sys.path.append('..')
import parameters.control_parameters as CTRL
import parameters.simulation_parameters as SIM
import parameters.sensor_parameters as SENSOR
#from tools.rotations import Euler2Rotation
from tools.wrap import wrap
import parameters.aerosonde_parameters as MAV

from message_types.msg_state import msg_state

class observer:
    def __init__(self, ts_control):
        # initialized estimated state message
        self.estimated_state = msg_state()
        # use alpha filters to low pass filter gyros and accels
        self.lpf_gyro_x = alpha_filter(alpha=0.5)
        self.lpf_gyro_y = alpha_filter(alpha=0.5)
        self.lpf_gyro_z = alpha_filter(alpha=0.5)
        self.lpf_accel_x = alpha_filter(alpha=0.5)
        self.lpf_accel_y = alpha_filter(alpha=0.5)
        self.lpf_accel_z = alpha_filter(alpha=0.5)
        # use alpha filters to low pass filter static and differential pressure
        self.lpf_static = alpha_filter(alpha=0.9)
        self.lpf_diff = alpha_filter(alpha=0.5)
        # ekf for phi and theta
        self.attitude_ekf = ekf_attitude()
        # ekf for pn, pe, Vg, chi, wn, we, psi
        self.position_ekf = ekf_position()

    def update(self, measurements):

        # estimates for p, q, r are low pass filter of gyro minus bias estimate
        self.estimated_state.p = self.lpf_gyro_x.update(measurements.gyro_x) - self.estimated_state.bx
        self.estimated_state.q = self.lpf_accel_y.update(measurements.gyro_y) - self.estimated_state.by
        self.estimated_state.r = self.lpf_accel_z.update(measurements.gyro_z) - self.estimated_state.bz

        # sensor model to get altitude and airspeed
        static_pressure = self.lpf_static.update(measurements.static_pressure)
        diff_prssure = self.lpf_diff.update(measurements.diff_pressure)
        self.estimated_state.h = static_pressure / (MAV.rho * MAV.gravity)
        self.estimated_state.Va = np.sqrt( (2*diff_prssure) / MAV.rho )

        # estimate phi and theta with simple ekf
        self.attitude_ekf.update(self.estimated_state, measurements)

        # estimate pn, pe, Vg, chi, wn, we, psi
        self.position_ekf.update(self.estimated_state, measurements)

        # not estimating these
        self.estimated_state.alpha = self.estimated_state.theta
        self.estimated_state.beta = 0.0
        self.estimated_state.bx = 0.0
        self.estimated_state.by = 0.0
        self.estimated_state.bz = 0.0
        return self.estimated_state

class alpha_filter:
    # alpha filter implements a simple low pass filter
    # y[k] = alpha * y[k-1] + (1-alpha) * u[k]
    def __init__(self, alpha=0.5, y0=0.0):
        self.alpha = alpha  # filter parameter
        self.y = y0  # initial condition

    def update(self, u):
        self.y = self.alpha*self.y + (1-self.alpha)*u
        return self.y

class ekf_attitude:
    # implement continous-discrete EKF to estimate roll and pitch angles
    def __init__(self):
        self.Q =
        self.Q_gyro =
        self.R_accel =
        self.N =   # number of prediction step per sample
        self.xhat =  # initial state: phi, theta
        self.P =
        self.Ts = SIM.ts_control/self.N

    def update(self, state, measurement):
        self.propagate_model(state)
        self.measurement_update(state, measurement)
        state.phi = self.xhat.item(0)
        state.theta = self.xhat.item(1)

    # what x use for?
    def f(self, x, state):
        # system dynamics for propagation model: xdot = f(x, u)
        # x = (phi, theta)
        # u = (p, q, r, Va)

        p = state.p
        q = state.q
        r = state.r
        phi = x.item(0)
        theta = x.item(1)

        a1 = p + q*np.sin(phi)*np.tan(theta) + r*np.cos(phi)*np.tan(theta)
        a2 = q*np.cos(phi) - r*np.sin(phi)


        _f = np.array([a1,a2])
        return _f


    def h(self, x, state):
        # measurement model y
        p = state.p
        q = state.q
        r = state.r
        Va = state.Va
        phi =  x.item(0)
        theta =  x.item(1)

        a1 = q*Va*np.sin(theta) + MAV.gravity*np.sin(theta)
        a2 = r*Va*np.cos(theta) - p*Va*np.sin(theta) - MAV.gravity*np.cos(theta)*np.sin(phi)
        a3 = -q*Va*np.cos(theta) - MAV.gravity*np.cos(theta)*np.cos(phi)

        _h = np.array([a1,a2,a3])
        return _h

    def propagate_model(self, state):
        # model propagation
        for i in range(0, self.N):
             # propagate model
            self.xhat = self.xhat + self.Ts*self.f(self.xhat, state)
            # compute Jacobian
            A = jacobian(self.f, self.xhat, state)
            # compute G matrix for gyro noise
            phi = self.xhat.item(0)
            theta = self.xhat.item(1)
            G = np.array([[1.0, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta), 0.0],[0.0, np.cos(phi), -np.sin(phi), 0.0]])
            # update P with continuous time model
            # self.P = self.P + self.Ts * (A @ self.P + self.P @ A.T + self.Q + G @ self.Q_gyro @ G.T)
            # convert to discrete time models
            A_d = np.eye(2) + A*self.Ts + (A@A*self.Ts**2)/2
            G_d =
            # update P with discrete time model
            self.P = A_d @ self.P @ A_d.T + self.Ts**2 * Q

    def measurement_update(self, state, measurement):
        # measurement updates
        threshold = 2.0
        h = self.h(self.xhat, measurement, state)
        C = jacobian(self.h, self.xhat, state)
        y = np.array([measurement.accel_x, measurement.accel_y, measurement.accel_z])
        for i in range(0, 3):
            if np.abs(y[i]-h[i,0]) < threshold:
                Ci =
                L =
                self.P =
                self.xhat =

class ekf_position:
    # implement continous-discrete EKF to estimate pn, pe, chi, Vg
    def __init__(self):
        self.Q =
        self.R =
        self.N =   # number of prediction step per sample
        self.Ts = (SIM.ts_control / self.N)
        self.xhat =
        self.P =
        self.gps_n_old = 9999
        self.gps_e_old = 9999
        self.gps_Vg_old = 9999
        self.gps_course_old = 9999


    def update(self, state, measurement):
        self.propagate_model(state)
        self.measurement_update(state, measurement)
        state.pn = self.xhat.item(0)
        state.pe = self.xhat.item(1)
        state.Vg = self.xhat.item(2)
        state.chi = self.xhat.item(3)
        state.wn = self.xhat.item(4)
        state.we = self.xhat.item(5)
        state.psi = self.xhat.item(6)

    def f(self, x, state):
        # system dynamics for propagation model: xdot = f(x, u)
        pn = x.item(0)
        pe = x.item(1)
        Vg = x.item(2)
        chi = x.item(3)
        wn = x.item(4)
        we = x.item(5)
        psi = x.item(6)
        Va = state.Va
        q = state.q
        r =state.r
        phi = state.phi
        theta = state.theta
        psi_dot = ( (q*np.sin(phi) ) / np.cos(theta)) + ( (r*np.cos(phi)) / np.cos(theta) )
        # p.38
        _f =np.array([Vg*np.cos(chi),
                      Vg*np.sin(chi),
                      ((Va*np.cos(psi)+wn)*(-Va*psi_dot*np.sin(psi)) + (Va*np.sin(psi)+we)*(Va*psi_dot*np.cos(psi)))/Vg,
                      (MAV.gravity*np.tan(phi)*np.cos(chi-psi)),
                      0,
                      0,
                      psi_dot])
        return _f

    def h_gps(self, x, state):
        # measurement model for gps measurements
        _h =
        return _h

    def h_pseudo(self, x, state):
        # measurement model for wind triangale pseudo measurement
        _h =
        return _h

    def propagate_model(self, state):
        # model propagation
        for i in range(0, self.N):
            # propagate model
            self.xhat =
            # compute Jacobian
            A = jacobian(self.f, self.xhat, state)
            # update P with continuous time model
            # self.P = self.P + self.Ts * (A @ self.P + self.P @ A.T + self.Q + G @ self.Q_gyro @ G.T)
            # convert to discrete time models
            A_d =
            # update P with discrete time model
            self.P =

    def measurement_update(self, state, measurement):
        # always update based on wind triangle pseudu measurement
        h = self.h_pseudo(self.xhat, state)
        C = jacobian(self.h_pseudo, self.xhat, state)
        y = np.array([0, 0])
        for i in range(0, 2):
            Ci =
            L =
            self.P =
            self.xhat =

        # only update GPS when one of the signals changes
        if (measurement.gps_n != self.gps_n_old) \
            or (measurement.gps_e != self.gps_e_old) \
            or (measurement.gps_Vg != self.gps_Vg_old) \
            or (measurement.gps_course != self.gps_course_old):

            h = self.h_gps(self.xhat, state)
            C = jacobian(self.h_gps, self.xhat, state)
            y = np.array([measurement.gps_n, measurement.gps_e, measurement.gps_Vg, measurement.gps_course])
            for i in range(0, 4):
                Ci =
                L =
                self.P =
                self.xhat =
            # update stored GPS signals
            self.gps_n_old = measurement.gps_n
            self.gps_e_old = measurement.gps_e
            self.gps_Vg_old = measurement.gps_Vg
            self.gps_course_old = measurement.gps_course

def jacobian(fun, x, state):
    # compute jacobian of fun with respect to x
    f = fun(x, state)
    m = f.shape[0]
    n = x.shape[0]
    eps = 0.01  # deviation
    J = np.zeros((m, n))
    for i in range(0, n):
        x_eps = np.copy(x)
        x_eps[i][0] += eps
        f_eps = fun(x_eps, state)
        df = (f_eps - f) / eps
        J[:, i] = df[:, 0]
    return J