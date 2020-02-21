"""
compute_trim 
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        2/5/2019 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
import math
from scipy.optimize import minimize
import parameters.aerosonde_parameters as MAV
from tools.tools import Euler2Quaternion

def compute_trim(mav, Va, gamma):
    # define initial state and input
    state0 = np.array([[MAV.pn0],  # (0)
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
    delta0 = np.array([[0],   # (e)
                       [0],   # (t)
                       [0],   # (a)
                       [0]])  # (r)
    x0 = np.concatenate((state0, delta0), axis=0)
    # define equality constraints
    cons = ({'type': 'eq',
             'fun': lambda x: np.array([
                                x[3]**2 + x[4]**2 + x[5]**2 - Va**2,  # magnitude of velocity vector is Va
                                x[4],  # v=0, force side velocity to be zero
                                x[6]**2 + x[7]**2 + x[8]**2 + x[9]**2 - 1.,  # force quaternion to be unit length
                                x[7], # e1=0  - forcing e1=e3=0 ensures zero roll and zero yaw in trim
                                x[9], # e3=0
                                x[10], # p=0  - angular rates should all be zero
                                x[11], # q=0
                                x[12], # r=0
                                ]),
             'jac': lambda x: np.array([
                                [0., 0., 0., 2*x[3], 2*x[4], 2*x[5], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 2*x[6], 2*x[7], 2*x[8], 2*x[9], 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
                                ])
             })
    # solve the minimization problem to find the trim states and inputs
    res = minimize(trim_objective, x0, method='SLSQP', args = (mav, Va, gamma),
                   constraints=cons, options={'ftol': 1e-10, 'disp': True})
    # extract trim state and input and return
    trim_state = np.array([res.x[0:13]]).T
    trim_input = np.array([res.x[13:17]]).T
    return trim_state, trim_input

# objective function to be minimized
def trim_objective(x, mav, Va, gamma):
    state = x[0:13]
    delta = x[13:17]
    xdot = np.array([
                     [0],   # (0)
                     [0],   # (1)
                     [-Va * np.sin(gamma)],   # (2)
                     [0],    # (3)
                     [0],    # (4)
                     [0],    # (5)
                     [0],    # (6)
                     [0],    # (7)
                     [0],    # (8)
                     [0],    # (9)
                     [0],    # (10)
                     [0],    # (11)
                     [0]])   # (12)
    mav._state = state
    mav._update_velocity_data()
    forces_moments = mav._forces_moments(delta)
    f = mav._derivatives(state, forces_moments)
    tmp = xdot[2:13] - f[2:13]
    # tmp = xdot - f

    J = np.linalg.norm(tmp)

    return J

