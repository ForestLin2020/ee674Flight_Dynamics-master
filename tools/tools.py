
import math
import numpy as np
cos = math.cos
sin = math.sin


def Quaternion2Euler(e):
	e0 = e[0]
	e1 = e[1]
	e2 = e[2]
	e3 = e[3]

	phi = np.arctan2(2*(e0*e1+e2*e3),(e0**2+e3**2-e1**2-e2**2))
	theta = np.arcsin(2*(e0*e2-e1*e3))
	psi = np.arctan2(2*(e0*e3+e1*e2),(e0**2+e1**2-e2**2-e3**2))
	angles = [phi, theta, psi]
	return angles



def Euler2Quaternion(angles):	

	phi = angles[0]
	theta = angles[1]
	psi = angles[2]

	e0 = cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2)
	e1 = cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2)
	e2 = cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2)
	e3 = sin(psi/2)*cos(theta/2)*cos(phi/2) - cos(psi/2)*sin(theta/2)*sin(phi/2)

	e = [e0, e1, e2, e3]
	return e


def TransformationFormInertialToBody(wind_s, angles):
	wn_s = wind_s[0]
	we_s = wind_s[1]
	wd_s = wind_s[2]

	phi = angles[0]
	theta = angles[1]
	psi = angles[2]


	wn_b = cos(theta)*cos(psi)*wn_s + cos(theta)*sin(psi)*we_s + (-sin(theta))*wd_s
	we_b = (sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))*wn_s + (sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi))*we_s + sin(phi)*cos(theta)*wd_s
	wd_b = (cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi))*wn_s + (cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi))*we_s + cos(phi)*cos(theta)*wd_s
	

	wind_b = [wn_b, we_b, wd_b]
	return wind_b




