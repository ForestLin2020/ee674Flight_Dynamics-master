
import sys
sys.path.append('..')
import numpy as np
import math
from chap4.wind_simulation import wind_simulation
from tools.transfer_function import transfer_function

state = np.array([[1],  # (0)
                 [2],  # (1)
                 [3],  # (2)
                 [4],  # (3)
                 [5]])  # (12)
print(state[0:5])
e0=state[0]
e1=state[1]
e2=state[2]
e3=state[3]
e4=state[4]
print(e2)
# y = np.zeros((7,1))
# print(y)
# print("y[1]=",y[1])
# print("y.item(1)=",y.item(1))


# print("-------")
# z = np.zeros((5,))
# print(z[2])