import numpy as np


# print(np.eye(2))
# # print(np.eye(3))
# # print(np.eye(3, k=1))
# #
# # r = np.array([23,25])
# # print(r)
# # print(r[0])
# #
# # s = [34,45]
# # print(s)
# # print(s[0])
a = np.eye(2)*1.*10**(-9.)
# b = np.array([5,6,7]).T
c = 1e-9*np.diag([1.0, 1.0])

# b = np.array([[3,4,5,6],
#               [7,8,9,10]])
# a = np.diag([1.0,1.0])
# print(a+b)
print('a = ', a)
print('b = ', c)