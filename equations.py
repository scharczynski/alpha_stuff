import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import math

d = np.random.normal(loc=300, size=1000)


def big_F(s):
    f = np.zeros(1000)
    f[100] = 1

    # f[100:130] = 1
    # f[200:250] = 1
    # f[300:330] = 1
    # f[400:450] = 1
    F_arr = []
    F = 0
    a = 1
    for ind, t in enumerate(d):
        F_arr.append(F)
        F += a*(-s*F+f[ind])
    return np.array(F_arr)
        
def deriv_time(fun, s,k, time):
    deriv_matrix = np.zeros((11, 11))
    spacing = 0.1
    # for i in range(2, 10):
        # deriv_matrix[i,i-1] = -(big_F(i+1)[time]-big_F(i)[time])/(big_F(i)[time]-big_F(i-1)[time])/(big_F(i+1)[time] - big_F(i-1)[time])
    # deriv_matrix = ((fun(s+1)[time]-big_F(s)[time])/(big_F(i)[time]- big_F(i-spacing)[time])-(big_F(i)[time]-big_F(i-spacing)[time])/(big_F(i+spacing)[time]-big_F(i)[time]))/(big_F(i+spacing)[time] - big_F(i-spacing)[time])
        # deriv_matrix[i,i+1] = (big_F(i)[time]-big_F(i-1)[time])/(big_F(i+1)[time]-big_F(i)[time])/(big_F(i+1)[time] - big_F(i-1)[time])
    s_list = np.arange(s-k, s+k+1)/10
    f_of_s = list(map(lambda x: fun(x)[time], s_list))
    deriv = np.gradient(f_of_s)
    for ds in range(k-1):
        deriv = np.gradient(deriv)
    # f, f_p, f_m = fun(s), fun(s+1), fun(s-1)
    # for ds in range(k):
    #     deriv = (f_p - f_m) / 2
        
    return deriv[s-k]

k = 4
s = 10
tau = -k/s
f_tilde = []
f_capital = []

for time in range(0, 1000):
    F = big_F(s)
    # F_diff = np.dot(np.linalg.matrix_power(deriv_time(time, s), k), F[time])
    # f = np.array(s)**(k+1)*(deriv_time(time)**k)
    result = deriv_time(big_F, s, k, time)

    # F_diff = np.array(s)**(k+1)*(deriv_time(time, s))**k#*F[time])
    F_diff = np.array(s)**(k+1)*result


    L1 = (-1)**k*s**(k+1) # this can be taken out
    L2 = (F_diff/math.factorial(k))# dont knw what this is for *(self.Taustarlist**self.g)
    if L1*L2 < 0:
        f_tilde.append(0)
    else:
        f_tilde.append(L1*L2)
    f_capital.append(F[time])
    # plot_arr.append(til_f[:, time])
plt.plot(f_tilde)
# plt.plot(f_capital)
plt.show()
print("stop1")