import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import math
import autograd
from scipy.optimize import basinhopping
import json


def big_F(s_arr, f):

    # F = np.zeros_like(data)
    # F[0] = 0
    a=1
    F = 0
    F_arr = np.zeros((len(s_arr), len(f)))
    # F_arr[:,0] = 1
    # f = args[1]
    for ind, t in enumerate(f[1:]):
        F_arr[:,ind+1] = F_arr[:,ind] + a*(-s_arr*F_arr[:,ind]+f[ind])



    return np.array(F_arr)

def deriv_time(f_of_s, s,k, time):
    s_list = np.arange(s-k, s+k+1, 0.1)
    deriv = np.diff(f_of_s)
    for ds in range(k-1):
        deriv = np.diff(deriv)
    index = np.where(np.isclose(s_list,s))[0]
    return deriv[index]

def invert(x, *args):
    s = x
    k = args[2]
    f_of_s = args[3]
    f_tilde = []
    # s_list = np.arange(s-k, s+k+1, 0.1)
    # f_of_s = np.array(list(map(lambda x: big_F(x, *args), s_list)))
    s_list = args[4]
    s_index = np.where(np.isclose(s_list,s))[0]
    for time in range(0, 53092):
    # for time in range(0, 500):
        result = deriv_time(f_of_s[:,time], s, k, time)
        F_diff = np.array(s)**(k+1)*result
        L1 = (-1)**k*s**(k+1) # this can be taken out
        L2 = (F_diff/math.factorial(k))# dont knw what this is for *(self.Taustarlist**self.g)
        if L1*L2 <= 0:
            f_tilde.append(0)
        else:
            f_tilde.append(L1*L2)
    
    return f_tilde

def load_data(path):
    with open(path, 'r') as f:
        return np.array(json.load(f))

def likelihood(x, *args):
    # if not np.isnan(x):
    f_tilde = np.array(invert(x, *args), dtype=float)
    spikes = args[0]
    func = np.ma.masked_invalid(f_tilde)
    obj = np.sum(spikes * (-np.log(func)) +
        (1 - spikes) * (-np.log(1 - (func))))
    # else:
    #     return np.inf

    return obj

def gradient_respecting_bounds(bounds, fun, eps=1e-8):
    """bounds: list of tuples (lower, upper)"""
    def gradient(x):
        fx = fun(x)
        grad = np.zeros(len(x))
        for k in range(len(x)):
            d = np.zeros(len(x))
            d[k] = eps if x[k] + eps <= bounds[k][1] else -eps
            grad[k] = (fun(x + d) - fx) / d[k]
        return grad
    return gradient




def minimize(lb, ub):
    spikes = load_data("/Users/stevecharczynski/Desktop/cromer_spikes_1D.json")[:53092]
    starts = load_data("/Users/stevecharczynski/Desktop/cromer_trial_start_1D.json")[:53092]
    # spikes = np.zeros(500)
    # spikes[105:110] = 1
    # spikes[112:116] = 1
    # spikes[125:127] = 1
    # spikes[139] = 1
    # spikes[160] = 1
    # spikes[350:355] = 1
    # spikes[360:362] = 1
    # spikes[366:340] = 1
    # spikes[344] = 1


    # f = np.zeros(500)
    # f[100] = 1
    # f[300] = 1
    # f[300:330] = 1
    # starts = f
    k = 4
    s_list = np.arange(lb[0]-k, ub[0]+k+1, 0.1)
    # f_of_s1 = np.array(list(map(lambda x: big_F(x, starts), s_list)))
    f_of_s = np.array(big_F(s_list, starts))
    print("yesy")
    out = scipy.optimize.minimize(
        likelihood,
        x0 = [0.002],
        method="TNC",
        jac=None,
        bounds = [(lb[0], ub[0])],
        args = (spikes, starts, k, f_of_s, s_list)
        
    )
    plt.plot(invert(out.x, *((spikes, starts, k, f_of_s, s_list))))
    # plt.plot((spikes))
    plt.show()

    return out

    # stepper = RandomDisplacementBounds(lb, ub, stepsize=step)
    # accepter = AccepterBounds(ub, lb)
    # out = basinhopping(
    #         likelihood,
    #         x0=[0.1, 4],
    #         niter=niter,
    #         disp=True,
    #         accept_test=accepter,
    #         take_step=stepper,  
    #         stepsize=step,
    #         minimizer_kwargs={"method":"TNC", "bounds":[[0.01, 20], [1, 5]], "jac":False, "args":(spikes, starts)},
    #         interval=10,
    #     )
        
    # fit = out.x
    # fun = out.fun
    # return fit, fun

x = minimize([0.001], [5.0])
print("hi")







#   for time in range(0, 500):
#         result = deriv_time(f_of_s[:,time], s, k, time)
#         F_diff = np.array(s)**(k+1)*result
#         L1 = (-1)**k*s**(k+1) # this can be taken out
#         L2 = (F_diff/math.factorial(k))# dont knw what this is for *(self.Taustarlist**self.g)
#         if L1*L2 <= 0:
#             f_tilde.append(0)
#         else:
#             f_tilde.append(L1*L2)