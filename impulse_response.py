import numpy as np
# import autograd.numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import math
import autograd
from scipy.optimize import basinhopping
import json
import math


def func(a, s, k, f):
    f_nonzero = np.where(f!=0)[0] 
    f_tilde = np.zeros_like(f) + 10e-30
    t_s = -k/s
    for t, i in enumerate(f):
        for t_i in f_nonzero:
            if t_i < t:
                # f_tilde[t] += a*((1/t_s)*(k**(k+1)/math.factorial(k))* ((-t/t_s)**k) * ( np.exp(-s*(t-t_i))))
                # f_tilde[t] += a*((1/(t_i-t))*(k**(k+1)/math.factorial(k))* ((t_i-t/t_s)**(k+1)) * (np.exp(k*((t-t_i)/t_s))))
                p1 = (-1)**k / math.factorial(k)
                p2 = (-k/t_s)**(k+1)
                p3 = (t_i-t)**k
                p4 = np.exp(k*(1/t_s)*(t-t_i))
                f_tilde[t] += a*(p1 * p2 * p3 * p4)

    # t = np.arange(f.shape)

    return f_tilde
    



def load_data(path):
    with open(path, 'r') as f:
        return np.array(json.load(f))

def likelihood(x, *args):
    s, a = x
    k = 4
    spikes = args[0]
    f = args[1]
    fun = func(a,s,k,f)
    
    obj = np.sum(spikes * (-np.log(fun)) +
        (1 - spikes) * (-np.log(1 - (fun))))
    return obj

def minimize(spikes, f):
    spikes = np.zeros(500)
    # spikes[40:90] =  1
    # spikes[200:250] = 1
    spikes[112:116] = 1
    # spikes[410:490] = 1
    spikes[139] = 1
    # # spikes[160] = 1
    # spikes[350:355] = 1
    # spikes[360:362] = 1
    # spikes[320:360] = 1
    spikes[400:440] = 1
    spikes[470:495] = 1
    # spikes[520:570] = 1
    # spikes[820:930] = 1
    # spikes[1800:1875] = 1
    # spikes[2800:2900] = 1
    # spikes[3400:3480] = 1
    # spikes[3900:4976] = 1
    # spikes[344] = 1

    # spikes = spikes[300:]

    f = np.zeros(500)
    f[1] = 1
    f[300] = 1
    
    # f = f[300:]

    # f[500] = 1
    # # f[300] = 1
    # # f[300:330] = 1
    x0 = [0.1, 1e-5]

    
            
    # minimizer_kwargs = {"method":solver_params["method"], "bounds":self.bounds, "jac":autograd.jacobian(self.objective)}
    minimizer_kwargs = {"method":"TNC", "bounds":[[0.001, 0.5], [0, 1]], "jac":False, "args":(spikes, f)}
        # second_pass_res = basinhopping(
        #     self.objective,
        #     self.x0,
        #     disp=solver_params["disp"],
        #     T=solver_params["T"],
        #     niter=solver_params["niter"],
        #     accept_test=accepter,
        #     take_step=stepper,  
        #     stepsize=solver_params["stepsize"],
        #     minimizer_kwargs=minimizer_kwargs,
        #     interval=solver_params["interval"],
        # )
    res = scipy.optimize.basinhopping(
        likelihood, 
        x0,
        niter=1,
        minimizer_kwargs=minimizer_kwargs
        
    )
    test_arr = []

    output = func(1e-10, res.x[0], 4, f)
    print(res.x, res.fun)
    # plt.xlim(0, 500)
    plt.plot(output, label=4/res.x[0])
    
    plt.plot(f*max(output))
    plt.plot(spikes*max(output)/2)
    # for s in [0.001, 0.002, 4, 0.8]:
    #     p = func(1e-10, s, 4, f)
    #     plt.plot(p, label=s)
    plt.legend()
    plt.show()


spikes = load_data("/Users/stevecharczynski/workspace/data/cromer/spikes_1d/2_spikes_1D.json")[0:4000]
starts = load_data("/Users/stevecharczynski/workspace/data/cromer/spikes_1d/2_starts_1D.json")[0:4000]
# x = minimize(spikes, starts)
# spikes = load_data("/Users/stevecharczynski/workspace/data/cromer/spikes_1d/3_spikes_1D.json")
# starts = load_data("/Users/stevecharczynski/workspace/data/cromer/spikes_1d/3_starts_1D.json")
# y = minimize(spikes, starts)

minimize(spikes, starts)