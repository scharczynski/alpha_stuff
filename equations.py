import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import math

d = np.random.normal(loc=300, size=1000)


def big_F(s):
    f = np.zeros(500)
    f[100] = 1

    # f[100:130] = 1
    # f[200:250] = 1
    # f[300:330] = 1
    # f[400:450] = 1
    F_arr = []
    F = 0
    a = 1
    for ind, t in enumerate(f):
        F_arr.append(F)
        F += a*(-s*F+f[ind])
    return np.array(F_arr)

def update(time_index, f=0, alpha=1):
    f = np.zeros(500)
    f[100] = 1
    F[:,time_index] = F[:,time_index-1]+(alpha*(-s_list.T*F[:,time_index-1]+f[time_index-1])*0.01)
    # to avoid numerical errors if small t grows too high
    if np.max(F[:,time_index]) > 0.01:
        F[:,time_index] = np.zeros(N)
    invert(time_index)

def invert(time):
    F_diff = np.dot(np.linalg.matrix_power(_DerivMatrix, k), F[:,time])

    L1 = (-1)**k*s_list**(k+1) # this can be taken out
    L2 = (F_diff/math.factorial(k))
    f_tilde[:,time] = L1.T*L2.T
    f_tilde[f_tilde[:,time]<0,time] = 0

        
def deriv_time(f_of_s, s,k, time):
    # spacing = 0.1
    # num_s = int((np.abs(s-k)+s+k+1)/spacing)
    deriv = np.diff(f_of_s)
    for ds in range(k-1):
        deriv = np.diff(deriv)
    index = np.where(np.isclose(s_list,s))[0]
    return deriv[index]

k = 4
tau = -400

f_capital = []
tau_of_t = []
taus = []
tstr_max = 2
tstr_min = 0.1
buff_len=10
N = buff_len+2*k

a = (tstr_max/tstr_min)**(1./buff_len)-1
pow_vec = np.arange(-k,buff_len + k) #-1
tau_list = (tstr_min * (1 + a)**pow_vec)
s_list = k/tau_list


# tau_list = np.geomspace(0.1, 200, num=N)
# s_list = (k/tau_list)

# for t in tau_list[k:len(tau_list)-k]:
# for t in np.arange(500, 540, 5):
# for t in Taustarlist[k:len(Taustarlist)-k]:

# s = k/t

_DerivMatrix = np.zeros((N,N))
for i in range(1, N-1):
    _DerivMatrix[i, i-1] = -(s_list[i+1]-s_list[i])/(s_list[i]-s_list[i-1])/(s_list[i+1] - s_list[i-1])
    _DerivMatrix[i, i] = ((s_list[i+1]-s_list[i])/(s_list[i]- s_list[i-1])-(s_list[i]-s_list[i-1])/(s_list[i+1]-s_list[i]))/(s_list[i+1] - s_list[i-1])
    _DerivMatrix[i, i+1] = (s_list[i]-s_list[i-1])/(s_list[i+1]-s_list[i])/(s_list[i+1] - s_list[i-1])
# s_list = np.arange(s-k, s+k+1, 0.1)
# s_list = np.geomspace(0.01, 10, len(np.arange(-k,buff_len + k)))

# s_list_log = np.log(s_list)
# s_list = np.logspace(s-k, s+k, num=25)
# f_of_s = np.array(list(map(lambda x: big_F(x), s_list)))
f_tilde = np.zeros((N, 500))
F = np.zeros((N,500))

for time in range(0, 500):
    update(time)
    # result = deriv_time(f_of_s[:,time], s, k, time)
    # F_diff = np.array(s)**(k+1)*result

for i, tau in enumerate(tau_list[1::2]):
    plt.plot(np.array(f_tilde[i, :]/sum(f_tilde[i,:])).T, label=tau)
plt.legend()
plt.show()

