import numpy as np
from scipy.special import comb
#import matplotlib.pyplot as plt

def smoothstep(x, t_stamp, trans_dur,y=[0,1], N=1):
    x = np.clip((x - t_stamp+trans_dur/2) / trans_dur,0,1)
    
    result = 0
    for n in range(0, N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n

    result *= x ** (N + 1)
    return result*(y[1]-y[0])+y[0]
 
def generateSmoStepData(xs, t_stamp, trans_dur,y=[0,1], N=0,label=0):
    print(label,end='\r')
    return [smoothstep(x,t_stamp,trans_dur,y,N) for x in xs]


# xs = np.linspace(-0.5, 0.5, 1000)

# test=generateSmoStepData(xs, 0, 0.5,y=[-1,2], N=1,label=0)

# plt.plot(xs, test)

# plt.legend()
# plt.show()