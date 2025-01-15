
import numpy as np 
import matplotlib.pyplot as plt

r = 10.24
a = 1
c = 1


H_list = []
P_list = []

def NB(H, P, r, a, c):
    H_t = r * H *np.exp(-a*np.sqrt(P))
    P_t = c*H*(1 - np.exp(-a*np.sqrt(P)))
    return np.array([H_t, P_t])

H0 = 6
P0 = 1
H, P = NB(H0, P0, r, a, c)

t = np.arange(1, 500, 1)

for i in range(len(t)):
    H, P = NB(H, P, r, a, c)
    H_list.append(H)
    P_list.append(P)

fixed_x = r/(r-1) * (np.log(r)/a)**2
fixed_y = (np.log(r)/a)**2



plt.figure()
plt.plot(t, H_list, label = 'Host population')
plt.plot(t, P_list, label = 'Parastoid population')
plt.title('General Solutions')
plt.xlabel('time')
plt.ylabel('population size')
plt.grid(True)
