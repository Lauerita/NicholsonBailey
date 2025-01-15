import numpy as np
import matplotlib.pyplot as plt

n = 2000
m = 2000
r = np.linspace(0, 12, 400)

def NB(H, P, r, a = 1, c = 1):
    H_t = r * H * np.exp(-a*np.sqrt(P))
    P_t = c*H*(1 - np.exp(-a*np.sqrt(P)))
    return np.array([H_t, P_t])

R = []
H = []

for i in r:
    # Changed to generate random numbers between 1 and 10
    x_val = 1 + 9 * np.random.random()
    y_val = 1 + 9 * np.random.random()
        
    for _ in range(m):
        x_val, y_val = NB(x_val, y_val, i)
        
    for _ in range(m):
        x_val, y_val = NB(x_val, y_val, i)
        
        H.append(y_val)
        R.append(i)

plt.figure(figsize = [10, 8])
plt.scatter(R, H, color = 'b', marker = '.', s = 0.1)
plt.xlabel('Parameter values')
plt.ylabel('H value at the fixed point')
plt.title('Bifurcation Diagram for r')
plt.grid(True)
plt.axvline(x = 4.92155, linestyle = '--')
#plt.savefig('Bif1.png', dpi = 300)

plt.figure(figsize = [10, 8])
plt.scatter(R, H, color = 'b', marker = '.', s = 1)
plt.xlabel('Parameter values')
plt.ylabel('H value at the fixed point')
plt.title('Bifurcation Diagram for r')
plt.xlim([4.9, 5.3])
plt.ylim([0,10])
plt.grid(True)
plt.axvline(x = 4.92155, linestyle = '--')
#plt.savefig('Bif2.png', dpi = 300)



