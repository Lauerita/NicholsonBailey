import numpy as np
import matplotlib.pyplot as plt



def NB(H, P, r, a = 1, c = 1):
    H_t = r * H * np.exp(-a*np.sqrt(P))
    P_t = c*H*(1 - np.exp(-a*np.sqrt(P)))
    
    return np.array([H_t, P_t])



fig,ax = plt.subplots() 
    
def NB_run(H0, P0):
    
    steps = 500
    storage = np.zeros([steps, 2])
    # r = 5.18
    r = 11.3
    #r = 6.06
    #r = 17.2
    storage[0,0] = H0
    storage[0,1] = P0

    for i in range(1, steps):
        
        storage[i, 0], storage[i, 1] = NB(storage[i-1, 0], storage[i-1, 1], r)
        
    return storage
        
    

H0 = np.linspace(0, 10, 200)
P0 = H0

initials_H = []
initials_P = []
orbit_H = [] # fixed points 
orbit_P = [] # fixed points 



for i in range(len(H0)):
    
    for j in range(len(P0)):
        
        y = NB_run(H0[i], P0[j])
        #print(abs(y[-2,0] - y[-1,0]))

        #print(y[-2,:] - y[-1,:])

        
        if abs(y[-2,0] - y[-1,0]) < 1e-6 and abs(y[-2,1] - y[-1,1]) < 1e-6 : 
            initials_P.append(P0[j])
            initials_H.append(H0[i])
            orbit_H.append(y[-1,0])
            orbit_P.append(y[-1,1])
            
            
plt.xlabel('initial values of  H')
plt.ylabel('initial values of P')
plt.scatter(initials_H, initials_P, color = 'black', s = 1)
plt.scatter(orbit_H, orbit_P, color = 'r', s = 1)
#plt.title('Basin of Attraction for r = 5.18')
#plt.savefig('BasinofAttraction2.png', dpi = 300)

