import numpy as np
import matplotlib.pyplot as plt

def NB(X, r, a, c):
    H_t = r * X[0] * np.exp(-a * np.sqrt(X[1]))
    P_t = c * X[0] * (1 - np.exp(-a * np.sqrt(X[1])))
    return np.array([H_t, P_t])

# Create figure with two subplots in one column
plt.figure(figsize=(10, 8))

# Bifurcation plot
plt.subplot(2, 1, 1)
r_range = np.linspace(10, 12, 500)
n = 2000
m = 2000

R = []
H = []
for i in r_range:
    x_val = 1 + 9 * np.random.random()
    y_val = 1 + 9 * np.random.random()
    
    for _ in range(m):
        x_val, y_val = NB(np.array([x_val, y_val]), i, 1, 1)
    
    for _ in range(m):
        x_val, y_val = NB(np.array([x_val, y_val]), i, 1, 1)
        
        H.append(x_val)
        R.append(i)

plt.scatter(R, H, color='b', marker='.', s=0.1)
plt.title('Bifurcation Diagram')
plt.xlabel('r')
plt.ylabel('H value')
plt.ylim([0,10])

# Lyapunov Exponent plot
plt.subplot(2, 1, 2)
r_sw = np.linspace(10, 12, 500)
LE1 = []
LE2 = []

for i in r_sw:
    # Parameters
    Ntrans = 1000  # transients
    Nit = 1000     # number of iterations
    
    # Map parameters
    r = i
    a = 1
    c = 1
    
    # Separation for contraction
    epsilon = 1e-12
    
    # Initial conditions
    x0 = 0.1
    y0 = 0.1
    H = np.array([x0, y0])
    
    # Initialize variables
    sumLog = 0
    
    # Transients
    for _ in range(Ntrans):
        H = NB(H, r, a, c)
    
    # Initial basis
    p0 = H
    p1 = H + epsilon * np.array([1, 0])
    p2 = H + epsilon * np.array([0, 1])
    
    # Main loop for Lyapunov exponent computation
    for _ in range(Nit):
        # Apply map to all points
        H_p0 = NB(p0, r, a, c)
        H_p1 = NB(p1, r, a, c)
        H_p2 = NB(p2, r, a, c)
        
        # Create matrix of separated points
        M = np.column_stack((H_p1 - H_p0, H_p2 - H_p0))
        
        # QR decomposition
        Q, R = np.linalg.qr(M)
        
        # Accumulate log of scaled diagonal elements
        sumLog += np.log(np.abs(np.diag(R)) / epsilon)
        
        # Update points for next iteration
        p0 = H_p0
        p1 = H_p0 + epsilon * Q[:, 0]
        p2 = H_p0 + epsilon * Q[:, 1]
    
    # Compute Lyapunov Exponent
    E, F = sumLog / Nit
    LE1.append(E)
    LE2.append(F)

plt.plot(r_sw, LE1, color='b')
#plt.plot(r_sw, LE2, color='orange', label='Lyapunov Exponent 2')
plt.title('Lyapunov Exponent')
plt.xlabel('r')
plt.ylabel('Lyapunov Exponent')
#plt.xlim([4.9, 5.3])
plt.axvline(x = 10.24, linestyle='--', color='r')
plt.legend()

plt.tight_layout()
plt.savefig('CombinedBifandLE.png', dpi = 300)
plt.show()
