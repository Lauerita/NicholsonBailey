import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# System parameters
r = 5  # growth rate
a = 1.0  # searching efficiency
c = 1.0  # conversion efficiency
Nt = 300          # Number of iterates
Nx = 700          # Grid resolution
Lx = 10.0          # Plot range for H
Ly = 10.0          # Plot range for P
distp = 1e-2      # Threshold for convergence
distM = 1e+5      # Cap before considering divergence
plottransient = 20  # Number of transient plots to show

# Create custom colormap (gray for undecided, red for divergence, black for convergence)
mycolmap = ListedColormap([[0.8, 0.8, 0.8], [1, 0, 0], [0, 0, 0]])

# Find fixed points numerically
# For simple fixed point, solve H* = r*H*exp(-a*sqrt(P*)) and P* = c*H*(1-exp(-a*sqrt(P*)))
def find_fixed_points():
    # Create a grid of initial guesses
    H_guess = np.linspace(0.1, 3.0, 20)
    P_guess = np.linspace(0.1, 3.0, 20)
    fixed_points = []
    
    for h in H_guess:
        for p in P_guess:
            # Use fixed point iteration to find equilibria
            H, P = h, p
            for _ in range(100):
                H_new, P_new = NB(H, P, r, a, c)
                if abs(H_new - H) < 1e-6 and abs(P_new - P) < 1e-6:
                    # Check if this fixed point is already found
                    is_new = True
                    for fp in fixed_points:
                        if abs(fp[0] - H_new) < 1e-4 and abs(fp[1] - P_new) < 1e-4:
                            is_new = False
                            break
                    if is_new and H_new > 0 and P_new > 0:
                        fixed_points.append((H_new, P_new))
                    break
                H, P = H_new, P_new
                
    return fixed_points

# Nicholson-Bailey map function
def NB(H, P, r, a=1, c=1):
    H_t = r * H * np.exp(-a*np.sqrt(P))
    P_t = c*H*(1 - np.exp(-a*np.sqrt(P)))
    return H_t, P_t

# Create mesh grid for initial conditions
x = np.linspace(0, Lx, Nx)  # Start from 0 since H and P are populations
y = np.linspace(0, Ly, Nx)
H0, P0 = np.meshgrid(x, y)
H, P = H0.copy(), P0.copy()

# Find fixed points
fixed_points = find_fixed_points()
print(f"Found {len(fixed_points)} fixed points:")
for i, fp in enumerate(fixed_points):
    print(f"Fixed point {i+1}: H* = {fp[0]:.4f}, P* = {fp[1]:.4f}")

# Create figure for transient plots
if plottransient:
    fig_trans = plt.figure(figsize=(10, 8))

# Iterate the map
print('Iterating... ', end='')
for ii in range(1, Nt+1):
    print(f'{ii},', end='')
    H_new, P_new = NB(H, P, r, a, c)
    
    # Handle potential numerical issues
    mask = np.isnan(H_new) | np.isnan(P_new) | (H_new > distM) | (P_new > distM)
    H_new[mask] = distM
    P_new[mask] = distM
    
    H, P = H_new, P_new
    
    if ii <= plottransient:
        plt.figure(fig_trans.number)
        plt.clf()
        
        # Calculate minimum distance to any fixed point
        if fixed_points:
            D = np.inf * np.ones_like(H)
            for fp in fixed_points:
                D = np.minimum(D, (H - fp[0])**2 + (P - fp[1])**2)
            D = np.sqrt(D)
            
            plt.pcolormesh(H0, P0, np.log10(D), shading='auto')
            plt.colorbar(label='Log10 distance to nearest fixed point')
            
            # Plot fixed points
            for fp in fixed_points:
                plt.plot(fp[0], fp[1], 'r*', markersize=10)
        
        plt.axis('equal')
        plt.xlabel('Host Population (H)')
        plt.ylabel('Parasitoid Population (P)')
        plt.title(f'Distance to fixed points after {ii} iterates')
        plt.pause(0.01)

print('\n')

# Create final plot
fig_final = plt.figure(figsize=(10, 8))

# Calculate final distances to fixed points
if fixed_points:
    D = np.inf * np.ones_like(H)
    for fp in fixed_points:
        D = np.minimum(D, (H - fp[0])**2 + (P - fp[1])**2)
    
    B = np.zeros_like(D)
    B[D > distM] = 1     # Diverged
    B[D < distp] = -1    # Converged to fixed point
    
    plt.pcolormesh(H0, P0, B, cmap=mycolmap, shading='auto')
    
    # Plot fixed points
    for fp in fixed_points:
        plt.plot(fp[0], fp[1], 'b*', markersize=10, label='Fixed Points')

plt.axis('equal')
plt.xlabel('Host Population (H)')
plt.ylabel('Parasitoid Population (P)')
plt.title('Basins of Attraction for Nicholson-Bailey Model')
plt.colorbar(label='Basin Classification')
plt.legend()
plt.show()
