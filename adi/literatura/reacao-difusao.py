'''
Backward method to solve 1D reaction-diffusion equation:
    u_t = D * u_xx + alpha * u
    
with Neumann boundary conditions 
at x=0: u_x = sin(pi/2)
at x=L: u_x = sin(3*pi/4) with L=1
and initial condition u(x,0) = 4*x - 4*x**2
'''

import numpy as np
from scipy import sparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


M = 10 # GRID POINTS on space interval
N = 100 # GRID POINTS on time interval

x0 = 0
xL = 1

# ----- Spatial discretization step -----
dx = 0.04

t0 = 0
tF = 100.0

# ----- Time step -----
dt = 0.05

D = 1.0  # Diffusion coefficient
alpha = 1 # Reaction rate

r = dt*D/dx**2
s = dt*alpha
a = 1 + 2*r - s


xspan = np.linspace(x0, xL, M)
tspan = np.linspace(t0, tF, N)

main_diag = (1 + 2*r - s)*np.ones((1,M))
off_diag = -r*np.ones((1, M-1))

a = main_diag.shape[1]

diagonals = [main_diag, off_diag, off_diag]

A = sparse.diags(diagonals, [0,-1,1], shape=(a,a)).toarray()
A[0,1] = -2*r
A[M-1,M-2] = -2*r

# ----- Initializes matrix U -----
U = np.zeros((M, N))

#----- Initial condition -----
def CondicaoInicial(x):
    if(x<=10.0):
        U[:,0] =-1.0
    if(x>10.0):
        U[:,0]=1.0;
for i in range (M+1):
    xi=i*dx
    CondicaoInicial(xi)

for k in range(1, N):
    c = np.zeros((M-2,1)).ravel()
    b = b1 + b2  # Right hand side
    U[0:M, k]=(dt*U[0:M,k-1]*(ONE[0:M,k-1]-U[0:M,k-1]**2))+ U[0:M,k-1]

# ----- Checks if the solution is correct:
gc = np.allclose(np.dot(A,U[0:M,N-1]), b)
print(gc)

# ----- Surface plot -----
X, T = np.meshgrid(tspan, xspan)

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, T, U, linewidth=0,
                       cmap=cm.coolwarm, antialiased=False)

#ax.set_xticks([0, 0.05, 0.1, 0.15, 0.2])

#ax.set_xlabel('Time')
ax.set_ylabel('Space')
ax.set_zlabel('U')
plt.tight_layout()
plt.show()