import numpy as np
import sympy as sp
from scipy import sparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.integrate import quad

M = 100 # GRID POINTS on space interval
N = 100 # GRID POINTS on time interval

x0 = 0.0
xL =1.0

# ----- Spatial discretization step -----
dx = 0.04

t0 = 0
tF = 1.0

# ----- Time step -----
dt = 0.5

D = 0.0001  # Diffusion coefficient
alpha = 1 # Reaction rate
tau=1.0
#w0=0.00
#w1=0.00
#f=6.0*(w0-w1)
r = dt*D/(dx*dx)
s = dt*alpha
a = 1 + 2*r

alpha2=2.0

xspan = np.linspace(x0, xL, M)
tspan = np.linspace(t0, tF, N)

main_diag = (1 + 2*r)*np.ones((1,M))
off_diag = -r*np.ones((1, M-1))

a = main_diag.shape[1]

diagonals = [main_diag, off_diag, off_diag]

A = sparse.diags(diagonals, [0,-1,1], shape=(a,a)).toarray()
A[0,1] = -2*r
A[M-1,M-2] = -2*r

# ----- Initializes matrix U -----
U = np.zeros((M, N))
ONE= np.ones((M, N))
#----- Initial condition -----
for i in range(0,M):
    if(xspan[i]>=0.4 and xspan[i]<=0.6):
        U[i,0] = 1.0
    else:
        U[i,0]= 0.0
def h(x):
    return (x*x)*(3-(2*x))
v0=0.6
I0 = quad(h, xspan[0], xspan[0])
vu=I0[0]
print('vu',vu)
f0=alpha2*(v0-vu)
print('f0',f0)
for k in range(1, N):
    
    #print('vu',vu)
    #print('U',U[0:M,k-1])
    b1 = np.asarray(((dt*U[0:M,k-1]*(ONE[0:M,k-1]-U[0:M,k-1]))*(U[0:M,k-1]-0.5+f0))/tau , dtype='float')
    b2 = np.array (U [0: M, k-1])
    #print('b2',b2)
    b = b1 + b2 # lado direito
    U [0: M, k] = np.linalg.solve (A, b) # Solve x = A \ b
    I0 = quad(h, xspan[0], xspan[k])
    vu=I0[0]
    print('vu',vu)
    f0=alpha2*(v0-vu)
    print('f0',f0)
 
#for i in range(0,M):
plt.plot(xspan,U)
#plt.plot(xspan,U[0:M,50])
#plt.plot(xspan, U[0:M,99])
#plt.legend(['U[0]','U[50]','U[99]'])
plt.title('v0=0.6')
plt.show()

# ----- Surface plot -----
X, T = np.meshgrid(tspan, xspan)

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, T, U, linewidth=0,
                       cmap=cm.coolwarm, antialiased=False)

ax.set_xlabel('Time')
ax.set_ylabel('Space')
ax.set_zlabel('U')
plt.tight_layout()
plt.show()
