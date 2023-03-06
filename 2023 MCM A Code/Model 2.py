import numpy as np
import matplotlib.pyplot as plt   
from gekko import GEKKO

def ode_RK4(f, X_0, dt, T): 
    '''RK4 method to solve ODE'''   
    N_t = int(round(T/dt))
    # Initial conditions
    usol = [X_0]
    u = np.copy(X_0)

    tt = np.linspace(0, N_t*dt, N_t + 1)
    # RK4
    for t in tt[:-1]:
        u1 = f(u + 0.5*dt* f(u, t), t + 0.5*dt)
        u2 = f(u + 0.5*dt*u1, t + 0.5*dt)
        u3 = f(u + dt*u2, t + dt)
        u = u + (1/6)*dt*( f(u, t) + 2*u1 + 2*u2 + u3)
        usol.append(u)
    return usol, tt


#read initial data
a=np.genfromtxt('a.csv', delimiter=',')
N=np.genfromtxt('N.csv', delimiter=',')
L=np.genfromtxt('L.csv', delimiter=',')
K=np.genfromtxt('K.csv', delimiter=',')
x0=np.genfromtxt('x0.csv', delimiter=',')

di=np.zeros(20)
def f(u,t):
    x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20=u
    xi=np.array([x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20])
    for i in range(20):
        h=49.1*np.exp(-0.00104375*(t-229.5)**2)
        r_t= L[i]/(1+np.exp(-K[i]*h))
        di[i]= r_t*xi[i]*(1-xi[i]/N[i]-sum(a[i,:]*xi[:]/N[:]))
    return di

ode_RK4(f, x0, 0.1, 365)

def fig1():
    global u1,u2,u3,u4,u5,u6,u
    u, t = ode_RK4(f, x0, 0.01, 365)
    u1 = [a[0] for a in u]
    u2 = [b[1] for b in u]
    u3 = [c[2] for c in u]
    u4 = [d[3] for d in u]
    u5 = [e[4] for e in u]
    u6 = [f[5] for f in u]
    u7 = [g[6] for g in u]
    u8 = [h[7] for h in u]
    u9 = [i[8] for i in u]
    u10 = [j[9] for j in u]
    u11 = [K[10] for k in u]
    u12 = [l[11] for l in u]
    u13 = [m[12] for m in u]
    u14 = [n[13] for n in u]
    u15 = [o[14] for o in u]
    u16 = [p[15] for p in u]
    u17 = [q[16] for q in u]
    u18 = [r[17] for r in u]
    u19 = [s[18] for s in u]
    u20 = [t[19] for t in u]
    plt.plot(t, u1)
    plt.plot(t, u2)
    plt.plot(t, u3)
    plt.plot(t, u4)
    plt.plot(t, u5)
    plt.plot(t, u6)
    plt.plot(t, u7)
    plt.plot(t, u8)
    plt.plot(t, u9)
    plt.plot(t, u10)
    plt.plot(t, u11)
    plt.plot(t, u12)
    plt.plot(t, u13)
    plt.plot(t, u14)
    plt.plot(t, u15)
    plt.plot(t, u16)
    plt.plot(t, u17)
    plt.plot(t, u18)
    plt.plot(t, u19)
    plt.plot(t, u20)
    plt.show()


av = []
for i in range(20):
    av.append(sum(u[i])/len(u[i]))

m = GEKKO()
x = m.Array(m.Var,[20,1],integer=True,lb=0,ub=1) # set decision variables
for i in range(0,20):
        x[i,0].value=1 #initial guess


m.Equation(sum(x[:,0]*av[:])<=0.9*sum(N[:])/3) 
for i in range(0,20):
    m.Equation(x[i,0]*av[i] <= N[i]) 

a1=52 #bio mass for class 1
a2=0.224 #bio mass for class 2
a3=0.028 #bio mass for class 3
a4=0.016 #bio mass for class 4


m.options.SOLVER = 1 # sole the system
m.Minimize(-(sum(a1*x[0:3,0]*av[0:3])+sum(a2*x[3:8,0]*av[3:8])+sum(a3*x[6:14,0]*av[6:14])+sum(a4*x[14:20,0]*av[14:20])))
m.solve(disp=False)
print(x)

np.savetxt('av', av, delimiter=",")
