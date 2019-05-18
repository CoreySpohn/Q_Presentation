#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:47:24 2019

@author: corey
"""

###############################################################################
# NOTES
# To get Latex with dots: print(vlatex(L))
###############################################################################
from sympy import *
from sympy.physics import *
from sympy.physics import vector
from sympy.physics.vector import time_derivative as dt
from sympy.physics.vector import dot
from sympy.physics.vector.printing import vpprint, vlatex
import numpy as np
import control
import matplotlib.pyplot as plt

vector.init_vprinting()

###############################################################################
# DYNAMICS
###############################################################################
# Variable creation
l, l1, l2, m, m1, m2, M, g, u = symbols('l l_1 l_2 m m_1 m_2 M g u')
t1, t2, x = vector.dynamicsymbols('theta_1 theta_2 x')

# Create reference frame and unit vectors
N = vector.ReferenceFrame('N')

er1 = -sin(t1)*N.x+cos(t1)*N.y
et1 = -cos(t1)*N.x-sin(t1)*N.y

er2 = sin(t2)*N.x+cos(t2)*N.y
et2 = cos(t2)*N.x-sin(t2)*N.y

# Position vectors
r1 = x*N.x+l1*er1
r2 = x*N.x+l2*er2

# Velocity vectors
v1 = dt(r1,N)
v2 = dt(r2,N)

# Kinetic and potential energies
T = trigsimp((1/2)*M*dot(dt(x*N.x,N),dt(x*N.x,N)) + (1/2)*m1*dot(v1,v1) + (1/2)*m2*dot(v2,v2))
U = g*(m1*l1*cos(t1)+m2*l2*cos(t2))

# Lagrangian
L = T-U
F = u*N.x
rf = x*N.x

# Equations of motion
eqx = dt(L.diff(dt(x,N)),N)-L.diff(x)-dot(rf.diff(x,N),F)
eqt1 = dt(L.diff(dt(t1,N)),N)-L.diff(t1) - dot(rf.diff(t1,N),F)
eqt2 = dt(L.diff(dt(t2,N)),N)-L.diff(t2) - dot(rf.diff(t2,N),F)

# Solve for variables
sol = solve([eqx, eqt1, eqt2],[dt(x,N,2), dt(t1,N,2), dt(t2,N,2)], dict=True)
xdd = sol[0][dt(x,N,2)]
t1dd = sol[0][dt(t1,N,2)]
t2dd = sol[0][dt(t2,N,2)]

# Linearize
lsubs = [(sin(t1),t1), (sin(2.0*t1),2*t1), (cos(t1),1), (sin(t2),t2), 
         (sin(2.0*t2),2*t2), (cos(t2), 1), (t1*t1, 0), (t1*t2, 0), (t2*t2, 0),
         (dt(t1,N)*dt(t1,N), 0), (dt(t2,N)*dt(t2,N), 0), (1.0,1)]
xddL = expand_trig(xdd).subs(lsubs) # Lienarized solutions
t1ddL = expand_trig(t1dd).subs(lsubs)
t2ddL = expand_trig(t2dd).subs(lsubs)

###############################################################################
# CONTROLS
###############################################################################
state = [x, dt(x,N), t1, dt(t1,N), t2, dt(t2,N)]
stated = [dt(x,N), xddL, dt(t1,N), t1ddL, dt(t2,N), t2ddL]
outputs = [x, t1, t2]
inputs = [u]

n = len(state)
m = len(inputs)
p = len(outputs)

# A matrix
avals = []
for i in range(0,n):
    for j in range(0, n):
        avals.append(stated[i].diff(state[j]))
A = Matrix(n, n, avals)
    
# B matrixu
bvals = []
for i in range(0,n):
    for j in range(0, m):
        bvals.append(stated[i].diff(inputs[j]))
B = Matrix(n, m, bvals)

# C matrix
cvals = []
for i in range(0,p):
    for j in range(0, n):
        cvals.append(outputs[i].diff(state[j]))
C = Matrix(p, n, cvals)

# D matrix
D = np.zeros([int(C.shape[0]),1])

# Controllability matrix
# This is messier than I'd like but it works
Qi = []
for i in range(0,n):
    Qi.append((A**i * B))

Qcvals = []
for i in range(0,n):
    for j in range(0, n):
        Qcvals.append(Qi[j][i])

Qc = simplify(Matrix(n, n, Qcvals))

# Observability matrix
# This is messier than I'd like but it works
Qio = []
for i in range(0,n):
    Qio.append((C * A**i))

Qovals = []
for i in range(0,n):
    for j in range(0, n):
        Qovals.append(Qio[j][i])

Qo = simplify(Matrix(n, n, Qovals))

# Conditions for stability
l1stabcond = solve(Qc.det(),l1)[0]

# Numerical values for the matrices
l1n = 1
l2n = 1
m1n = 5
m2n = 1
Mn  = 5
gn  = 1
t10 = 1 # Initial angle
dt10 = 0
t20 = 1
dt20 = 0
x10 = 0
dx10 = 0
nSubs = [(l1,l1n), (l2,l2n), (m1,m1n), (m2,m2n), (M,Mn), (g,gn)]

An = np.array(A.subs(nSubs)).astype(float)
Bn = np.array(B.subs(nSubs)).astype(float)
Cn = np.array(C.subs(nSubs)).astype(float)
Qcn = np.array(Qc.subs(nSubs), dtype='float')
# Stability check
if np.linalg.matrix_rank(Qcn) == n:
    print('System is controllable\n')
else:
    print('System is uncontrollable\n')
vector.vprint(Qcn)

# Control system
sys = control.ss(An,Bn,Cn,D)
x0 = [x10, dx10, t10, dt10, t20, dt20]
tspan = np.linspace(0,5,1000)


# Step response
#t, yout = control.step_response(ss, tspan, x0)
#
#plt.rc('text', usetex=True)
#plt.plot(t,yout[0], label=r'$\theta_1$')
#plt.plot(t,yout[1], label=r'$\theta_2$')
#plt.xlabel(r'Time, t')
#plt.ylabel(r'Angle, $\theta$')
#plt.legend(loc='best')


# Sin input
#usin = 15*np.cos(tspan*2*np.pi)
#t, yout, xout = control.forced_response(sys, tspan, usin, x0)
#plt.rc('text', usetex=True)
#plt.plot(t, xout[0], label=r'$\theta_1$')
#plt.plot(t, xout[2], label=r'$\theta_2$')
#plt.xlabel(r'Time, t')
#plt.ylabel(r'Angle, $\theta$')
#plt.legend(loc='best')

# Calculate position x
#x = np.zeros(len(tspan))
#v = np.zeros(len(tspan))
#deltat = tspan[1]-tspan[0]
#for i in range(0,len(t)-1):
#    if i == 0:
#        x[i+1] = 0.5*usin[i]*(deltat)**2
#    else:
#        x[i+1] = x[i] +v[i]*deltat + 0.5*usin[i]*(deltat)**2
#        v[i+1] = (x[i+1]-x[i])/deltat

# Control with LQR
# Starting with just Q = R = I
Q = np.dot(np.transpose(Cn), Cn)
R = np.eye(m)
K, S, E = control.lqr(sys, Q, R)
print('Anything')

Ac = An-Bn*K

x0 = [x10, dx10, t10, dt10, t20, dt20]
tspan = np.linspace(0,20,1000)

sysc = control.ss(Ac,Bn,Cn,D)
tc, youtc, xoutc = control.impulse_response(sysc, tspan, x0, return_x=True)
plt.rc('text', usetex=True)
plt.plot(tc, xoutc[0], label=r'$x$')
plt.plot(tc, xoutc[2], label=r'$\theta_1$')
plt.plot(tc, xoutc[4], label=r'$\theta_2$')
plt.xlabel(r'Time, t')
plt.ylabel(r'Angle, $\theta$')
plt.legend(loc='best')

###############################################################################
# Animation stuff
###############################################################################
#f = open("textAn1.py", "w+")
#f.write("class BasicEquations(Scene):\n" + 
#"    def construct(self): \n" +
#"        eq1=TextMobject(\"$\\vec{X}_0 \\cdot \\vec{Y}_1 = 3$\") \n" +
#"        eq1.shift(2*UP) \n" +
#"        eq2=TexMobject(r\"\vec{F}_{net} = \sum_i \vec{F}_i\") \n" +
#"        eq2.shift(2*DOWN) \n" + 
#"        self.play(Write(eq1))  \n" +
#"        self.play(Write(eq2))")
#f.close()
#Qctext = print(vlatex(Qc))