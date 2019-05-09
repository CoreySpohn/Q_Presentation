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
from sympy.physics.vector import time_derivative as dt
from sympy.physics.vector import dot
from sympy.physics.vector.printing import vpprint, vlatex
import numpy as np
import control
import matplotlib.pyplot as plt
physics.vector.init_vprinting()

###############################################################################
# DYNAMICS
###############################################################################
# Variable creation
l1, l2, m1, m2, M, g = symbols('l_1 l_2 m_1 m_2 M g')
t1, t2, x = vector.dynamicsymbols('theta_1 theta_2 x')

# Constant values for substitution later
l1n = 5
l2n = 6
m1n = 1
m2n = 1
Mn  = 2
gn  = 1
t10 = 1 # Initial angle
dt10 = 0
t20 = 1
dt20 = 0

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
T = (1/2)*M*dot(dt(x*N.x,N),dt(x*N.x,N)) + (1/2)*m1*dot(v1,v1) + (1/2)*m2*dot(v2,v2)
U = M*g*(l1*cos(t1)+l2*cos(t2))

# Lagrangian
L = T-U

# Equations of motion
eqx = dt(L.diff(dt(x,N)),N)-L.diff(x)
eqt1 = dt(L.diff(dt(t1,N)),N)-L.diff(t1)
eqt2 = dt(L.diff(dt(t2,N)),N)-L.diff(t2)

# Solve for variables
xdd = solve(eqx,dt(x,N,2))
t1dd = solve(eqt1,dt(t1,N,2))
t2dd = solve(eqt2,dt(t2,N,2))

# Linearize
xdd = xdd[0].subs([(sin(t1),t1), (cos(t1),1), (sin(t2),t2), (cos(t2), 1)])
t1dd = t1dd[0].subs([(sin(t1),t1), (cos(t1),1), (sin(t2),t2), (cos(t2), 1)])
t2dd = t2dd[0].subs([(sin(t1),t1), (cos(t1),1), (sin(t2),t2), (cos(t2), 1)])

###############################################################################
# CONTROLS
###############################################################################
state = [t1, dt(t1,N), t2, dt(t2,N)]
stated = [dt(t1,N), t1dd, dt(t2,N), t2dd]
outputs = [t1, t2]
inputs = [dt(x,N,2)]

# A matrix
avals = []
for i in range(0,len(stated)):
    for j in range(0, len(state)):
        avals.append(stated[i].diff(state[j]))
A = Matrix(len(stated), len(stated), avals)
    
# B matrix
bvals = []
for i in range(0,len(stated)):
    for j in range(0, len(inputs)):
        bvals.append(stated[i].diff(inputs[j]))
B = Matrix(len(stated), len(inputs), bvals)

# C matrix
cvals = []
for i in range(0,len(outputs)):
    for j in range(0, len(state)):
        cvals.append(outputs[i].diff(state[j]))
C = Matrix(len(outputs), len(state), cvals)

# D matrix
D = np.zeros([int(C.shape[0]),1])

# Controllability matrix
# This is messier than I'd like but it works
Qi = []
for i in range(0,len(state)):
    Qi.append((A**i * B))

Qcvals = []
for i in range(0,len(state)):
    for j in range(0, len(state)):
        Qcvals.append(Qi[j][i])

Qc = Matrix(len(state), len(state), Qcvals)


# Numerical values for the matrices
nSubs = [(l1,l1n), (l2,l2n), (m1,m1n), (m2,m2n), (M,Mn), (g,gn)]

An = np.array(A.subs(nSubs))
Bn = np.array(B.subs(nSubs))
Cn = np.array(C.subs(nSubs))
Qcn = np.array(Qc.subs(nSubs), dtype='float')
# Stability check
if np.linalg.matrix_rank(Qcn) == len(state):
    print('System is controllable\n')
else:
    print('System is uncontrollable\n')
vector.vprint(Qcn)

# Control system
ss = control.ss(An,Bn,Cn,D)
x0 = [t10, dt10, t20, dt20]
tspan = np.linspace(0,2,1000)


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
usin = 15*np.cos(tspan*2*np.pi)
t, yout, xout = control.forced_response(ss, tspan, usin, x0)
plt.rc('text', usetex=True)
plt.plot(t,yout[0], label=r'$\theta_1$')
plt.plot(t,yout[1], label=r'$\theta_2$')
plt.xlabel(r'Time, t')
plt.ylabel(r'Angle, $\theta$')
plt.legend(loc='best')

# Calculate position x
x = np.zeros(len(tspan))
v = np.zeros(len(tspan))
deltat = tspan[1]-tspan[0]
for i in range(0,len(t)-1):
    if i == 0:
        x[i+1] = 0.5*usin[i]*(deltat)**2
    else:
        x[i+1] = x[i] +v[i]*deltat + 0.5*usin[i]*(deltat)**2
        v[i+1] = (x[i+1]-x[i])/deltat

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