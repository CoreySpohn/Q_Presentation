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
l1n = 2
l2n = 1
m1n = 0.1
m2n = 0.1
Mn  = 1
gn  = 10
# Initial conditions for initial response
t10 = 1 # Initial angle
dt10 = 3
t20 = -5
dt20 = 0
x10 = 0
dx10 = 0
nSubs = [(l1,l1n), (l2,l2n), (m1,m1n), (m2,m2n), (M,Mn), (g,gn)]

An = np.array(A.subs(nSubs)).astype(float)
Bn = np.array(B.subs(nSubs)).astype(float)
Cn = np.array(C.subs(nSubs)).astype(float)
Qcn = np.array(Qc.subs(nSubs), dtype='float')
Qon = np.array(Qo.subs(nSubs), dtype='float')

# Stability check
if np.linalg.matrix_rank(Qcn) == n:
    print('System is controllable\n')
else:
    print('System is uncontrollable\n')


# Control system
sys = control.ss(An,Bn,Cn,D)
x0 = [x10, dx10, t10, dt10, t20, dt20]
x0_stable = [0, 0, 0, 0, 0, 0]
tspan = np.linspace(0,5,1000)


# Control with LQR
# Starting with just Q = R = I
rho = 100
Q = rho*np.diag([0.1, 0, 1, 0, 1, 0])
# Matrix(np.diag([0.1, 0, 1, 0, 1, 0])).subs([(0.0, 0), (1.0, 1)])
R = np.eye(m)
K, S, E = control.lqr(sys, Q, R)


Ac = An-Bn*K

x0 = [x10, dx10, t10, dt10, t20, dt20]
tspan = np.linspace(0,10,1000)

sysc = control.ss(Ac,Bn,Cn,D) 


#xlims = [0, 10]
condition = '$l_1 = 2 l_2$'
# Initial Response
tcic, youtic, xoutic = control.initial_response(sysc, tspan, x0, return_x=True)

plt.style.use('dark_background')
fig0, ax0 = plt.subplots()
plt.rc('text', usetex=True)
ax0.plot(tcic, xoutic[2], label=r'$\theta_1$', color='#39ff14')
ax0.plot(tcic, xoutic[4], label=r'$\theta_2$', color='#4DA4B5')
ax0.set_xlim(right=xlims[1])
ax0.set_xlabel(r'Time, seconds')
ax0.set_ylabel(r'$\theta$, degrees')
impulsetitle = condition + r' Initial Response'
ax0.set_title(impulsetitle)
ax0.legend(loc='best')
fig0.savefig(condition+'initial_response.png', dpi=300)

# Impulse Response
tcimpulse, youtimpulse, xoutimpulse = control.impulse_response(sysc, tspan, x0_stable, return_x=True)
plt.style.use('dark_background')
fig1, ax1 = plt.subplots()
plt.rc('text', usetex=True)
ax1.plot(tcimpulse, xoutimpulse[2], label=r'$\theta_1$', color='#39ff14')
ax1.plot(tcimpulse, xoutimpulse[4], label=r'$\theta_2$', color='#4DA4B5')
ax1.set_xlim(right=xlims[1])
ax1.set_xlabel(r'Time, seconds')
ax1.set_ylabel(r'$\theta$, degrees')
impulsetitle = condition+ ' Impulse Response'
ax1.set_title(impulsetitle, usetex=True)
ax1.legend(loc='best')
fig1.savefig(condition+'impulse_response.png', dpi=300)

# Step Response
#tcstep, youtstep, xoutstep = control.step_response(sysc, tspan, x0_stable, return_x=True)
##plt.style.use('dark_background')
#fig2, ax2 = plt.subplots()
#plt.rc('text', usetex=True)
#ax2.plot(tcimpulse, xoutstep[2], label=r'$\theta_1$', color='#39ff14')
#ax2.plot(tcimpulse, xoutstep[4], label=r'$\theta_2$', color='#4DA4B5')
#ax2.set_xlim(right=xlims[1])
#ax2.set_xlabel(r'Time, seconds')
#ax2.set_ylabel(r'$\theta$, degrees')
#impulsetitle = condition+' Step Response'
#ax2.set_title(impulsetitle)
#ax2.legend(loc='best')
#fig2.savefig(condition+'step_response.png', dpi=300)

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