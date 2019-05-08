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
init_vprinting()

###############################################################################
# DYNAMICS
###############################################################################
# Variable creation
l1, l2, m1, m2, M, g = symbols('l_1 l_2 m_1 m_2 M g')
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

