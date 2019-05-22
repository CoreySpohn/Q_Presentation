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
import matplotlib.animation as animation

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
l1n = 2 # 
l2n = 1.8
m1n = 0.1
m2n = 0.1
Mn  = 1 # kg
gn  = 10 # m/s^2
# Initial conditions for initial response
t10 = 2.5 # Initial angle in degrees
dt10 = 0
t20 = 2.5
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
Q = rho*np.diag([0.01, 0, 1, 0, 1, 0])
# Matrix(np.diag([0.01, 0, 1, 0, 1, 0])).subs([(0.0, 0), (1.0, 1)])
R = np.eye(m)
K, S, E = control.lqr(sys, Q, R)


Ac = An-Bn*K

x0 = [x10, dx10, t10, dt10, t20, dt20]
tspan = np.linspace(0,15,1000)

sysc = control.ss(Ac,Bn,Cn,D) 


#xlims = [0, 10]
l2l1ratio = l2n/l1n
m2m1ratio = m2n/m1n
if l2l1ratio == 1:
    l2l1ratiostr= ''
else:
    l2l1ratiotr = str(round(l2l1ratio,1))
if m2m1ratio == 1:
    m2m1ratiostr= ''
else:
    m2m1ratiostr = str(round(m2m1ratio,1))
    
condition = '$l_2 = ' + l2l1ratiotr + 'l_1, m_1=' + m2m1ratiostr + 'm_2$'
filecondition = condition.replace('$','').replace('_', '').replace(',', '_').replace(' ', '')

# Colors
t1color = '#39ff14'
t2color = '#4DA4B5'
xcolor = '#F9C80E'

# Initial Response
tcic, youtic, xoutic = control.initial_response(sysc, tspan, x0, return_x=True)

plt.style.use('dark_background')
fig0, ax0 = plt.subplots()
ax0 = plt.subplot(1,1,1)
plt.rc('text', usetex=True)
ax0.plot(tcic, xoutic[2], label=r'$\theta_1$', color=t1color)
ax0.plot(tcic, xoutic[4], label=r'$\theta_2$', color=t2color)
ax0.set_xlabel(r'Time, seconds')
ax0.set_ylabel(r'$\theta$, degrees')
ictitle = condition + r' Initial Response'
ax0.set_title(ictitle)
ax0.legend(loc='best')
# With x
ax0R = ax0.twinx()
ax0R.yaxis.tick_right()
ax0R.yaxis.set_label_position('right')
ax0R.set_ylabel('Position')
ax0R.plot(tcic, xoutic[0], label=r'x', color=xcolor)

#fig0.savefig(condition+'initial_response.png', dpi=300)

# Impulse Response
timpulse, youtimpulse, xoutimpulse = control.impulse_response(sysc, tspan, x0_stable, return_x=True)
plt.style.use('dark_background')
fig1, ax1 = plt.subplots()
plt.rc('text', usetex=True)
ax1.plot(timpulse, xoutimpulse[2], label=r'$\theta_1$', color=t1color)
ax1.plot(timpulse, xoutimpulse[4], label=r'$\theta_2$', color=t2color)

ax1.set_xlabel(r'Time, seconds')
ax1.set_ylabel(r'$\theta$, degrees')
steptitle = condition+ ' step Response'
ax1.set_title(steptitle, usetex=True)
ax1.legend(loc='best')
#fig1.savefig(condition+'step_response.png', dpi=300)

####
# Animations
###
fileloc = 'media/plots/'
dpinum = 200
fpsnum = 200
codecname = 'libx264'
bitratenum = 1000

def align_yaxis_np(ax1, ax2):
    """Align zeros of the two axes, zooming them out by same ratio"""
    axes = np.array([ax1, ax2])
    extrema = np.array([ax.get_ylim() for ax in axes])
    tops = extrema[:,1] / (extrema[:,1] - extrema[:,0])
    # Ensure that plots (intervals) are ordered bottom to top:
    if tops[0] > tops[1]:
        axes, extrema, tops = [a[::-1] for a in (axes, extrema, tops)]

    # How much would the plot overflow if we kept current zoom levels?
    tot_span = tops[1] + 1 - tops[0]

    extrema[0,1] = extrema[0,0] + tot_span * (extrema[0,1] - extrema[0,0])
    extrema[1,0] = extrema[1,1] + tot_span * (extrema[1,0] - extrema[1,1])
    [axes[i].set_ylim(*extrema[i]) for i in range(2)]    

##
# Initial condition response with x
##
tic = tcic
xic = xoutic[0]
t1ic = xoutic[2]
t2ic = xoutic[4]
fig3, ax3 = plt.subplots()
ax3.set_xlabel(r'Time')
ax3.set_ylabel(r'$\theta$, degrees')
ax3.set_title(r'')
plt.rc('text', usetex=True)
t1icline, = ax3.plot(tcic,xoutic[2], label=r'$\theta_1$', color=t1color)
t2icline, = ax3.plot(tcic,xoutic[4], label=r'$\theta_2$', color=t2color)
ax3.legend(loc='upper right')
iclims = [ax0.get_xlim()[0], ax0.get_xlim()[1], ax0.get_ylim()[0], ax0.get_ylim()[1]]
# x stuff
ax3R = ax3.twinx()
ax3R.yaxis.tick_right()
ax3R.yaxis.set_label_position('right')
ax3R.set_ylabel('Position, m')
xicline, = ax3R.plot(tcic, xic, label = r'x', color=xcolor)
ax3Rlims = [ax0R.get_xlim()[0], ax0R.get_xlim()[1], ax0R.get_ylim()[0], ax0R.get_ylim()[1]]
# legend
lines = [t1icline, t2icline, xicline]
labels = [l.get_label() for l in lines]
ax3.legend(lines, labels, loc = 'upper right', prop={'size': 13})
ictitle = condition + r' Initial Response'
ax3.set_title(ictitle)
# align the zeros
#align_yaxis(ax3, 0, ax3R, 0)
align_yaxis_np(ax3, ax3R)

ax3Rlims = [ax3R.get_xlim()[0], ax3R.get_xlim()[1], ax3R.get_ylim()[0], ax3R.get_ylim()[1]]
xicline.axes.axis(ax3Rlims)
fig3.tight_layout()
plt.show()
fig3.savefig(fileloc+'initial_response_wx_' + filecondition+'.png', dpi=300)
#def updateic(num, tic, xic, xicline):
#    xicline.set_data(tic[:num], xic[:num])
#    xicline.axes.axis(ax3Rlims)
#    return xicline,
#
#
#ani3 = animation.FuncAnimation(fig3, updateic, len(tic), fargs=[tic, xic, xicline],
#                              interval=5, blit=True)
#ani3.save(fileloc+'initial_response_wx_' + filecondition + '.mp4', writer='ffmpeg', dpi=dpinum, extra_args=['-pix_fmt', 'yuv420p'])

##
# Initial condition response w/o x
##
tic = tcic
t1ic = xoutic[2]
t2ic = xoutic[4]
fig2, ax2 = plt.subplots()
ax2.set_xlabel(r'Time, seconds')
ax2.set_ylabel(r'$\theta$, degrees')
plt.rc('text', usetex=True)
t1icline, = ax2.plot(tcic,t1ic, label=r'$\theta_1$', color = '#39ff14')
t2icline, = ax2.plot(tcic,t2ic, label=r'$\theta_2$', color = '#4DA4B5')
ax2.legend(loc='upper right', prop={'size': 13})
iclims = [ax0.get_xlim()[0], ax0.get_xlim()[1], ax3.get_ylim()[0], ax3.get_ylim()[1]]
t1icline.axes.axis(iclims)
ictitle = condition + r' Initial Response'
ax2.set_title(ictitle)
fig2.tight_layout()

fig2.savefig(fileloc+'initial_response_' + filecondition+'.png', dpi=300)

#def updateic(num, tic, t1ic, t2ic, t1icline, t2icline):
#    t1icline.set_data(tic[:num], t1ic[:num])
#    t2icline.set_data(tic[:num], t2ic[:num])
#    t1icline.axes.axis(iclims)
#    return t1icline, t2icline,
#
#ani2 = animation.FuncAnimation(fig2, updateic, len(tic), fargs=[tic, t1ic, t2ic, t1icline, t2icline],
#                              interval=5, blit=True)
#ani2.save(fileloc+'initial_response_' + filecondition + '.mp4', writer='ffmpeg', dpi=dpinum, extra_args=['-pix_fmt', 'yuv420p'])


##
# Impulse response with x
##
ximpulse = xoutimpulse[0]
t1impulse = xoutimpulse[2]
t2impulse = xoutimpulse[4]
fig5, ax5 = plt.subplots()
ax5.set_xlabel(r'Time')
ax5.set_ylabel(r'$\theta$, degrees')
ax5.set_title(r'')
plt.rc('text', usetex=True)
t1impulseline, = ax5.plot(timpulse,xoutimpulse[2], label=r'$\theta_1$', color=t1color)
t2impulseline, = ax5.plot(timpulse,xoutimpulse[4], label=r'$\theta_2$', color=t2color)
ax5.legend(loc='upper right')
impulselims = [ax1.get_xlim()[0], ax1.get_xlim()[1], ax1.get_ylim()[0], ax1.get_ylim()[1]]
# x stuff
ax5R = ax5.twinx()
ax5R.yaxis.tick_right()
ax5R.yaxis.set_label_position('right')
ax5R.set_ylabel('Position, m')
ximpulseline, = ax5R.plot(timpulse, ximpulse, label = r'x', color=xcolor)
ax5Rlims = [ax1.get_xlim()[0], ax1.get_xlim()[1], ax1.get_ylim()[0], ax1.get_ylim()[1]]
# legend
lines = [t1impulseline, t2impulseline, ximpulseline]
labels = [l.get_label() for l in lines]
ax5.legend(lines, labels, loc = 'upper right', prop={'size': 13})
impulsetitle = condition + r' Impulse Response'
ax5.set_title(impulsetitle)
# align the zeros
#align_yaxis(ax5, 0, ax5R, 0)
align_yaxis_np(ax5, ax5R)
ax5Rlims = [ax5R.get_xlim()[0], ax5R.get_xlim()[1], ax5R.get_ylim()[0], ax5R.get_ylim()[1]]
ximpulseline.axes.axis(ax5Rlims)
fig5.tight_layout()
fig5.savefig(fileloc+'impulse_response_wx_' + filecondition + '.png', dpi=300)
#def updateimpulse(num, timpulse, ximpulse, ximpulseline):
#    ximpulseline.set_data(timpulse[:num], ximpulse[:num])
#    ximpulseline.axes.axis(ax5Rlims)
#    return ximpulseline,
#
#ani5 = animation.FuncAnimation(fig5, updateimpulse, len(tic), fargs=[timpulse, ximpulse, ximpulseline],
#                              interval=5, blit=True)
#
#ani5.save(fileloc+'impulse_response_wx_' + filecondition + '.mp4', writer='ffmpeg', codec='h264', dpi=dpinum)

##
# Impulse response w/o x
##
t1impulse = xoutimpulse[2]
t2impulse = xoutimpulse[4]
fig4, ax4 = plt.subplots()
ax4.set_xlabel(r'Time')
ax4.set_ylabel(r'$\theta$, degrees')
plt.rc('text', usetex=True)
t1impulseline, = ax4.plot(tcic, t1impulse, label=r'$\theta_1$', color=t1color)
t2impulseline, = ax4.plot(tcic, t2impulse, label=r'$\theta_2$', color=t2color)
ax4.legend(loc='upper right', prop={'size': 13})
impulselims = [ax1.get_xlim()[0], ax1.get_xlim()[1], ax5.get_ylim()[0], ax5.get_ylim()[1]]
t1impulseline.axes.axis(impulselims)
impulsetitle = condition + r' Impulse Response'
ax4.set_title(impulsetitle)
fig4.tight_layout()
fig4.savefig(fileloc+'impulse_response_' + filecondition + '.png', dpi=300)
#def updateic(num, timpulse, t1impulse, t2impulse, t1impulseline, t2impulseline):
#    t1impulseline.set_data(timpulse[:num], t1impulse[:num])
#    t2impulseline.set_data(timpulse[:num], t2impulse[:num])
#    t1impulseline.axes.axis(impulselims)
#    return t1impulseline, t2impulseline,
#ani4 = animation.FuncAnimation(fig4, updateic, len(tic), fargs=[timpulse, t1impulse, t2impulse, t1impulseline, t2impulseline],
#                              interval=5, blit=True)
#
#ani4.save(fileloc+'impulse_response_' + filecondition + '.mp4', writer='ffmpeg', codec='h264', dpi=dpinum)







# Step response with x
#tic = tcic
#xic = xoutic[0]
#
#fig2, ax2 = plt.subplots()
#ax2.set_xlabel(r'Time')
#ax2.set_ylabel(r'$\theta$, degrees')
#ax2.set_title(r'')
#plt.rc('text', usetex=True)
#t1icline, = ax2.plot(tcic,xoutic[2], label=r'$\theta_1$', color=t1color)
#t2icline, = ax2.plot(tcic,xoutic[4], label=r'$\theta_2$', color=t2color)
#ax2.legend(loc='upper right')
#iclims = [ax0.get_xlim()[0], ax0.get_xlim()[1], ax0.get_ylim()[0], ax0.get_ylim()[1]]
#
## x stuff
#ax2R = ax2.twinx()
#ax2R.yaxis.tick_right()
#ax2R.yaxis.set_label_position('right')
#ax2R.set_ylabel('Position')
#xicline, = ax2R.plot(tcic, xic, label = r'x', color=xcolor)
#ax2Rlims = [ax0R.get_xlim()[0], ax0R.get_xlim()[1], ax0R.get_ylim()[0], ax0R.get_ylim()[1]]
## legend
#lines = [t1icline, t2icline, xicline]
#labels = [l.get_label() for l in lines]
#ax2.legend(lines, labels, loc = 'upper right', prop={'size': 13})
#steptitle = condition + r' Step Response'
#ax2.set_title(steptitle)
#
## align the zeros
#def align_yaxis(ax1, v1, ax2, v2):
#    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
#    _, y1 = ax1.transData.transform((0, v1))
#    _, y2 = ax2.transData.transform((0, v2))
#    inv = ax2.transData.inverted()
#    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
#    miny, maxy = ax2.get_ylim()
#    ax2.set_ylim(miny+dy, maxy+dy)
#
#align_yaxis(ax2, 0, ax2R, 0)
#ax2Rlims = [ax2R.get_xlim()[0], ax2R.get_xlim()[1], ax2R.get_ylim()[0], ax2R.get_ylim()[1]]
#def updateic(num, tic, xic, xicline):
#    xicline.set_data(tic[:num], xic[:num])
#    xicline.axes.axis(ax2Rlims)
#    return xicline,
#
#ani = animation.FuncAnimation(fig2, updateic, len(tic), fargs=[tic, xic, xicline],
#                              interval=5, blit=True)
#ani.save(condition+'initial_response_wx.mp4', writer='ffmpeg', codec='h264', dpi=300)


#plt.show()

# Step Response
#tcstep, youtstep, xoutstep = control.step_response(sysc, tspan, x0_stable, return_x=True)
##plt.style.use('dark_background')
#fig2, ax2 = plt.subplots()
#plt.rc('text', usetex=True)
#ax2.plot(tcstep, xoutstep[2], label=r'$\theta_1$', color=t1color)
#ax2.plot(tcstep, xoutstep[4], label=r'$\theta_2$', color=t2color)
#ax2.set_xlim(right=xlims[1])
#ax2.set_xlabel(r'Time, seconds')
#ax2.set_ylabel(r'$\theta$, degrees')
#steptitle = condition+' Step Response'
#ax2.set_title(steptitle)
#ax2.legend(loc='best')
#fig2.savefig(condition+'step_response.png', dpi=300)

