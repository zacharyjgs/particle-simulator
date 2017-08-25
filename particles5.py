"""
This program is essentially a physics engine that
allows one to caluclate and visualize the behaviors
of interacting particles.
"""

#import various packages
from __future__ import division
import numpy as np
from visual import *

#sets up axes
xaxis = np.array([1,0,0])
yaxis = np.array([0,1,0])
zaxis = np.array([0,0,1])

#defines universal gravitation
def gravity(rad, m):
    a = -(G*m*rad)/(mag(rad)**3)
    return a

#defines hookes law/point of breakage
def tensile(rad, m, k, l, maxrad):
    radl = l*(rad/mag(rad))
    if (mag(rad) > maxrad) or (mag(rad) < l): a = [0.,0.,0.]
    else: a = (k*(radl - rad))/m
    return a

def compressive(rad, m, k, l, minrad):
    radl = l*(rad/mag(rad))
    if (mag(rad) < minrad) or (mag(rad) > l): a = [0.,0.,0.]
    else: a = (k*(radl - rad))/m
    return a

#finds a list of accelerations on each particle
def asum(p, x):
    k = []
    for b in range(len(p.l)):
        a = np.array([0., 0., 0.])
        if not p.l[b].fixed:
            for ob in range(len(p.l)):
                #if not p.l[b] is p.l[ob]:
                    #a += gravity(rad, m)
                if (b + 1) == ob or (b - 1) == ob:
                    rad = x[b] - x[ob]
                    m = p.m[b]
                    a += tensile(rad, m, 100000, 0.025, 100)
                    #a += compressive(rad, m, 5, 0.5, 0.2)
            a += g
        k.append(a)
    return np.array(k)

#runge kutta 4 integrator
def rk4(p, dt):
    x1 = np.array(p.x)
    v1 = np.array(p.v)
    k1 = asum(p, x1)
    x2 = x1 + v1*dt*0.5
    v2 = v1 + k1*dt*0.5
    k2 = asum(p, x2)
    x3 = x1 + v2*dt*0.5
    v3 = v1 + k2*dt*0.5
    k3 = asum(p, x3)
    x4 = x1 + v3*dt
    v4 = v1 + k3*dt
    k4 = asum(p, x4)
    xf = x1 + (dt/6)*(v1 + 2*v2 + 2*v3 + v4)/6
    vf = v1 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    #print "pos:"
    #print xf
    return xf, vf

#sets up particle class
class p(object):
    l = []; x = []; v = []; m = []; i = 0
    def __init__(self, x, v, m=0, fixed=False):
        self.fixed = fixed
        self.i = p.i
        p.l.append(self)
        p.x.append(x)
        p.v.append(v)
        p.m.append(m)
        p.i += 1

#makes spheres for each particle in vpython
def spheres(p):
    s = []
    for b in p.l: s.append(sphere(pos = p.x[b.i], radius = (p.m[b.i]**(1/3))/15))
    return s

#animates particles using time incrementation
def animate(p, dt, t0, tf):
    #creates spheres in vpython
    s = spheres(p)
    t = t0
    if tf < t0: dt = -dt; forwards = False
    else: forwards = True
    while (t < tf) == forwards:
        rk4vals = rk4(p, dt)
        p.x = rk4vals[0]
        p.v = rk4vals[1]
        p.x[fixed.i] = [0, 0.2*cos(t/25), 0.2*sin(t/25)]
        #updates spheres in vpython
        for i in range(len(s)):
            s[i].pos = p.x[i]
        t += dt/k

#constants
G = 6.67e-11
g = [0.,0.,-9.8]
k = 0.0172020989
#1.990983669e-7

#creates a set of particles
fixed = p([0.,0.2,0.], [0.,0.,0.], fixed=True)
for i in range(1,20):
    p([0,0.2,-i/40], [0.,0.,0.], m=0.1)
p([0,0.2,-0.5], [0.,0.,0.], m=0.4)

#sets up vpython display
scene = display(width = 900, height = 900, up = zaxis, forward = -xaxis)
arrow(pos = (0,0,0), axis = xaxis, shaftwidth = 0.01, color = color.red)
arrow(pos = (0,0,0), axis = yaxis, shaftwidth = 0.01, color = color.green)
arrow(pos = (0,0,0), axis = zaxis, shaftwidth = 0.01, color = color.blue)

animate(p, 0.0015, 0, 10000)
