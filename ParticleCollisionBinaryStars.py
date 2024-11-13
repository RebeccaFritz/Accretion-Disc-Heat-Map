# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:45:08 2024

@author: rgfri
"""

import vpython as vp

autoscale = 1

# Define Globals
G  = 6.67E-11 # gravitational constant
P = ... # orbital period
AMBIENT_COLOR # color of particles when ambient temp 
""" Eventually:
HOT_COLOR # color of particles when hot
COLD_COLOR # color of particles when cold
"""
h = # time step
L1 # distance

# create Objects
starA """ If we don't want to animate the stars we can just put their mass and position and such in variables"""
starB
# legend """ Eventually """
starA.trail """ do we want trails"""
starB.trail
accretionDisk = """ do we still want to do this? could do two transparent cylinders compounded together"""


# a function to create particles
## vector, vector --> sphere
""" I already created this """

# ancillaries
particles = [] # list of particles
t = 0 # time

# a function to determine the acceleration of a star
## ... --> ...
def accStar(...):
    return

# a function to run the Runge-Kutta alg on a star
## ... --> ...
def rkStar(...):
    return

# a function to check if the particles have collided
## ... --> ...
def particleCol(...):
    return

# a function to run the Runge-Kutta alg on particles
"""
can do all the particles here like in Dwarf Nova code a 
single particle here and loop in the animation we can do 
it whichever way we prefer
"""
## ... --> ...
def rkParticle(...):
    return

""" Eventually:
# a function to check particle temperature and change their color
## ... --> ...
def particleTemp(...):
    return
"""

# run animation
while True:
    rate(...)