# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 18:44:56 2024

@author: rgfri
"""

import vpython as vp
from random import randrange 

#function to create particles 
# vec, vec --> sphere
def create_particle(pos, vel):
    return vp.sphere(pos = pos, vel = vel, mass = 1, radius = 1, color = vp.color.white)


#create 1000 particles and add them to a list
particles = []
for i in range(0, 1000):
    x = randrange(0,100) # pick random number
    y = randrange(0,100)
    z = randrange(0,100)
    vx = 0 
    vy = 0
    vz = 0
    particles.append(create_particle(vp.vec(x,y,z), vp.vec(vx,vy,vz)))
    
particles[3].color = vp.color.red # proof of concept on how to access a particle

# ancillaries
vp.scene.autoscale = 1

# animate
while True:
    vp.rate(800)