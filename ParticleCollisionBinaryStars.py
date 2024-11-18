Web VPython 3.2
from random import randrange

# create canvases
scene.autoscale = 1
scene = canvas(title = 'Particle Collisions in a Binary Star System', color = color.black, align='left')

# Define Globals
G  = 6.67E-11 # gravitational constant
AU = 1.5E11 # astronomical unit
YEAR = 365.25*24*60*60 # year in seconds
MS = 1.988400e20 # mass of the sun 
P = 7.8719E3 # orbital period                  # from sample code
MID_COLOR = vec(2.50, 1.50, 0.50) # color of particles when ambient temp 
HIGH_COLOR = vec(2.00, 0.00, 0.00) # color of particles when hot
LOW_COLOR = vec(2.50, 2.00, 0.00) # color of particles when cold
STARA_COLOR = vec(2.55, 2.50, 2.00)
STARB_COLOR = vec(2.55, 2.55, 2.55)

h = 2.0 # time step                            # from sample code
L1 = 3.69E8 # distance to Lagrange Point 1            # from sample code

# create scene objects
starA = sphere(pos=vec(0,0,0), radius= 6.6E6, mass=1.5E30, color=vec(0,1,0), visible=True) # numbers from sample code
starB = sphere(pos=vec(5.80E8,0,0), radius= 1.5E8, mass=3.75E29, color=vec(0,0.5,0), visible=True) # numbers from sample code
# legend = # Eventually 
starA.trail = curve(pos=starA.pos, color=starA.color) # do we want trails
starB.trail = curve(pos=starB.pos, color=starB.color)
# accretionDisk = # do we still want to do this? could do two transparent cylinders compounded together


# a function to create particles
## vector, vector --> sphere
def create_particle(pos, vel):
    return sphere(pos = pos, vel = vel, mass = 2.0e10, radius = 4E6, color = MID_COLOR) # numbers from sample code

# ancillaries
t = 0 # time

#create 1000 particles and add them to a list
particles = []
for i in range(0, 1000):
    x = randrange(0,100)*4E6 # pick random number                        # need to update these
    y = randrange(0,100)*4E6
    z = randrange(0,100)*4E6
    vx = 0                                                               # need to update these 
    vy = 0
    vz = 0
    particles.append(create_particle(vec(x,y,z), vec(vx,vy,vz)))

# a function to determine the acceleration of a star
## ... --> ...
"""
def accStar(...):
    return
"""

# a function to run the Runge-Kutta alg on a star
## ... --> ...
"""
def rkStar(...):
    return
"""

# a function to check if the particles have collided
## ... --> ...
"""
def particleCol(...):
    return
"""

# a function to run the Runge-Kutta alg on particles
"""
can do all the particles here like in Dwarf Nova code a 
single particle here and loop in the animation we can do 
it whichever way we prefer
"""
## ... --> ...
"""
def rkParticle(...):
    return
"""

""" Eventually:
# a function to check particle temperature and change their color
## ... --> ...
def particleTemp(...):
    return
"""

# run animation
while True:
    rate(30) # allows the program to look the same on every computer
