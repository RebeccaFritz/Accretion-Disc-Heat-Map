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
u = 1.66e-27 # atomic mass in kg
H_MASS = 1.01*u # mass of a hydrogen atom
MID_COLOR = vec(0.98, 0.59, 0.20) # color of particles when ambient temp : rgb(250,150, 0.20)
HIGH_COLOR = vec(0.78, 0.00, 0.00) # color of particles when hot: rgb(200, 000, 000)
LOW_COLOR = vec(0.98, 0.78, 0.00) # color of particles when cold: rgb(250, 200, 000)
STARA_COLOR = vec(0.98, 0.98, 0.78) #rgb(250,250,200)
STARB_COLOR = vec(0.98, 0.94, 0.59 ) #rgb(250,240,150)
P_RADIUS = 4e6 # particle radius

h = 2.0 # time step                            # from sample code
L1 = 3.69E8 # distance to Lagrange Point 1            # from sample code
L2= 2.11e8                                             # from sample code

# create scene objects
starA = sphere(pos=vec(0,0,0), radius= 6.6E6, mass=1.5E30, color=STARA_COLOR, visible=True) # numbers from sample code
starB = sphere(pos=vec(5.80E8,0,0), radius= 1.5E8, mass=3.75E29, color=STARB_COLOR, visible=True) # numbers from sample code
# legend = # Eventually 
starA.trail = curve(pos=starA.pos, color=starA.color) # do we want trails
starB.trail = curve(pos=starB.pos, color=starB.color)
# accretionDisk = # do we still want to do this? could do two transparent cylinders compounded together

# set initial velocities
starB.vel = 4.0*vector(0,+6.78e4,0)             # from sample code
starA.vel = 4.0*vector(0,-6.78e4,0)*0.25        # from sample code

# ancillaries
t = 0 # time
particle_list = []

# a function to determine the acceleration of a star
## ... --> ...
def accStar(starpos):
    return -G*starother.mass*(starpos-starother.pos)/r**3

# a function to run the Runge-Kutta alg on a star
## ... --> ...
# borrowed from sample code
def rkStar(star):
    k1v = h*accStar(star.pos)
    k1x = h*star.vel

    k2v = h*accStar(star.pos + k1x/2.0)
    k2x = h*(star.vel + k1v/2.0)

    k3v = h*accStar(star.pos + k2x/2.0)
    k3x = h*(star.vel + k2v/2.0)

    k4v = h*accStar(star.pos + k3x)
    k4x = h*(star.vel + k3v)

    star.vel += (k1v + 2.0*k2v + 2.0*k3v + k4v)/6.0
    star.pos += (k1x + 2.0*k2x + 2.0*k3x + k4x)/6.0
    

# a function to remove particles from the particle_list and to increase starB's mass
## ... --> ...
"""
def removeParticle(...):
    return
"""

# a function to add particles to the particle_list and to decrease starA's mass
## void --> void
def AddParticle():
    starA.mass = starA.mass - H_MASS*1e3 
    # This may not be how we want to create particles
    Omega = 2.0*pi/P                                                            # from sample code
    for i in range(0, 10):
        particle = sphere(mass = H_MASS*1e3, radius = P_RADIUS, color = MID_COLOR)
        particle.pos = starA.pos +0.6367*(starB.pos-starA.pos)                      # from sample code
        # set initial velocities
        particle.vel = vector(0,0,0)                                           # from sample code
        particle.vel.x = -particle.pos.x*Omega*(1.0 + random())                # from sample code
        particle.vel.y = particle.pos.y*Omega*(1.0 + random())                 # from sample code
        particle.vel.z = particle.pos.z*Omega*(1.0 + random())                 # from sample code
        ##
        particle_list.append(particle)
        
#Final velocities of particles in collisions
# int, int --> void
def velfinal(particleidx, otheridx):
    v1f = ((particle_list[particleidx].mass - particle_list[otheridx].mass) / (particle_list[particleidx].mass + particle_list[otheridx].mass))*particle_list[particleidx].vel \
        + ((2*particle_list[otheridx].mass) / (particle_list[particleidx].mass + particle_list[otheridx].mass))*particle_list[otheridx].vel
    particle_list[particleidx].vel = v1f
    
    #other particle final vel
    v2f = ((particle_list[otheridx].mass - particle_list[particleidx].mass) / (particle_list[particleidx].mass + particle_list[otheridx].mass))*particle_list[otheridx].vel \
        + ((2*particle_list[particleidx].mass) / (particle_list[particleidx].mass + particle_list[otheridx].mass))*particle_list[particleidx].vel
    particle_list[otheridx].vel = v2f
    return

# a function to check if the particles have collided
## ... --> ...
def particleCol():
    particleDiameter = P_RADIUS*2
    collided = [] # list of indicies of the particles for which the collisions have been calculated
    
    # this does check each particle pair for collisions, but it is so slow it breaks the animation
    for i in range(0, len(particle_list)):         
        if i not in collided: # make sure that the particle at this index has not already had a collision run on it
            particle1 = particle_list[i]
        
            collisions = [] # list of indicies of currently colliding particles  
            collisions.append(i) # the index of colliding particle1
            for j in range(i+1, len(particle_list)):
                particle2 = particle_list[j]
                    
                if (particle2.pos.x - particleDiameter <= particle1.pos.x <= particle2.pos.x + particleDiameter) \
                        and (particle2.pos.y - particleDiameter <= particle1.pos.y <= particle2.pos.y + particleDiameter) \
                            and (particle2.pos.z - particleDiameter <= particle1.pos.z <= particle2.pos.z + particleDiameter):
                    collisions.append(j) # the index of colliding particle2
                    
            # collisions now has the indicies of all the particles colliding with particle1
            if len(collisions) > 1:
                # calculate particle collision for all particles at the indicies in the collisions list
                
                #velocity change in elastic particle collisions
                for k in range(0,len(collisions), 2):
                    velfinal(collisions[k], collisions[k+1])
                
                # when done colliding add collisions indicies to collided indices list
                collided.extend(collisions) 
        i += 1
        
    


# a function to run the Runge-Kutta alg on particles       # from the sample code
## void --> void
def rkParticles():
    for particle in particle_list:  
        rA = mag(particle.pos - starA.pos)
        rB = mag(particle.pos - starB.pos)
        if rA > L1/100 and rB > L2/100:
            acc_p = G*starA.mass*(starA.pos - particle.pos)/rA**3 + G*starB.mass*(starB.pos - particle.pos)/rB**3
            k1v = h*acc_p
            k1x = h*particle.vel
            k2v = h*(G*starA.mass*(starA.pos - (particle.pos+k1x/2.0))/rA**3 + G*starB.mass*(starB.pos - (particle.pos+k1x/2.0))/rB**3)
            k2x = h*(particle.vel + k1v/2.0)
            k3v = h*(G*starA.mass*(starA.pos - (particle.pos+k2x/2.0))/rA**3 + G*starB.mass*(starB.pos - (particle.pos+k2x/2.0))/rB**3)
            k3x = h*(particle.vel + k2v/2.0)
            k4v = h*(G*starA.mass*(starA.pos - (particle.pos+k3x/2.0))/rA**3 + G*starB.mass*(starB.pos - (particle.pos+k3x/2.0))/rB**3)            
            k4x = h*(particle.vel + k3v)
            particle.vel += (k1v + 2.0*k2v + 2.0*k3v + k4v)/6.0
            particle.pos += (k1x + 2.0*k2x + 2.0*k3x + k4x)/6.0
        



""" Eventually:
# a function to check particle temperature and change their color
## ... --> ...
def particleTemp(...):
    return
"""

done = 0

# run animation
while True:
    r = mag(starA.pos - starB.pos)
    
    # run Runge-Kutta on the stars
    starother = starB
    rkStar(starA)
    starother = starA
    rkStar(starB)
    
    
    # run Runge-Kutta on the particles
    rkParticles()
    if len(particle_list) < 2000:#2000: # this is not ultimately how we will do this
        AddParticle()
    
    starA.trail.append(pos=starA.pos)
    starB.trail.append(pos=starB.pos)

    #while not done:
     #   particleCol()
      #  done = 1
    
    particleCol()
    
    rate(500) 
    