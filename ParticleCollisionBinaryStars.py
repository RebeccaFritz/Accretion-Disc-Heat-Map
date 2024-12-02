Web VPython 3.2
from random import randrange

# create canvases
scene.autoscale = 1
scene = canvas(title = 'Particle Collisions in a Binary Star System', color = color.black, align='left')

# Define Globals
G  = 6.67E-11 # gravitational constant
AU = 1.5E11 # astronomical unit
SR = 0.0046524726*AU # solar radius
YEAR = 365.25*24*60*60 # year in seconds
MS = 1.988400e20 # mass of the sun 
MID_COLOR = vec(0.98, 0.59, 0.20) # color of particles when ambient temp : rgb(250,150, 0.20)
HIGH_COLOR = vec(0.78, 0.00, 0.00) # color of particles when hot: rgb(200, 000, 000)
LOW_COLOR = vec(0.98, 0.78, 0.00) # color of particles when cold: rgb(250, 200, 000)
STARA_COLOR = vec(0.76, 0.81, 0.95) #rgb(193,207,242)
STARB_COLOR = vec(0.70, 0.04, 0.20) #rgb(179,9,52)
A_RADIUS = 6.6E6*15 # accretion disk radius
h = 2.0*1e6 # time step                            

# convert hydrogen atom to particle 
u = 1.66e-27 # atomic mass in kg
H_MASS = 1.01*u # mass of a hydrogen atom 
H_RADIUS = 120e-10 # van der Waals radius of a hydrogen atom 
Hscaler = 1e14*4 # one particle is 1e14*4 atoms (200 000 000 000 000)
P_RADIUS = H_RADIUS*Hscaler # particle radius
P_MASS = H_MASS*Hscaler # particle mass

# create scene objects
starA = sphere(color=STARA_COLOR, visible=1) # white dwarf (primary star)
starB = sphere(color=STARB_COLOR, opacity=0.5, visible=1) # secondary star 

# EPIC 220615486 star system details from Montgomery et. al. paper
# Note: these are from section 3.4 "... Computational Modeling" unless otherwise stated
P = 0.065837*24*60*60 # orbital period in seconds
#r0 = 0.590*SR # from section 3.3 "... Analytical Models"  # distance between starA and starB
r0 = 0.5
starA.mass = 0.062*MS # assumed
starA.radius = 0.0081*SR*r0 # scaled down by a
starB.mass = 0.040*MS # estimate
starB.radius = 0.64*SR*r0 # scaled down by a # from section 3.3 #... Analytical Models" 
starA.pos = vec(0,0,0)
starB.pos = vec(r0*SR,0,0) 


# legend = # Eventually 
starA.trail = curve(pos=starA.pos, color=starA.color) # do we want trails
starB.trail = curve(pos=starB.pos, color=starB.color)
# accretionDisk = # do we still want to do this? could do two transparent cylinders compounded together

# set initial velocities
starB.vel = vector(0,sqrt(G*(starA.mass)/(r0*SR)),0)*0.8 # vis-viva equation     #0.922 if a is semi-major       
starA.vel = -starB.vel*(starB.mass/starA.mass)

# ancillaries
t = 0 # time
particle_list = []

# a function to determine the acceleration of a star
## vector --> vector
def accStar(starpos):
    return -G*starother.mass*(starpos-starother.pos)/r**3

# a function to run the Runge-Kutta alg on a star
## object --> void
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

# a function to 'remove' particles from the particle_list if it is out of range 
#   OR to 'remove' them and to increase starB's mass if they run into starB
#   NOTE: this really just makes the particle invisible
## ... --> ...
def removeParticleOriginal(): # works with sample system
    for particle in particle_list:
        # to check if the particle is outside of an imaginary circle with:
        #   center: centered halfway between the starting position of the donor at <5.80E8/2, 0, 0>
        outerCircle = (particle.pos.x - 0)**2 + (particle.pos.y - 0)**2
        # to check if the particle has joined the accreting star
        particleCircle = (particle.pos.x - starA.pos.x)**2 + (particle.pos.y - starA.pos.y)**2 
        
        if outerCircle > 5.80E8**2: # if the particle leaves the system
            particle.visible = False
            #particle.color = color.blue # for testing
        else if particleCircle < starA.radius**2: # if the particle joins the accreting star
            particle.visible = False
            starA.mass += particle.mass
#        else:     # for testing
#            particle.color = MID_COLOR
            
# a function to 'remove' particles from the particle_list if it is out of range 
#   OR to 'remove' them and to increase starB's mass if they run into starB
#   NOTE: this really just makes the particle invisible
## ... --> ...
def removeParticle():
    for particle in particle_list:
        # to check if the particle is outside of an imaginary circle around the system
        outerCircle = (particle.pos.x - 0)**2 + (particle.pos.y - 0)**2
        # to check if the particle has joined the accreting star
        innerCircle = (particle.pos.x - starA.pos.x)**2 + (particle.pos.y - starA.pos.y)**2 
        
        if outerCircle > (r0*SR)*2: # if the particle leaves the system
            #particle.visible = False
            particle.color = color.blue # for testing
        else if innerCircle < starA.radius**2: # if the particle joins the accreting star
            particle.visible = False
            starA.mass += particle.mass
#        else:     # for testing
#            particle.color = MID_COLOR

# a function to add particles to the particle_list and to decrease starA's mass
## void --> void
def AddParticle():
    # This may not be how we want to create particles
    Omega = 2.0*pi/P # angular velocity                                                   # from sample code
    for i in range(0, 2):
        # particle mass and raidus are split into the number being released and given a random additional hydrogen atoms within the range 0 to 10
        pScaler = random()*10
        particle_mass = P_MASS + pScaler*H_MASS 
        particle_radius = P_RADIUS + pScaler*H_RADIUS
        
        particle = sphere(mass = particle_mass, radius = particle_radius, color = MID_COLOR)
        #particle.pos = starA.pos + 0.6367*(starB.pos-starA.pos) + vec(i*particle_radius*3, 0, 0)  # adjusted from sample code
        particle.pos = starA.pos + 0.45*(starB.pos-starA.pos) + vec(i*particle_radius*3, 0, 0)
        # set initial velocities
        particle.vel = vector(0,0,0)                                           # from sample code
        particle.vel.x = -particle.pos.x*Omega*(1.0 + random())                # from sample code # v = radius*angular velocity
        particle.vel.y = particle.pos.y*Omega*(1.0 + random())                 # from sample code
        particle.vel.z = particle.pos.z*Omega*(1.0 + random())                 # from sample code
        ##
        particle_list.append(particle)
        starA.mass = starA.mass - particle_mass # remove the mass of the particle from the star

# print particle velocites before and after a collison (for testing purposes)
# vec, vec, vec, vec --> print()
counter = 200 # for printing purposes
def printCollVels(v1i, v2i, v1f, v2f):
    global counter
    if counter >= 200:
        print('particle : vel-init <x,y,z>, vel-final <x,y,z> ')
        print('one:   {0},   {1}'.format(v1i, v1f))
        print('two:   {0},   {1}'.format(v2i, v2f))
        print()
        counter = 0
    counter += 1

#Final velocities of particles in collisions
# int, int --> void
def velfinal(particleidx, otheridx):
    # get masses, positions, inital velocites
    m1 = particle_list[particleidx].mass
    m2 = particle_list[otheridx].mass
    r1 = particle_list[particleidx].pos
    r2 = particle_list[otheridx].pos
    v1i = particle_list[particleidx].vel
    v2i = particle_list[otheridx].vel
    ## 
    
    v1f = v1i - (2*m2)/(m1+m2)*( (dot((v1i-v2i), (r1-r2))) / mag(r1-r2)**2 ) * (r1-r2) 
    particle_list[particleidx].vel = v1f
    
    v2f = v2i - (2*m1)/(m1+m2)*( (dot((v2i-v1i), (r2-r1))) / mag(r2-r1)**2 ) * (r2-r1)
    particle_list[otheridx].vel = v2f
    
    # print velocities (for testing purposes)
    #printCollVels(v1i, v2i, v1f, v2f)
        
    return        

# a function to check if the particles have collided
## ... --> ...
def particleColBounded():
    particleDiameter = P_RADIUS*2
    collisions = [] # list of indicies of the particles for which the collisions have been calculated
    
    # find every colliding pair and add it to collisions
    for i in range(0, len(particle_list)):
        particle1 = particle_list[i]

        particle1Circle = (particle1.pos.x - starA.pos.x)**2 + (particle1.pos.y - starA.pos.y)**2 
        if particle1.visible and particle1Circle <= A_RADIUS**2: #mag(L1point.pos)**2:  # is particle1 in the accreation disk
            #particle1.color = color.red # for testing
            for j in range(i+1, len(particle_list)):
                particle2 = particle_list[j]
                
                particle2Circle = (particle2.pos.x - starA.pos.x)**2 + (particle2.pos.y - starA.pos.y)**2 
                if particle2.visible and particle2Circle <= A_RADIUS**2: #mag(L1point.pos)**2: # is particle2 in the accreation disk
                    if mag(particle2.pos - particle1.pos) <= particleDiameter: 
                        collisions.append([i,j]) # the index of colliding pair [particle1, particle2]
        #else: particle1.color = color.blue # for testing
        
    # collisions now has the indicies of all the particles colliding with particle1
    if len(collisions) > 1:
        # calculate velocity chage due to elastic particle collision for all particles at the indicies in the collisions list
        for ij in collisions:
            velfinal(ij[0], ij[1])
            
# a function to check if the particles have collided
## ... --> ...
def particleCol():
    particleDiameter = P_RADIUS*2
    collisions = [] # list of indicies of the particles for which the collisions have been calculated

    # find every colliding pair and add it to collisions
    for i in range(0, len(particle_list)):
            particle1 = particle_list[i]
            
            #particle1Circle = (particle1.pos.x - starA.pos.x)**2 + (particle1.pos.y - starA.pos.y)**2 
            if particle1.visible: #and particle1Circle <= A_RADIUS**2: # is particle1 in the accreation disk
                for j in range(i+1, len(particle_list)):
                    particle2 = particle_list[j]

                    #particle2Circle = (particle2.pos.x - starA.pos.x)**2 + (particle2.pos.y - starA.pos.y)**2 
                    if particle2.visible: # and particle2Circle <= A_RADIUS**2: # is particle2 in the accreation disk
                        if mag(particle2.pos - particle1.pos) <= particleDiameter: 
                            collisions.append([i,j]) # the index of colliding pair [particle1, particle2]
#                            particle1.color = color.blue # for testing
#                            particle2.color = color.blue

    # collisions now has the indicies of all the particles colliding with particle1
    if len(collisions) > 1:
        # calculate velocity chage due to elastic particle collision for all particles at the indicies in the collisions list
        for ij in collisions:
            velfinal(ij[0], ij[1])


## a function to run the Runge-Kutta alg on particles       # from the sample code
### void --> void
#def rkParticles():
#    for particle in particle_list:  
#        rA = mag(particle.pos - starA.pos)
#        rB = mag(particle.pos - starB.pos)
#        if rA > L1/100 and rB > L2/100 and particle.visible: # I'm pretty sure this statement is only false once the particle's shoot out
#            acc_p = G*starA.mass*(starA.pos - particle.pos)/rA**3 + G*starB.mass*(starB.pos - particle.pos)/rB**3
#            k1v = h*acc_p
#            k1x = h*particle.vel
#            k2v = h*(G*starA.mass*(starA.pos - (particle.pos+k1x/2.0))/rA**3 + G*starB.mass*(starB.pos - (particle.pos+k1x/2.0))/rB**3)
#            k2x = h*(particle.vel + k1v/2.0)
#            k3v = h*(G*starA.mass*(starA.pos - (particle.pos+k2x/2.0))/rA**3 + G*starB.mass*(starB.pos - (particle.pos+k2x/2.0))/rB**3)
#            k3x = h*(particle.vel + k2v/2.0)
#            k4v = h*(G*starA.mass*(starA.pos - (particle.pos+k3x/2.0))/rA**3 + G*starB.mass*(starB.pos - (particle.pos+k3x/2.0))/rB**3)            
#            k4x = h*(particle.vel + k3v)
#            particle.vel += (k1v + 2.0*k2v + 2.0*k3v + k4v)/6.0
#            particle.pos += (k1x + 2.0*k2x + 2.0*k3x + k4x)/6.0
            
# a function to determine the acceleration of a particle
## vector --> vector
def accParticle(particlepos):
    return G*starA.mass*(starA.pos - particlepos)/mag(particlepos - starA.pos)**3 + G*starB.mass*(starB.pos - particlepos)/mag(particlepos - starB.pos)**3

# a function to run the Runge-Kutta alg on particles       
## void --> void
def rkParticles():
    for particle in particle_list: 
        if particle.visible:
            rk1v = h*accParticle(particle.pos)
            rk1x = h*particle.vel
            
            rk2v = h*accParticle(particle.pos + rk1x/2.0)
            rk2x = h*(particle.vel + rk1v/2.0)
            
            rk3v = h*accParticle(particle.pos + rk2x/2.0)
            rk3x = h*(particle.vel + rk2v/2.0)
            
            rk4v = h*accParticle(particle.pos + rk3x)  
            rk4x = h*(particle.vel + rk3v)
            
            particle.vel += (rk1v + 2.0*rk2v + 2.0*rk3v + rk4v)/6.0
            particle.pos += (rk1x + 2.0*rk2x + 2.0*rk3x + rk4x)/6.0
        
# Eventually:
# a function to check particle temperature and change their color
## ... --> ...
#def particleTemp(...):
 #   return

        
# run animation
counter2 = 30
while True:
    r = mag(starA.pos - starB.pos)
    
    # run Runge-Kutta on the stars
    starother = starB
    rkStar(starA)
    starother = starA
    rkStar(starB)
    
    starA.trail.append(pos=starA.pos)
    starB.trail.append(pos=starB.pos)
    


#    # run Runge-Kutta on the particles
#    rkParticles()
    if len(particle_list) < 2000:: 
        AddParticle()
        
#    # start calculating collisions
    #particleCol()
#    particleColBounded()
#    
#    removeParticle() # does not work
    removeParticleOriginal() # from sample system
        
    rate(50) #500
