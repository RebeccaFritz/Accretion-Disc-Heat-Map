Web VPython 3.2
#from vpython import*
from random import randrange

# create canvases
scene = canvas(title = 'Particle Collisions in a Binary Star System', color = color.black, align='left')

# Define Globals
G  = 6.67E-11 # gravitational constant
AU = 1.5E11 # astronomical unit in meters
SR = 6.95700E8 # solar radius in meters
YEAR = 365.25*24*60*60 # year in seconds
MS = 1.988400E30 # mass of the sun in kg
MID_COLOR = vec(0.98, 0.59, 0.20) # color of particles when ambient temp : rgb(250,150, 0.20)
HIGH_COLOR = vec(0.78, 0.00, 0.00) # color of particles when hot: rgb(200, 000, 000)
LOW_COLOR = vec(0.98, 0.78, 0.00) # color of particles when cold: rgb(250, 200, 000)
STARA_COLOR = vec(0.76, 0.81, 0.95) #rgb(193,207,242)
STARB_COLOR = vec(0.70, 0.04, 0.20) #rgb(179,9,52)
k = 1.380649E-23    #Boltzmann constant in Joules per Kelvin

# convert hydrogen atom to particle 
u = 1.66e-27 # atomic mass in kg
H_MASS = 1.01*u # mass of a hydrogen atom 
H_RADIUS = 120e-10 # van der Waals radius of a hydrogen atom in m
Hscaler = 1e14*4 #1e14*4 # one particle is 1e14*4 atoms (400 000 000 000 000)
P_RADIUS = H_RADIUS*Hscaler # particle radius
P_MASS = H_MASS*Hscaler # particle mass
   
########## EPIC 220615486 star system 
# create scene objects
starA = sphere(color=STARA_COLOR, visible=1) # white dwarf (primary star)
starB = sphere(color=STARB_COLOR, visible=1) # secondary star
RocheLobeB = sphere(color=STARB_COLOR, opacity=0.5, visible=1)
AccDisc = cylinder(pos=vec(0,0,0), axis=vec(0,0,0.5), color=vec(1.00, 0.88, 0.49), opacity=0.35)

P_days = 0.065837 # orbital period in days                                                          # from Montgomery et. al. paper
P = P_days*24*60*60 # orbital period in seconds
starA.mass = 0.62*MS # assumed                                                                      # from Montgomery et. al. paper
starB.mass = 0.040*MS # estimate                                                                    # from Montgomery et. al. paper

### calculating the distance between the stars using equation 18 in Montgomery et. al. paper
# which is Newton's version of Kepler's Third Law
# - period is in years
# - distance is in AU
# - G is in terms of solar masses not kg: 39.46 AU**3 / (yr**2 * MS)
a = ((P_days/(365.25))**2 * (39.46/MS) * (starA.mass+starB.mass) / (4*pi**2) )**(1/3) # in Astronomical units 
a = a*214.93946938 # convert to number of SR
### 

starA.radius = 0.0081*SR
starB.radius = (starB.mass/MS)**(13/15)*SR                                                          # equation 19 from Montgomery et. al.
starA.pos = vec(0,0,0)
starB.pos = vec(a*SR,0,0)                                                                           # from Montgomery et. al. paper
h = P/2000 # time step 
q = starB.mass/starA.mass #ratio of star masses                                                     # from Montgomery et. al. paper
L1 = a*(0.500 - 0.227*log10(q)) * SR # approximate distance of L1 from starA
L2 = L1
A_RADIUS = 0.6*a/(1 + q)*SR  #average accretion disc radius                                         # from Montgomery et. al. paper
RocheLobeB.radius = a*(0.500 + 0.227*log10(q)) * SR # approximate distance of L1 from starB         # from Montgomery et. al. paper
AccDisc.radius = A_RADIUS

# set initial velocities
starB.vel = vector(0,sqrt(G*(starA.mass)/(a*SR)),0)   #m/s     
starA.vel = -starB.vel*(starB.mass/starA.mass)

scene.range = a*SR*1.35
##########

# adjusting positions relative to the center of mass
com = (starA.mass*starA.pos+starB.mass*starB.pos)/(starA.mass+starB.mass) # center of mass <2.51741e7,0,0>
starA.pos = starA.pos-com
starB.pos = starB.pos-com
RocheLobeB.pos = starB.pos
AccDisc.pos = starA.pos

starA.trail = curve(pos=starA.pos, color=starA.color, retain=2500) # do we want trails
starB.trail = curve(pos=starB.pos, color=starB.color, retain=2500)
L1point = sphere(pos = starA.pos + L1*vec(1,0,0), radius = P_RADIUS, color = color.blue, visible=0)

# ancillaries
t = 0 # time
particle_list = []
invisible_list = [] # contains indicies of the particles that are not visible

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
    
# a function to update the position of a lagrange point
## void --> void
def posLagrange(): 
    approx_dist = mag(starB.pos-starA.pos)*(0.500 - 0.227*log10(starB.mass/starA.mass)) 
    
    ang = diff_angle(vec(starA.pos.x, 0,0), starB.pos-vec(0, starA.pos.y, 0)) # angle between two vectors
    
    if starB.pos.x >= 0 and starB.pos.y >= 0: # quadrant 1
        L1point.pos = starA.pos - vec((approx_dist)*cos(-ang), (approx_dist)*sin(-ang), 0)
    elif starB.pos.x < 0 and starB.pos.y < 0: # quadrant 3
        L1point.pos = starA.pos + vec((approx_dist)*cos(-ang), (approx_dist)*sin(-ang), 0)
    elif starB.pos.x >= 0 and starB.pos.y < 0: # quadrant 4
        L1point.pos = starA.pos - vec((approx_dist)*cos(ang), (approx_dist)*sin(ang), 0) 
    elif starB.pos.x < 0 and starB.pos.y >= 0: # quadrant 2
        L1point.pos = starA.pos + vec((approx_dist)*cos(ang), (approx_dist)*sin(ang), 0)
            
# a function to 'remove' particles from the particle_list if it is out of range 
#   OR to 'remove' them and to increase starB's mass if they run into starB
#   NOTE: this really just makes the particle invisible
## ... --> ...
def removeParticle():
    for particle in particle_list:
        # to check if the particle is outside of an imaginary circle around the system
        outerCircle = (particle.pos.x - 0)**2 + (particle.pos.y - 0)**2 + (particle.pos.z - 0)**2
        # to check if the particle has joined the accreting star
        innerCircle = (particle.pos.x - starA.pos.x)**2 + (particle.pos.y - starA.pos.y)**2 + (particle.pos.z - starA.pos.z)**2
        # to check if the particle has rejoined starB
        starBCircle = (particle.pos.x - starB.pos.x)**2 + (particle.pos.y - starB.pos.y)**2 + (particle.pos.z - starB.pos.z)**2
        
        if outerCircle >= (1.5*a*SR)**2: # if the particle leaves the system #mag(L1point.pos)**2
            particle.visible = False
            invisible_list.append(particle.index)
        elif innerCircle <= (starA.radius)**2: # if the particle joins the dwarf star 
            particle.visible = False
            starA.mass += particle.mass
            invisible_list.append(particle.index)
        elif starBCircle <= (mag(starB.pos - L1point.pos)-P_RADIUS)**2:  # if the particle rejoins the donor star
            particle.visible = False
            invisible_list.append(particle.index)
            starB.mass += particle.mass

# a function to add position, visibility, index, and velocity to a particle and to remove its mass from the starB
# object, int --> void
def alterParticle(particle, i):
    sim_P = 0.050427*24*60*60 # period of the simulation in seconds (slightly lower than then period of the actual system)
    Omega = 2.0*pi/sim_P # angular velocity #mag(starB.vel)/mag(starB.pos)
    # particle mass and raidus are split into the number being released and given a random additional hydrogen atoms within the range 0 to 10
#   pScaler = random()*10
#   particle_mass = P_MASS + pScaler*H_MASS 
#   particle_radius = P_RADIUS + pScaler*H_RADIUS
    
    particle_mass = P_MASS # equal as in Montgomery et. al. paper
    particle_radius = P_RADIUS # equal as in Montgomery et. al. paper
    
    particle.mass = particle_mass
    particle.radius = particle_radius 
    particle.visible = 1
    particle.temp = 0 # kelvin
    
    # releasing particles at L1 point 
    if L1point.pos.x >= 0:
        particle.pos = L1point.pos - vec(i*particle_radius*2.5, 0, 0) # is so that if more than one particle is released at a time they are not released on top of each other
    else:
        particle.pos = L1point.pos + vec(i*particle_radius*2.5, 0, 0) 
    
    # set initial velocities
    particle.vel = vector(0,0,0)                                           # from Dr. Olenick's sample code
    particle.vel.x = -particle.pos.x*Omega*(1.0 + random())                # from Dr. Olenick's sample code
    particle.vel.y = -particle.pos.y*Omega*(1.0 + random())                             
    particle.vel.z = particle.pos.z*Omega*(1.0 + random())                 # from Dr. Olenick's sample code 
    
    starB.mass = starB.mass - particle_mass # remove the mass of the particle from the donor star
    
# a function to add particles to the particle_list
## void --> void
def AddParticle():                                              # from sample code
    for i in range(0, 1):
        
        particle = sphere(color = MID_COLOR)
        alterParticle(particle, i)
        particle_list.append(particle)
        particle.index = len(particle_list)-1

# a function to rerelease a removed particle
## void --> void
def rereleaseParticles():
    
    stop = 1
    if len(invisible_list) < stop:
        stop = len(invisible_list)
        
    for i in range(0, stop):
        index = invisible_list[i]
        particle = particle_list[index]
        if not particle.visible:
            
            alterParticle(particle, i)
            particle.color = MID_COLOR
            
    del invisible_list[0:stop]
            
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
    
    v1f = v1i - ((2*m2)/(m1+m2)) * ( (dot((v1i-v2i), (r1-r2))) / mag(r1-r2)**2 ) * (r1-r2) 
    particle_list[particleidx].vel = v1f
    
    v2f = v2i - ((2*m1)/(m1+m2)) *( (dot((v2i-v1i), (r2-r1))) / mag(r2-r1)**2 ) * (r2-r1)
    particle_list[otheridx].vel = v2f
 

# a function to check if the particles have collided
## void --> void
def particleColBounded():
    particleDiameter = P_RADIUS*2
    collisions = [] # list of indicies of the particles for which the collisions have been calculated
    
    # find every colliding pair and add it to collisions
    for i in range(0, len(particle_list)):
        particle1 = particle_list[i]

        particle1Circle = (particle1.pos.x - starA.pos.x)**2 + (particle1.pos.y - starA.pos.y)**2 + (particle1.pos.z - starA.pos.z)**2 
        if particle1.visible and particle1Circle <= A_RADIUS**2: # is particle1 in the accreation disk
            for j in range(i+1, len(particle_list)):
                particle2 = particle_list[j]
                
                particle2Circle = (particle2.pos.x - starA.pos.x)**2 + (particle2.pos.y - starA.pos.y)**2 + (particle2.pos.z - starA.pos.z)**2 
                if particle2.visible and particle2Circle <= A_RADIUS**2: # is particle2 in the accreation disk
                    if mag(particle2.pos - particle1.pos) <= particleDiameter: 
                        collisions.append([i,j]) # the index of colliding pair [particle1, particle2]
        
    # collisions now has the indicies of all the particles colliding with particle1
    if len(collisions) > 1:
        # calculate velocity chage due to elastic particle collision for all particles at the indicies in the collisions list
        for ij in collisions:
            velfinal(ij[0], ij[1])
          
# a function to run the Runge-Kutta alg on particles (with gravitational softening)    
## void --> void
def rkParticles():
    for particle in particle_list:  
        alpha = 0.5e-10
        beta = 0.5e-10    
        rA = mag(particle.pos - starA.pos)
        rB = mag(particle.pos - starB.pos)
        hp = h # particle time step
        if rA > L1/100 and rB > L2/100 and particle.visible: 
            #### implement grav softening
            particle_circ = (particle.pos.x - starA.pos.x)**2 + (particle.pos.y - starA.pos.y)**2 + (particle.pos.z - starA.pos.z)**2 # the radius of a circle made with the particle (x,y) and starA (center)
            epsilon = 1000*SR
#            if particle_circ >= (A_RADIUS)**2:
#                epsilon = 0
            
            acc_p = G*starA.mass*(starA.pos - particle.pos)/(rA**3 + epsilon**2) + G*starB.mass*(starB.pos - particle.pos)/(rB**3 + epsilon**2)
            k1v = hp*(acc_p)
            k1x = hp*particle.vel
            
            acc_2 = (G*starA.mass*(starA.pos - (particle.pos+k1x/2.0))/(rA**3 + epsilon**2) + G*starB.mass*(starB.pos - (particle.pos+k1x/2.0))/(rB**3 + epsilon**2))
            k2v = hp*(acc_2)
            k2x = hp*(particle.vel + k1v/2.0)
            
            acc_3 = (G*starA.mass*(starA.pos - (particle.pos+k2x/2.0))/(rA**3 + epsilon**2) + G*starB.mass*(starB.pos - (particle.pos+k2x/2.0))/(rB**3 + epsilon**2))
            k3v = hp*(acc_3)
            k3x = hp*(particle.vel + k2v/2.0)
            
            acc_4 = (G*starA.mass*(starA.pos - (particle.pos+k3x/2.0))/(rA**3 + epsilon**2) + G*starB.mass*(starB.pos - (particle.pos+k3x/2.0))/(rB**3 + epsilon**2))
            k4v = hp*(acc_4)            
            k4x = hp*(particle.vel + k3v)
            
            particle.vel += (k1v + 2.0*k2v + 2.0*k3v + k4v)/6.0
            particle.pos += (k1x + 2.0*k2x + 2.0*k3x + k4x)/6.0    
            
            # updating the heatmap and the particle temperature
            particleTemp(particle)


# a function to check particle temperature and change their color
## object --> void
def particleTemp(particle):
    Temp = particle.mass*mag2(particle.vel)/3*k # particle's temperature in kelvin
    
    hotTemp = particle.mass*(2766183.366**2)/3*k ## approx max vel: 2766183.366 
    warmTemp = particle.mass*(2330811.3152**2)/3*k
    chillTemp = particle.mass*(1895439.2644**2)/3*k
    subChillTemp = particle.mass*(1460067.2136**2)/3*k
    lessColdTemp = particle.mass*(1024695.1628**2)/3*k
    coldTemp = particle.mass*(589323.112**2)/3*k ## approx min vel: 589323.112 
    
    if Temp >= hotTemp:
        particle.color = vec(0.88, 1.00, 1.00) # white blue   #hottest Temp 
    elif Temp >= warmTemp and Temp < hotTemp:
        particle.color = vec(1.00, 0.98, 0.80) # white yellow 
    elif Temp >= chillTemp and Temp < warmTemp:
        particle.color = vec(0.99, 0.97, 0.37) # light yellow
    elif Temp >= subChillTemp and Temp < chillTemp:
        particle.color = vec(1.00, 0.73, 0.00) # yellow
    elif Temp >= lessColdTemp and Temp < subChillTemp:
        particle.color = vec(1.00, 0.46, 0.09) # orange
    elif Temp >= coldTemp and Temp < lessColdTemp:
        particle.color = vec(1.00, 0.22, 0.00) # red
    else: # Temp < coldTemp
        particle.color = vec(0.60, 0.00, 0.00) # dark red
        
    particle.temp = Temp
    

# more ancillaries
rmin = 0
rmax = SR

releaseCounter = 10

vmin = 9e90
vmax = 0

#define position vs. time plots
plot1 = graph(title='Particle Distance as a function of Time', xtitle='Time [s]', ytitle='Distance [SR]', width=800, height=450, align = 'right') 
p1curvePos = gcurve(color=color.red, label = 'particle 0', legend = True, width=2)
p2curvePos = gcurve(color=color.cyan, label = 'particle 74', legend = True, width=2)
p3curvePos = gcurve(color=color.green, label = 'particle 149', legend = True, width=2)


#define temp vs. time plots
plot2= graph(title='\nChange in Particle Temperature as a function of Time', xtitle='Time [s]', ytitle='Change in Temp [K]', width=820, height=450, align = 'right') 
p1curveTemp = gcurve(color=color.red, label = 'particle 0', legend = True, width=2)
p2curveTemp = gcurve(color=color.cyan, label = 'particle 74', legend = True, width=2)
p3curveTemp = gcurve(color=color.green, label = 'particle 149', legend = True, width=2)

#define velocity vs. time plots
plot2= graph(title='\nParticle Velocity as a function of Time', xtitle='Time [s]', ytitle='Velocity [SR/s]', width=815, height=450, align = 'right') 
p1curveVel = gcurve(color=color.red, label = 'particle 0', legend = True, width=2)
p2curveVel = gcurve(color=color.cyan, label = 'particle 74', legend = True, width=2)
p3curveVel = gcurve(color=color.green, label = 'particle 149', legend = True, width=2)


# run animation
while True:
    r = mag(starA.pos - starB.pos)
    
    # run Runge-Kutta on the stars
    starother = starB
    rkStar(starA)
    starother = starA
    rkStar(starB)
    
    # update star attachment positions
    RocheLobeB.pos = starB.pos
    AccDisc.pos = starA.pos        
    starA.trail.append(pos=starA.pos)
    starB.trail.append(pos=starB.pos)
    
    # position lagrange 1 and 2
    posLagrange() 

    # run Runge-Kutta on the particles
    rkParticles()

    # release or rerelease one particle every 10 time steps, prioritizing rereleasing particles 0, 74, and 149 which are needed for graphing
    if releaseCounter == 10:
        if 0 in invisible_list:
            alterParticle(particle_list[0], 0) # release particle 0 if invisible for graphing purposes
            invisible_list.remove(0)
        elif 74 in invisible_list:
            alterParticle(particle_list[74], 0) # release particle 74 if invisible for graphing purposes
            invisible_list.remove(74)
        elif 149 in invisible_list:
            alterParticle(particle_list[149], 0) # release particle 149 if invisible for graphing purposes
            invisible_list.remove(149)
        elif len(particle_list) < 2000: 
            AddParticle()
        elif len(invisible_list) > 0: 
            rereleaseParticles()  
        releaseCounter = 0
    releaseCounter += 1
    
    # start calculating collisions
    particleColBounded()
    
    # remove particles
    removeParticle()
    
    if t>5000 and t<15000: 
        if particle_list[0].visible: 
            p1curvePos.plot(t/50, mag(starA.pos-particle_list[0].pos)/SR)
            p1curveTemp.plot(t/50, particle_list[0].temp)
            p1curveVel.plot(t/50, mag(particle_list[0].vel/SR))
        if len(particle_list) >= 75 and particle_list[74].visible:
            p2curvePos.plot(t/50, mag(starA.pos-particle_list[74].pos)/SR)
            p2curveTemp.plot(t/50, particle_list[74].temp)
            p2curveVel.plot(t/50, mag(particle_list[74].vel/SR))
        if len(particle_list) >= 150 and particle_list[149].visible:
            p3curvePos.plot(t/50, mag(starA.pos-particle_list[149].pos)/SR)
            p3curveTemp.plot(t/50, particle_list[149].temp)
            p3curveVel.plot(t/50, mag(particle_list[149].vel/SR))
               
    t += h
    
    rate(500) #500
