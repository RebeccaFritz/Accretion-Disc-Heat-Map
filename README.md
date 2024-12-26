# Accretion Disc Heat-Map Simulation
We created this project as the final project for the Computational Physics class at the University of Dallas in Fall 2024. 

## Creators: 
**Rebecca Fritz**
<br>
Department of Computer Science, University of Dallas, 1845 E. Northgate Drive, Irving, 75062
<br>
E-mail: rgfritz@udallas.edu
<br>
**Elizabeth Jacoby**
<br>
Department of Physics, University of Dallas, 1845 E. Northgate Drive, Irving, 75062
<br>
E-mail: ejacoby@udallas.edu
<br>
<br>
*December 14, 2024*
<br>
<br>
<br>

## ABSTRACT
This simulation models the cataclysmic variable system EPIC 220615486. Using VPython, a binary star system composed of an accreting white dwarf and companion, mass donor star, was simulated. The matter lost by the donor-star was modeled as individual particles rather than as a fluid. The gravity of the white dwarf pulls these particles through Lagrange Point 1 (L1) and into its orbit: forming an accretion disc. Once the particles pass through L1, the interactions of these particles are tracked. Assuming elastic collisions, a fourth-order Runge-Kutta algorithm was implemented to measure the velocity changes from collisions between the particles in the accretion disc. To obtain a stable accretion disc with collisions, the change of velocity (in the event of a collision) was scaled to ten percent. A heat-map of the accretion disc was then constructed with the information gained about the velocities from this simulation.
<br>
<br>
**Key words:** Cataclysmic Variable -- Accretion Disc -- Temperature-Map -- Simulation
<br> <br> <br>

## 1. INTRODUCTION
Binary star systems evolve in many different ways depending on multiple variables. For example, the stars' mass, proximity to each other, and age are a few such variables. Our simulation investigates the mechanics of a Cataclysmic Variable binary star system (CV). Our CV—EPIC 220615486—is composed of a white dwarf and a red giant, companion star, orbiting each other in a close binary. The dwarf star is denser and its gravitational potential is significantly stronger than the companion star. Further, the companion star is of lower mass and, due to its age, has expanded and filled its Roche lobe. A Roche lobe is the three-dimensional area around a star in a binary in which objects inside are gravitationally bound to that star. As the stars orbit each other, the strength of the gravity of the dwarf extrudes matter from the companion star and pulls this matter towards itself. The matter enters into the Roche lobe of the dwarf and orbits the dwarf in a disc (visually similar to the rings of Saturn) until the matter spirals in to join the dwarf star. This process is known as accretion and the disc is known as an accretion disc (Cataclysmic Variables).
<br>
<br>
The conditions in an accretion disc vary throughout the disc. In our simulation, we primarily focus on investigating the heat changes in differing areas of the disc. The disc is modeled by individual particles being released from the companion star and drawn into the orbit of the dwarf. These particles are not only affected by the gravitational forces of the dwarf and giant, but they are also colliding with each other. By measuring the velocities of the particles from their interactions, the temperature of each particle is calculated.
<br>
<br>
This simulation aims to model a CV and its accretion disc. The matter transfer will be simulated with point particles and the interaction of these particles will be tracked and used to determine a map of the heats in different areas of the disc. 
<br> <br> <br>

## 2. PROCEDURES
To start, we adapted code provided by Dr. Richard Olenick that simulated a accretion disk in a unnamed system to the system parameters provided in a paper he collaborated on with Montgomery et al. in 2024. After struggling to get a stable orbit for the EPIC 220615486 dwarf and donor stars with the values provided in the Montgomery paper, we made the observation that the listed radius of the secondary star was larger than the average distance separating the centers of the two stars. In order to verify these numbers, we recalculated the average distance between the stars and the radius of the secondary using data and equations from that same paper. The data used from the paper can be found in table 1 (Montgomery et al. 2024). We used Newton’s version of Kepler’s Third Law as expressed in equation 1 for the distance between the stars. This was also used by Montgomery's team in their paper. In this equation, the masses are in solar masses, the period is in years, and the distance—a— comes out in AU, it was converted to solar radii to match the measurements we already had. We found that our number, 0.597 $R_{\odot}$, was very close to Montgomery's value, 0.590 $R_{\odot}$. It is expected that the values are not exactly the same since we used rounded values found in the computation section of the paper rather than the raw values from the analytic section. 
<br><br>
$$Eq. 1$$
<br>
$$P^{2}_{orb}=\frac{4\pi^{2}a^3}{G(M_1+M_2)}$$
<br><br> 
![Data Table](https://github.com/user-attachments/assets/d921d83a-c672-4521-9eef-633b43775b54)
<br>
For the radius of the secondary star we used equation 2 which relates the radius of the star to the ratio of that star’s mass to the mass of the sun. The number we obtained was 0.064 $R_{\odot}$ which has a decimal point difference from the number provided in the paper, which is 0.64 $R_{\odot}$. Along with the data we had before, these values allowed us to calculate a stable orbit for the two stars.
<br><br>
$$Eq. 2$$
<br>
$$R_{2}=\left(\frac{M_2}{M\odot}\right)^{\frac{13}{15}}$$
<br><br>
In order to calculate a stable orbit we used Newton’s second law solved for orbital speed (Eq. 3) to get the velocity of the secondary star. We then multiplied that by the ratio of the masses and took the inverse to get the velocity of the dwarf star (Eq. 4). We then used a provided code snippet from Dr. Olenick in order to calculate the movement of the stars for each time step. This places a fourth-order Runge-Kutta algorithm inside a VPython function which takes a star object. 
<br><br>
$$Eq. 3$$
<br>
$$v_{orb} = \sqrt{\frac{GM}{r}}$$
<br><br>
$$Eq. 4$$
<br>
$$v_A = -v_B \frac{M_B}{M_A}$$
<br><br>
In the first step of the Runge-Kutta, a velocity is determined by taking the acceleration of that star at the start position and multiplying that by the time step, which in this program is 0.02 seconds. A position is then determined by multiplying the time step by the velocity of the star. For the second, third, and fourth steps this is repeated, but the position and velocity are modified by the previously calculated position or velocity. Finally, the star object’s velocity and position are each updated using the calculated values. See figure 1.
<br>
<br>
**Fig 1: The fourth order Runge-Kutta**
<br>
![Runge-Kutta](https://github.com/user-attachments/assets/bfddbae7-4cab-468a-802d-641ab1ff1ec1)
<br>
<br>
The approximate distance of the inner Lagrangian point (L1) from the dwarf star was calculated using equation 5. This was provided to us by Dr. Olenick. We created an “L1point” object to keep track of the location of that point at all times and used the angle between the vectors of the dwarf and donor’s positions to project an updated distance at that angle to give a new location of the L1 point for each iteration of the animation loop. We then used this location as the release point for the particles leaving the donor. 
<br><br>
$$Eq. 5$$
<br>
$$r_{L1}=a[0.500-0.227 \log_{10}\frac{M_2}{M_1}]$$
<br><br>
Each released particle has the approximate radius and mass of  $4 \cdot 10^{14}$ hydrogen atoms. The mass is removed from the donor star when a particle is released. Their release velocity is calculated by first using equation 6 to get the angular velocity, where P is the period of the orbit in seconds. In this case, we ran the simulation once without releasing particles in order to get the period of the simulation and manually entered it into equation 6. We calculated the period using Kepler's Third law. This was to avoid any errors from using the period of the actual orbit as they are not the same. The period of the simulation is 0.050426 days, while the period of the EPIC 220615486 star system is 0.065837 days. Once the angular velocity is found, it is translated to linear velocity by multiplying it by the negative x position, the negative y position, and the z position, respectively. Each of these is then multiplied with a random number between 1.0 and 2.0 to give each particle a slight difference in velocity.
<br><br>
$$Eq. 6$$
<br>
$$\omega=\frac{2\pi}{P}$$
<br><br>
The movement of each particle is calculated with a fourth-order Runge-Kutta in a similar manner to the stars. Since the mass of each particle is negligible, the stars are the only objects which contribute a force to the acceleration.  After the position and velocity are updated for each particle at the end of the Runge-Kutta, that particle object is passed into a function which updates the temperature of that particle based on its velocity and then changes the color of that particle accordingly. The multicolored particles display as a heat map in the accretion disk of the dwarf star (see figure 2). In order to calculate the temperatures, the particles in the accretion disk are treated like an ideal gas. Their temperatures in Kelvin are calculated using the kinetic energy equation for an ideal gas (eq. 7).
<br><br>
$$Eq. 7$$
<br>
$$\frac{3}{2}kT=\frac{1}{2}mv^2$$
<br><br>
<br>
**Fig 2: The Colors of the particles according to heat** (Coolors)
<br>
![TempRange](https://github.com/user-attachments/assets/b5ae4808-bba4-4c72-ba01-82f52f6a7769)
<br>
<br>
We implemented particle collisions by having a function check if a particle is inside the average radius of the dwarf star's accretion disk and, if so, whether that particle has collided with any other particle. The indices of any colliding pair are added to a list which is fed into a separate function that determines their velocity after collision once all the colliding pairs have been found. The change in velocity is found by treating the particle collisions as fully elastic. We used angle-free representations of the equations for final velocity in 2D collisions. Unfortunately, once the simulation is run the change in velocity due to a collision is so drastic that most of the particles shoot out of orbit. 
<br>
<br>
We tried various methods to fix this problem, but the only options that were somewhat successful were implementing gravitational softening or scaling the change in velocity due to a collision. Gravitational softening was introduced in the Runge-Kutta (eq. 8) with $\epsilon= 1000$ $R_\odot$ (since this was the value that had an effect on the particles). This made the particles orbit loosely around the star, but is not consistent with the expected heat-map results in which a particle is hotter the closer it is to the star.  In order to scale the final velocity we simply added ten percent of the change in velocity due to a collision to the initial velocity before a collision (eq. 9). This creates a much smaller change in the correct direction. This gets the expected results (see figures 3.5, 3.6 and 3.7), but needs more evidence to back it up before it can be used as a model.
<br><br>
$$Eq. 8$$
<br>
$$\overset{\rightharpoonup}{a}(t)=-\frac{GM_{star}\overset{\rightharpoonup}{r}}{||\overset{\rightharpoonup}{r}||^3+\epsilon^2}$$
<br><br>
$$Eq. 9$$
<br>
$$\overset{\rightharpoonup}{v}_{new} = \overset{\rightharpoonup}{v}_i+0.1\Delta\overset{\rightharpoonup}{v}$$
<br><br><br>

## 3 RESULTS

<br>

**Fig. 1: EPIC 220615486 orbit with gravitational smoothing**
<br>
![GravSmoothOrbit](https://github.com/user-attachments/assets/fec8ec8c-dd85-459e-a6a6-7a67f6ed539f)
<br>
<br>

**Fig. 2: Particle distance from dwarf star as a function of time (gravitational smoothing)**
<br>
![DistGraphGrav](https://github.com/user-attachments/assets/f1ef7b2f-547e-4816-8e03-826c909c4405)
<br>
<br>
**Fig. 3: Change in particle temperature as a function of time (gravitational smoothing)**
<br>
![TempGraphGrav](https://github.com/user-attachments/assets/af872373-13d2-4cde-8e01-4275303913aa)
<br>
<br>
**Fig. 4: Particle velocity as a function of time (gravitational smoothing)**
<br>
![VelGraphGrav](https://github.com/user-attachments/assets/6fc1530b-934e-4323-860c-c85896ae5837)
<br>
<br>
**Fig. 5: EPIC 220615486 orbit with scaled collision velocity**
<br>
![DeltaVelOrbit](https://github.com/user-attachments/assets/78370bcb-7c88-4f1e-8e92-c6409fc0f0e7)
<br>
<br>
**Fig. 6: Particle distance from dwarf star as a function of time (scaled collision velocity)**
<br>
![DistGraphVel](https://github.com/user-attachments/assets/04a3483a-573b-4e31-b972-c3baab312128)
<br>
<br>
**Fig. 7: Change in particle temperature as a function of time (scaled collision velocity)**
<br>
![TempGraphVel](https://github.com/user-attachments/assets/22f09230-8088-44f5-a723-87258ec931a7)
<br>
<br>
**Fig. 8: Particle velocity as a function of time (scaled collision velocity)**
<br>
![VelGraphVel](https://github.com/user-attachments/assets/2b65b308-a693-40cf-b795-acd754273ca1)
<br><br><br>

## 4 DISCUSSION AND CONCLUSION
Finding the temperature of the dwarf and giant in EPIC 220615486 would enhance the results of this simulation. Rather than figure 3.7  displaying the change in particle temperature as it depends on the velocity change, it would be able to show the change in temperature as a function of the heat radiation from each of the stars in addition to the temperature's dependence on the velocities of the particles.
<br>
<br>
Our simulation relies on scaling down the change in velocity during particle collisions. It would greatly improve the accuracy of this simulation if an alternate way to program collisions were implemented. Multiple methods were discussed during this project. Partially elastic and partially inelastic collisions were proposed. Additionally, viscosity was also considered as a method to slow the particles. This proved unsuccessful. One of last methods of scaling that was attempted was scaling the time-step in the program during the particle collisions. The only method which proved useful was the artificial scaling down of the velocities after collisions (this process is outlined in the procedure). While improvements need to be made, this program succeeds in simulating an accretion disk for EPIC 220615486 and provides the framework for creating a heat-map with accurate particle collisions.
<br> <br> <br>

## ACKNOWLEDGEMENTS
We would like to thank Dr. Olenick for letting us harass him with questions and providing the accretion disk code we built from.
<br> <br> <br>

## BIBLIOGRAPHY
Cataclysimic Variables, Cataclysimic Variables, https://imagine.gsfc.nasa.gov/science/objects/cataclysmic_variables.html
<br>
<br>
Fabrizio Bianchi, Coolors Color Pallet Generator, Coolors.co. Available at: https://coolors.co/. 
<br>
<br>
Montgomery M. M., Olenick R. P., Sweeny A., Lacomb J., Demasi M., Morris W., Smith N., 2024, Monthly Notices of the Royal Astronomical Society
<br> <br> <br>

## APPENDIX
### A. Links to Simulations
- Simulation with Gravitational Smoothing: https://www.glowscript.org/#/user/rgfritz/folder/ParticleCollisonBinaryStars/program/EPIC220615486GravSmoothing 
- Simulation with Scaled Velocity: https://www.glowscript.org/#/user/rgfritz/folder/ParticleCollisonBinaryStars/program/EPIC220615486DeltaVel
- Simulation with Viscocity: https://www.glowscript.org/#/user/rgfritz/folder/ParticleCollisonBinaryStars/program/EPIC220615486Viscocity
<br> <br> <br>
