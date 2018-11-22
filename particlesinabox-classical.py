# vpython simulation - particles in a box - classical collisions
#NOTES:
#Note that with deterministic momentum the collisions will repeat themselves, if there is an element of randomness, 
#an arbitary number of new states can be reached! The entire collision space could hypothetically be sampled!
from vpython import *
from scipy.stats import maxwell
import random
import numpy as np
   
def generate_particles(n,box_x,box_y,box_z,particle_radius):
	#Generate balls in a grid like pattern starting from the bottom
	#of the box so that none of the balls are overlapping at the 
	#start of the simulation
	array = []
	box_length = box_x # since all 3 directions are the same length
	e = box_x/100 # this is epsilon - the length that separates each ball
	num_balls = box_length / (2 * particle_radius + 2 * e)
	num_balls = int(num_balls) # round down to be on the safe side
	counter_y = 0
	counter_z = 0
	inc = particle_radius + e
	x_pos = -box_x/2
	y_pos = -box_x/2 + inc
	z_pos = -box_x/2 + inc
	x_vel_dist = maxwell.rvs(size=n)
	y_vel_dist = maxwell.rvs(size=n)
	z_vel_dist = maxwell.rvs(size=n)
	#The following loop generates particles that are evenly spaced
	for i in range(n):
		x_vel = x_vel_dist[i]
		y_vel = y_vel_dist[i]
		z_vel = z_vel_dist[i]
		x_pos += inc
		particle = sphere( pos=vector(x_pos,y_pos,z_pos),radius=particle_radius, color=color.white,
		velocity=vector(x_vel,y_vel,z_vel), mass=1) 	
		array.append(particle)
		counter_y += 1
		counter_z += 1
		if counter_z % num_balls == 0:	
			z_pos += 2 * inc
			x_pos = -box_x/2 - inc
		if counter_y % num_balls ** 2 == 0:	
			y_pos += 2 * inc
			z_pos = -box_x/2 + inc
		x_pos += inc
	return array
		
def boundary_collision(particle,particle_radius,box):
    if (abs(particle.pos.x) + particle_radius  >= box.size.x/2):
        particle.velocity.x = -particle.velocity.x
    elif (abs(particle.pos.y) + particle_radius >= box.size.y/2):
        particle.velocity.y = -particle.velocity.y
    elif (abs(particle.pos.z) + particle_radius >= box.size.z/2):
        particle.velocity.z = -particle.velocity.z

def particle_collision(particle1,particle2,dt):
    if ((particle1.pos-particle2.pos).mag < particle1.radius + particle2.radius):
        angle = random.uniform(0,2 * np.pi)
        rand_x = random.uniform(0,1)
        rand_y = random.uniform(0,1)
        rand_z = random.uniform(0,1)
        n = (particle1.pos - particle2.pos).norm() # m-axis
        a1 = particle1.velocity.dot(n)
        a2 = particle2.velocity.dot(n)
        P = (a1-a2)  
        particle1.velocity = particle1.velocity - P * n
        particle2.velocity = particle2.velocity + P *  n
			
scene = canvas(title='Particles',
     width=1430, height=710,
     center=vector(0,0,0), background=color.black) # configure the vypython scene

number_of_particles = 20 
box_x = 20 # length of box in x-direction
box_y = 20 # length of box in y-direction
box_z = 20 # length of box in z-direction
particle_radius = 1
box = box(pos=vector(0,0,0), size=vector(box_x,box_y,box_z), 
	color=color.white, opacity=0.1)			
			
#Generate array of particles
particles = generate_particles(number_of_particles,box_x,box_y,box_z,particle_radius)
len_particles = len(particles)
k = 0 # some value that determines when the simulation stops
dt = 0.01 # time increment
while (k >= 0):
    rate(1500)
    if k % 100 == 0:
        avg_x_vel = 0 # average squared x-velocity
        avg_y_vel = 0 # average squared y-velocity
        avg_z_vel = 0 # average squared z-velocity
        #Computing pressure and temperature for the box
        for p in particles:
            avg_x_vel += (p.velocity.x ** 2)
            avg_y_vel += (p.velocity.y ** 2)
            avg_z_vel += (p.velocity.z ** 2)
        avg_x_vel = avg_x_vel/number_of_particles
        avg_y_vel = avg_y_vel/number_of_particles
        avg_z_vel = avg_z_vel/number_of_particles
        tot_avg_square_velocity = avg_x_vel + avg_y_vel + avg_z_vel
        volume = box_x * box_y * box_z
        mass = 1
        pressure = (number_of_particles * mass / volume) * tot_avg_square_velocity
        temperature = (pressure * volume ) / (number_of_particles * 1.38 * (10 ** -23))
        print("The pressure is: " + str(pressure))
        print("The temperature is: " + str(temperature))
        print("The average x_vel squared is: " + str(avg_x_vel))
        print("The average y_vel squared is: " + str(avg_y_vel))
        print("The average z_vel squared is: " + str(avg_z_vel))
        print()
    #Checking for collisions and updating each particle's position
    for i in range(len_particles):
        particles[i].pos = particles[i].pos + particles[i].velocity * dt
        boundary_collision(particles[i],particle_radius,box)
        for j in range(0,len_particles):
            if j == i:
                continue
            particle_collision(particles[i],particles[j],dt)
            boundary_collision(particles[i],particle_radius,box)
            boundary_collision(particles[j],particle_radius,box)
    k = k + 1    
      
      
	
