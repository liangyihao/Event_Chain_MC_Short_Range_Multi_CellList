#!/usr/bin/env python
# encoding: utf-8

########################################################
# This code is to generate input ECMC code for spheres.
# Two types of particles: A and B
# Interactions: hard sphere & Gauss
# Written by: Boran Ma (boran.ma@northwestern.edu)
# Update: 06/25/19
# Note: integer division has changed to // in Python 3+
########################################################

import math

#################
diamA = 10.0
diamB = 2.0
lx=61.0
ly=61.0
lz=61.0

volFrac = 0.55
numA = 64

hsAA = diamA
hsBB = diamB
hsAB = 0.5*(diamA+diamB)


Vbox=lx*ly*lz
Veff=Vbox-numA*4.0*math.pi/3.0*((hsAB)**3)
numB = int(Veff*volFrac/(4.0*math.pi/3.0*((diamB/2)**3)))

epsilon = -5.0
sigma = 10.0
rcut = 10.0
timestep = 1e4
distance = 10.0
skipFrame = 1e3
stride = 1
expandCoeff = 1.02
expandCoeff2 = 1.0001
#################

# start writing input file
fileWrite = open('in_A.ecmc' , 'w')
fileWrite.write('Lx %f\n' %lx)
fileWrite.write('Ly %f\n' %ly)
fileWrite.write('Lz %f\n' %lz)
fileWrite.write('\n')
fileWrite.write('Type A\n')
fileWrite.write('Type B\n')
fileWrite.write('\n')
fileWrite.write('Interaction-Type-Type A A Hard %f\n' %hsAA)
fileWrite.write('Interaction-Type-Type A B Hard %f\n' %hsAB)
fileWrite.write('Interaction-Type-Type B B Hard %f\n' %hsBB)
fileWrite.write('Interaction-Type-Type A B Gauss %f %f %f\n' %(sigma, epsilon, rcut))
fileWrite.write('\n')

# position section
fileWrite.write('Positions\n')

# write type A beads first
x_A=[]
y_A=[]
z_A=[]
nx = int(lx/(diamA*expandCoeff))
ny = int(ly/(diamA*expandCoeff))
nz = int(lz/(diamA*expandCoeff))
hx = lx/nx
hy = ly/ny
hz = lz/nz
A_count = 0
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            if A_count<numA:
                x_coord=i*hx
                y_coord=j*hy
                z_coord=k*hz
                mark=1
                for l in range(len(x_A)):
                    dx=x_coord-x_A[l]
                    dy=y_coord-y_A[l]
                    dz=z_coord-z_A[l]
                    if dx>lx/2:
                        dx=dx-lx
                    if dx<-lx/2:
                        dx=dx+lx

                    if dy>ly/2:
                        dy=dy-ly
                    if dy<-ly/2:
                        dy=dy+ly
                    
                    if dz>lz/2:
                        dz=dz-lz
                    if dz<-lz/2:
                        dz=dz+lz
                    
                    dr=(dx*dx+dy*dy+dz*dz)**0.5
                    if dr<hsAA*expandCoeff2:
                        mark=0
                if mark==1:
                    A_count=A_count+1
                    print('A:',A_count)
                    x_A.append(x_coord)
                    y_A.append(y_coord)
                    z_A.append(z_coord)
                    fileWrite.write('A %f %f %f\n' %(x_coord, y_coord, z_coord))

# write type B beads next
x_B=[]
y_B=[]
z_B=[]
#nx = int(lx/(diamB*expandCoeff))
#ny = int(ly/(diamB*expandCoeff))
#nz = int(lz/(diamB*expandCoeff))
hx=diamB*expandCoeff
hy=0.86603*diamB*expandCoeff
hz=0.81650*diamB*expandCoeff
nx=int(lx/hx)
ny=int(ly/hy)
nz=int(lz/hz)
B_count = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            if B_count<numB:
                offset_x=0.5*hx*(j%2)
                offset_y=(k%2)*(2/3*hy)

                x_coord=offset_x+i*hx
                y_coord=offset_y+j*hy
                z_coord=k*hz
                mark=1
                for l in range(len(x_A)):
                    dx=x_coord-x_A[l]
                    dy=y_coord-y_A[l]
                    dz=z_coord-z_A[l]
                    if dx>lx/2:
                        dx=dx-lx
                    if dx<-lx/2:
                        dx=dx+lx

                    if dy>ly/2:
                        dy=dy-ly
                    if dy<-ly/2:
                        dy=dy+ly
                    
                    if dz>lz/2:
                        dz=dz-lz
                    if dz<-lz/2:
                        dz=dz+lz

                    dr=(dx*dx+dy*dy+dz*dz)**0.5
                    if dr<hsAB*expandCoeff2:
                        mark=0
                        break

                for l in range(len(x_B)):
                    if mark==0:
                        break
                    dx=x_coord-x_B[l]
                    dy=y_coord-y_B[l]
                    dz=z_coord-z_B[l]
                    if dx>lx/2:
                        dx=dx-lx
                    if dx<-lx/2:
                        dx=dx+lx

                    if dy>ly/2:
                        dy=dy-ly
                    if dy<-ly/2:
                        dy=dy+ly
                    
                    if dz>lz/2:
                        dz=dz-lz
                    if dz<-lz/2:
                        dz=dz+lz

                    dr=(dx*dx+dy*dy+dz*dz)**0.5
                    if dr<hsBB*expandCoeff2:
                        mark=0
                        break
                if mark==1:
                    B_count=B_count+1
                    print('B:',B_count,' of ',numB,' created')
                    print(i,',',j,',',k)
                    x_B.append(x_coord)
                    y_B.append(y_coord)
                    z_B.append(z_coord)
                    fileWrite.write('B %f %f %f\n' %(x_coord, y_coord, z_coord))

fileWrite.write('End-Positions\n')
fileWrite.write('\n')

# ECMC run section
fileWrite.write('loop %d\n' %timestep)
fileWrite.write('\tECMC +x %f\n' %distance)
fileWrite.write('\tECMC +y %f\n' %distance)
fileWrite.write('\tECMC +z %f\n' %distance)
fileWrite.write('\tECMC -x %f\n' %distance)
fileWrite.write('\tECMC -y %f\n' %distance)
fileWrite.write('\tECMC -z %f\n' %distance)
fileWrite.write('\tOut %d %d\n' %(skipFrame, stride))
fileWrite.write('End-loop\n')

