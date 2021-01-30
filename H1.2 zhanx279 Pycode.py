# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:29:47 2021

@author: sophi
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt


flexcode = """TITLE 'H1.2 zhanx279 flexpdePython'    
COORDINATES cartesian2  
VARIABLES      
vx(threshold=0.1)  !velocity in x
vy(threshold=0.1) 	!velocity in y
xd(threshold=0.1) 	!x-displacement
yd(threshold=0.1)   !y-displacement
SELECT     
ngrid = 1 !since we're not using spatial depenedence within the object we don't need a dense mesh
DEFINITIONS  
g  = 9.81									!gravity constant
theta_i = 50*pi/180				 !initial angle (50)
mi=4000									!mass (kg)
percentFuel= 0.85					 !percent of rocket mass that is fuel
mfi=percentFuel*mi			!initial fuel mass (80 of mass)
mrocketNotfuel=mi-mfi 	 !or just mi*(1-percentFuel)
vi= 80 										!initial velocity
q=%s 										!rate that rocket expends fuel (kg/s) (flow rate)
v=sqrt(vx^2+vy^2) 				!velocity calculated from components vx vy
tfuel=mfi/q
vfuel =%s 							!fuel velocity as it leaves the rocket (m/s)

!PART b
!vwindx = 	if (yd>4000) then 0 else -5				!wind vel x
vwindx = -5
vwindy = -10						!wind vel y
CD=0.5 								!drag coeff
A=2											!frontal area (m^2)
!rho=1.225
!rho=1.225-0.000113*yd							!Air density....this changes in relation to the rockets altitude..
rho= if (yd<1000) then
 1.225 else if (yd<2000) then 0.112 
else if (yd <3000) then 1.007
 else if (yd<4000) then 0.9093
 else 0.8194	!Air density....this changes in relation to the rockets altitude..

vrelx = vx-vwindx
vrely = vy-vwindy
vrel = sqrt(vrelx^2+vrely^2)

!PART b end

mfuel=if(t<tfuel) then mfi-q*t else 0

m=mfuel+mrocketNotfuel

Ft=if(t<tfuel) then q*vfuel else 0		!thurst force only if rocket still has fuel left

KE = 0.5*m*v^2 !kinetic energy

Fd=0.5*rho*CD*A*vrel^2	!drag froce calculation
!drag force broken down to components
Fdx=Fd*vrelx/vrel		!drag force x
Fdy=Fd*vrely/vrel 	!drag force y

Ftx = Ft*vx/v		!thrust force x
Fty = Ft*vy/v		!thrust force y

Fx= Ftx-Fdx
Fy= Fty -Fdy

ax = (Fx)/m 		!accel x with thrust force and drag force
ay = (Fy)/m -g		!accel y with thrust force and drag force



INITIAL VALUES
vx = vi*cos(theta_i)
vy = vi*sin(theta_i)
xd = 0
yd = 0
EQUATIONS       
vx: dt(vx) = ax
vy: dt(vy) = ay
xd: dt(xd) = vx
yd: dt(yd) = vy
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 200 halt (yd<-1000) ! if time dependent }
PLOTS           
for t = 0 by 1 to endtime
history(yd) at (0,0) vs xd
history(ax,ay) at (0,0)
history(vx,vy) at (0,0) 

SUMMARY
report val (xd, 0,0) as 'x range'
report val(KE,0,0) as 'Kinetic energy at impact'
report val(vx,0,0) as 'vx final'
report val(vy,0,0) as 'vy final'


history(xd,yd,KE,vx) at (0,0) PrintOnly Export Format '#t#b#1#b#2#b#3#b#4' file = 'filetext.txt'
END


"""


FlexfileName = 'H1.2 zhanx279.pde'
# angles = np.arange(10,21,10)
# veli = np.arange(30,9,-10)

flowrate = np.arange(70,0,-10)
indx = [0,1,2,3,4,5,6]
fuelvel = np.arange(2500,6801,300)


max_range = -1
best_flowrate = -1
best_fuelvel =-1

for i in indx:   
    with open(FlexfileName, "w") as f:
        print(flexcode%(flowrate[i],fuelvel[i]), file=f)
    subprocess.run(["FlexPDE6s","-S",FlexfileName], timeout=20)
    with open(FlexfileName, "r") as f:
        flexoutputrawdata = np.loadtxt('filetext.txt', skiprows = 7)
    t = flexoutputrawdata[:,0]
    xd = flexoutputrawdata[:,1]
    yd = flexoutputrawdata[:,2]
    KE = flexoutputrawdata[:,3]
    print('\nFlow rate of {flow} kg/s and initial velocity'.format(flow = flowrate[i]),
          ' of {veli} m/s landed after {t} sec at x = {x} m'.format(veli=fuelvel[i],
                                                                    t=t[-1],x=xd[-1]))
    #print('Kinetic energy on impact with angle of {angle}\N{DEGREE SIGN} is {ke} J'.format(angle=flwrte,ke=KE[-1]))
    #the -1 gives last elemtn in array
    
    x_range=xd[-1]
    
    if x_range > max_range:      
        max_range = x_range       
        best_flowrate = flowrate[i]
        best_fuelvel = fuelvel[i]


    plt.plot(xd,yd,label=[flowrate[i],fuelvel[i]])
        

plt.xlabel('x-displacement, [m]')
plt.ylabel('y-displacement, [m]')
plt.legend(loc='upper right', ncol = 1)
plt.show()

print("\nMaximum range: ", max_range, "m achieved with flow rate of", 
      best_flowrate, 'kg/s and fuel velocity of', best_fuelvel, 'm/s')  