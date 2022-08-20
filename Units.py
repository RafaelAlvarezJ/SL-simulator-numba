# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 15:23:47 2020

@author: GTB
"""

"""
 The units.m function returns a struct containing the SI values of 
 many common units. The example below demonstrates how to use 
 this struct to cleanly and effeciently perform matlab calculations 
 using commonly encountered units and physical constants. The code 
 is easily modified to include any non-standard unit or constant desired.

 To determine the exact syntax of a particular unit you would like to use,
 simply run the function with no semicolon and all the units will be
 displayed.

  Using the units struct:
  --------------------------------------------------------
    First create a units struct by including in your code the following
    line
          u = units;
    Then:
    
          %To enter a number in a given unit, MULTIPLY by the unit:
                L = 5*u.in   % matlab automatically displays L in SI units
  
          % To display in a desired unit, simply divide by that unit
                L/u.ft       % displays L in ft.
  
          % To convert between units, MULTIPLY by the starting unit, and
          % DIVIDE by the ending unit:
                u.mi^2/u.ft^2  %displays the number of feet^2 in one mile^2
  
          %More complicated units can be obtained through arithmatic
                mach1 = 340.29*u.m/u.s;  %speed of sound 
  
             %Note... to make the speed of sound available wherever your
             %units struct is defined, simply write:
                u.mach1 = 340.29*u.m/u.s;   %mach1 now part of units struct
 


 ------  BEGIN EXAMPLE CODE --------------------------------
 This is an example calculation that uses the units mfile to calculate the
 pressure at the bottom of a long vertically oriented pipe that is capped 
 at the bottom and filled with oil.

 u = units;
 pipeInnerDiameter = 4*u.in;     %4 inch inner diameter
 pipeHeight = 30*u.ft;           %pipe sticks 30 feet up into the air
 densityOfOil = 0.926*u.gm/u.cc; %density of oil as found on some random web site = .926 gm/cc
 pipeCrossSectionArea = pi*(pipeInnerDiameter/2)^2;  %cross sectional area of pipe bore
 volumeOfOil = pipeCrossSectionArea * pipeHeight;    %volume of oil that the pipe can hold
 pressurePipeBottom = densityOfOil * u.g * pipeHeight;  %pressure formula from physics: P = rho*g*h.
 forceOnPipeBottom = pressurePipeBottom * pipeCrossSectionArea;  %force exterted on bottom cap of the pipe.
 
 %Note that each variable holds its value as expressed in SI units.  To
 %express these values in different units, simply divide by the desired
 %unit as shown below.
 line1 = sprintf('A %2.3g inch diameter pipe sticking %3.3g meters into the air and filled',pipeInnerDiameter/u.in, pipeHeight/u.m);
 line2 = sprintf('with %3.3g fluid ounces of oil will have a pressure at the bottom of %4.4g psi.',volumeOfOil/u.floz, pressurePipeBottom/u.psi);
 line3 = sprintf('This will cause a total force of %5.5g lbs to press on the bottom cap of the pipe.',forceOnPipeBottom/u.lbf);
 
 textVal = sprintf('\n\n%s\n%s\n%s\n',line1,line2,line3);
 disp(textVal);
------  END EXAMPLE CODE --------------------------------

"""

#============ START THE ACTUAL CODE TO DEFINE THE UNITS STRUCT =========
#-------- UNITS ------------------------------
#------- length ----
m = 1
km = 1e3*m
cm = 1e-2*m
mm = 1e-3*m
inch = 2.54*cm
ft = 12*inch
yd = 3*ft
#------- Area -------
sqm=1
sqcm=(cm)**2
sqkm=(km)**2
sqft=(ft)**2
D=9.869233e-13*sqm
mD=1e-3*D
#------- Volume -------
cc = (cm)**3
L = 1000*cc
mL = cc
gal = 3.78541197*L
bl = 42*gal

#----- mass ---------
kg = 1
gm = 1e-3*kg
mg = 1e-3*gm
lb = 0.45359237*kg

#---- time -------
sec = 1
minu = 60*sec
hr = 60*minu
day = 24*hr
yr = 365.242199*day 

#---- force -------
N = 1
dyne = 1e-5*N
lbf = 4.44822*N

#----- energy -----
J = 1
MJ = 1e6*J
kJ = 1e3*J
mJ = 1e-3*J
uJ = 1e-6*J
nJ = 1e-9*J
eV = 1.6022e-19*J
BTU = 1.0550559e3*J
kWh = 3.6e6*J
cal = 4.1868*J
kCal = 1e3*cal

#---- temperature ---
K = 1
mK = 1e-3*K
uK = 1e-6*K
nK = 1e-9*K

#---- pressure -----
Pa = 1
MPa = 1e6*Pa
bar = 1e5*Pa
atm = 1.013e5*Pa
psi = 6.895e3*Pa

# Viscosity
PaS=1
Cp=1e-3*PaS

#----fundamental constants ----
g = 9.80665*m/sec**2;



