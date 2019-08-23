%Oersted field calculated using the Biot Savart law, according to A.J.Phys.
%68, 254 (2000). Stright wire 
clear all;
close all;
d=10*10^-6 % distance from the wire center to the point of interest in [m]
w=500*10^-6; %stripe lenght in [m]
I=0.06% current in [A]
betha=180/pi*atan((w/2)/d)+90
alpha=90-180/pi*atan((w/2)/d)
mi0=1.2566*10^-6% vaccum permeability in [T*m/A] 
B=mi0*I/(4*pi*d)*(cos(alpha*pi/180)-cos(betha*pi/180))
fprintf('Magnetic field [mT] is:%f', B*1000);