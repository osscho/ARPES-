%critical current for spin accumulation/SOT
clear all;
close all;
I_density_crit=6*10^9; %critical current density in [A/cm2]
%ro=1; %resistivity in [ohm*cm]
w=0.5*10^-4; %stripe width in [mm]
d=3.6; %film thickness in [nm]
A=w*10^-1*d*10^-7; %area in [cm2]
I_crit=I_density_crit*A; %critical current [A]
fprintf('Critical current [A] is:%f', I_crit);

R=10^2;%ohm;
U=R*I_crit;
fprintf('Critical voltage [V] is:%f', U);
