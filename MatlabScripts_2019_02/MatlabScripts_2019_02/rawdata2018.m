clear all;
% close all;

file = 283;         %CHANGE
region=0;          %CHANGE


name = sprintf('SPE%05d_%05d.txt',file,region);        
data = dlmread(name,'\t',57,0);

figure(12)
imagesc(data)


