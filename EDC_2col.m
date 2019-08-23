close all;
clear all;
hold on;

file_start=188;
file_end=188;
cmap = colormap(jet); 
for i=file_start:file_end
    f=i;

file_name = sprintf('SPE0%d_10.txt',f);
data_EDC(:,:) = dlmread(file_name,'\t',82,0); %1000X256

g=(i-1914)*0.0001;
plot(data_EDC(:,1),data_EDC(:,2),'r', 'LineWidth',1)
title(' '); xlabel('Kinetic energy [eV]'); ylabel('Counts [arb.u]')
hold on;
end;

