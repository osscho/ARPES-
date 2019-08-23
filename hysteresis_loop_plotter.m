clear all;

figure();

for i=1:2
    if i<10 
file_name = sprintf('s100-00%d.dat',i);
    else
file_name = sprintf('s100-0%d.dat',i);
    end
data(:,:) = dlmread(file_name);

plot(data(:,1),data(:,2),'LineWidth',1,'Color',[1 0 0]);
%H = subplot(4,6,i)


%xlabel('Magnetic field [T]');
%ylabel('Signal intensity[V]');
title(file_name);
xlim([-0.01 0.01]);
hold on;
grid on;
   
end