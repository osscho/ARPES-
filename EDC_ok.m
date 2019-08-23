close all;
clear all;
hold on;

file_start=2818;
file_end=2819;
cmap = colormap(jet); 
for i=file_start:file_end
    f=i;

file_name = sprintf('SPE%d_09.txt',f);
data_EDC(:,:) = dlmread(file_name,'\t',82,0); %1000X256
size_data=size(data_EDC);
EDC_sum=zeros(1,size_data(1));
for j=1:size_data(1)
 for i=2:256
 EDC_sum(1,j)=EDC_sum(1,j)+data_EDC(j,i);
 end;
end;


g=(i-1914)*0.0001;
plot(data_EDC(1,1),EDC_sum(1,1),data_EDC(2,1),EDC_sum(1,2),'r', 'LineWidth',1)
title('ANGLE 950'); xlabel('Kinetic energy [eV]'); ylabel('Counts [arb.u]');
legend('350','260');
hold on;
end;

