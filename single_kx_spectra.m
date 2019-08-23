addpath('Z:\Matlab scripts');
addpath('Z:\Matlab scripts\real2rgb');
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultAxesFontName', 'Arial')

clear all;
C = {'r','k','b','g','y','m','c',[.5 .6 .7],[.8 .2 .6],[.1 .2 .6],[.8 .2 .1],[.8 .8 .6]} % Cell array of colors.
E_photon = 21.22;
wf=4.423;


l=[1408];
%check: RIGHT2 5971
%l=[l1,l2,l3,l4];
s=size(l);
s1=s(2);
for i=1:s1
    f=l(i);
        
    file_name = sprintf('SPE%d_06.txt',f);       %CHANGE !    
    data_1(:,:,i) = dlmread(file_name,'\t',82,0);
    for j=1:256
        kinetic_energy(:,j,i) = data_1(:,1,i);
        binding_energy(:,j,i) = E_photon - wf - kinetic_energy(:,j,i);
        data_norm(:,:,i)=data_1(:,:,i)./max(max(data_1(:,:,i)));
        

    end;
end;
%smoothing data 
 
for i=1:s1
for j=2:256
      data_smooth(:,j,i) = fastsmooth(data_norm(:,j,i),2,1,0);
  end; 
end;
%EDC
size_data=size(data_1);
EDC_sum=zeros(s1,size_data(1));
for k=1:s1
for j=1:size_data(1)
 for i=1:230 %49-183 He; 74-158 Xe
 EDC_sum(k,j)=EDC_sum(k,j)+data_1(j,i,k);
 EDC_sum_p=EDC_sum';

 end;
end;
end;

norm=600;
       
  center_x_pixel=130; %change
  detector_angle=30.0; %0d8; to be checked
 %defining kx scale
 for j=1:s1 
 for i = 1:256
        k_x(:,i,j) = 0.5123*sqrt(kinetic_energy(:,1,j))*sin(((i - center_x_pixel)*detector_angle/256)*2*pi/360);    
  end;
 end;
 
figure();
 surf(k_x(:,:,1)*0,k_x(:,:,1),-binding_energy(:,:,1),real2rgb((data_smooth(:,:,1)),'temp',[0 0.8]),'edgecolor','none');view(90,0); ylim([-0.5,0.5]); box on;
 zlim([-0.3,0.02])
 ylabel('kx [1/A]');zlabel ('binding energy [eV]');title('Theta=962 deg');
 
 
