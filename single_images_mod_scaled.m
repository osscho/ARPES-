close all; clear all;
E_photon = 21.22;
l=3962:1:3962;
%l1=3109:1:3139;
%l2=3164:1:3180;
%l3=3206:1:3214;
%l4=3078:1:3108;
%l=[l1,l2,l3,l4];
s=size(l);
s1=s(2);
for i=1:s1
    f=l(i);
    file_name = sprintf('SPE%d_10.txt',f);       %CHANGE !    
    data_1(:,:) = dlmread(file_name,'\t',82,0);
    for j=1:256
        kinetic_energy(:,j) = data_1(:,1);
        binding_energy(:,j) = E_photon - 4.423 - kinetic_energy(:,j);
    end;
    figure(1);
    
       
  center_pixel=116; %change
  detector_angle=43.7; %0d8; to be checked
  k_range=0.5124*sqrt(16)*sin(detector_angle*pi/180);
  k_step=k_range/256;
  k=[1:1:256];
  k_scaled=(k-116)*k_step;
  image(k_scaled,data_1(:,1)+1.162,real2rgb(data_1(:,:),'temp',[20 100]));
       
 % file_save=sprintf('EM14_%d_scaled.png',f);
  
  % a=real2rgb(data_1(:,:),'temp',[0 80]);
  % imwrite(a,file_save,'Compression','None');
   
   %norm=400;
  
    % if f>3108 && f<3140 %l1
      %  norm=8000
     %end;
    
    % if f>3163 && f<3181 % l2
    %    norm=7000 
    % end;
    
    % if 3205<f && f<3215 %l3 ok
    %    norm=5000
   %  end;
    
   %  if f<3109 % l4
   %     norm=12000   
  %   end;
 % a=real2rgb(data_1(:,:),'temp',[10 norm]);  
 % if i==1
 % imwrite(a,'em14_FS_RT.tif','Compression','None');
 %  else
 %  imwrite(a,'em14_FS_RT.tif','Compression','None','WriteMode','append');
 %  end;
end;

