close all; clear all;

E_photon = 21.22;
f=297;


    file_name = sprintf('SPE0%d_10.txt',f);       %CHANGE !    
    data_1(:,:) = dlmread(file_name,'\t',82,0);
    for j=1:256
        kinetic_energy(:,j) = data_1(:,1);
        binding_energy(:,j) = E_photon - 4.423 - kinetic_energy(:,j);
    end;
    figure(1);
    image(real2rgb(data_1(:,:),'temp',[50 500]));