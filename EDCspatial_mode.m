clear all;
file_start=1915;
file_end=1926;
E_photon=21.22;
%Ep_Xe = 8.44;
%Ep_Ne = 16.848;
%Ep_HeI = 21.22;
%Ep_HeII = 40.8;

figure();
for i=file_start:file_end
    file_name = sprintf('SPE%d_11.txt',i);    
    data(:,:) = dlmread(file_name,'\t',82,0); %1000x256
   % ed=zeros(numel(data(:,1)),1);
   % for l=1:numel(data(:,1));
   % for j=2:256    
   % ed(l,1)=ed(l,1)+data(l,j);
    end;
    
    figure();
    kinetic_energy(:,1) = data(:,1);
    binding_energy(:,1) = -(E_photon - 4.423 - kinetic_energy(:,1));
    plot(kinetic_energy(:,1),data(:,1));
    title('2072 [grounded]');
    %xlim([-0.5 0.5]);
    xlabel('Kinetic energy [eV]');
    ylabel('Intensity [arb.u.]');
    
end;