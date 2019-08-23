
%clear all;

function [k_surf_x, E_surf, ret ] =  other_cuts (workfunction)
    

    
    fprintf('reading calculated Data...');
    switch nargin
        case 0
            %workfunction=4.433+1.12;
            workfunction=4.433;

    end

    %EF = 4.2;
    %Ep_Xe = 8.44;
    %Ep_Ne = 16.848;
    Ep_HeI = 21.21802;
    Ep_HeIIa = 40.81303;
    Ep_HeIIb = 40.81376;
    Ep_HeII = (Ep_HeIIa+Ep_HeIIb)/2;

    E_photon = Ep_HeI;


    [kx,ky,Ekin,V, angles_x, angles_y]=generate_volumedata('properties_he1_few.txt');

    
    
    %creating a surface

    
    %dimensions of the mesh (density)
    dim_A=50;
    dim_B=40;
    
    %define boundaries
    x_min = -1;
    x_max = 1;
    y_min = -1.2;
    y_max = -1.2;
    E_min = 11;
    E_max = 16;
    
    
    %!!create the first two meshes. the parameter that doesn't change (y in this case) cannot be here, or the mesh will be only a line
    [k_surf_x, E_surf] = meshgrid ( linspace(x_min, x_max, dim_A), linspace(E_min, E_max, dim_B));
    
    %the thid parameter
    k_surf_y = meshgrid ( linspace(y_min, y_max, dim_A), linspace(y_min, y_max, dim_B));


    subplot(1,2,1);
    hold on;
    
    
    %constant energy cuts:
    surf( squeeze(kx(:,:,1000)), squeeze(ky(:,:,1000)), squeeze(Ekin(:,:,1000)), squeeze(V(:,:,1000)));
    surf( squeeze(kx(:,:,500)), squeeze(ky(:,:,500)), squeeze(Ekin(:,:,500)), squeeze(V(:,:,500)));
    
    %the mesh itself just with a constant color shows, where your surface lays with respect to your data
    mesh(k_surf_x, k_surf_y, E_surf, E_surf*0 + 0.5); 
    shading interp;
    grid on;
    
    
    
    subplot(1,2,2)
    
    %making the cut:
    ret=dataslice(kx,ky,Ekin,V, angles_x, angles_y, k_surf_x,k_surf_y,E_surf);      
    
    surf(k_surf_x, k_surf_y, E_surf, ret); 
    colormap('hot');
    shading interp;


    fprintf(' done\n');
end
