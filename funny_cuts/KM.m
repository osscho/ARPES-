%close all;
%clear all;

function [x_new, E_surf, ret ] =  KM (workfunction)

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
    %workfunction=4.433;
    %workfunction=4.433+1.205-0.14;
    addpath('./data');
    addpath('./WIEN2k_maps');


    [kx,ky,Ekin,V, angles_x, angles_y]=generate_volumedata('properties_he1_few.txt');




    x1=0.1829; %M
    y1=-1.106;
    x2=-0.5332; %K
    y2=-1.201;



    a=(y2-y1)/(x2-x1);
    b=y1-a*x1;



    steps=1500;


    k_surf_x=meshgrid(linspace(-0.7,0.5,steps));

    %k_surf_y=k_surf_x' * a ;
    %k_surf_y= k_surf_x + b;
    %k_surf_y=k_surf_y';

    x_steps = size(k_surf_x,2);
    y_steps = size(k_surf_x,1);

    E_surf = meshgrid(linspace(9, 18, y_steps))';


    %E_surf=zeros(steps, steps);
    for itx = 1:steps
        for ity = 1:steps
            k_surf_y(itx,ity)=a*k_surf_x(itx,ity)+b;
        end
    end
    


    
    ret=dataslice(kx,ky,Ekin,V, angles_x, angles_y, k_surf_x,k_surf_y,E_surf);
    
    hold on;
    surf( squeeze(kx(:,:,1000)), squeeze(ky(:,:,1000)), squeeze(Ekin(:,:,1000)), squeeze(V(:,:,1000)));
    surf( squeeze(kx(:,:,500)), squeeze(ky(:,:,500)), squeeze(Ekin(:,:,500)), squeeze(V(:,:,500)));
    surf(k_surf_x, k_surf_y, E_surf, ret); 
    colormap('hot');
    shading interp;


    x_new=meshgrid( linspace(0, norm( [ k_surf_x(steps,steps)-k_surf_x(1,1)  k_surf_y(steps,steps)-k_surf_y(1,1) ] ) , steps) );
    E_surf=E_surf+workfunction-E_photon;
    fprintf(' done\n');
end
