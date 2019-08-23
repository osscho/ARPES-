function [kx, ky, energies, v, angles_x, angles_y] = generate_volumedata(properties_file);

workfunction=4.433;


fileproperties = dlmread(properties_file); 


%EF = 4.2;
%Ep_Xe = 8.44;
%Ep_Ne = 16.848;
Ep_HeI = 21.21802;
Ep_HeIIa = 40.81303;
Ep_HeIIb = 40.81376;
Ep_HeII = (Ep_HeIIa+Ep_HeIIb)/2;

E_photon = Ep_HeI;

startenergy= 10;
endenergy=15.5;
energystep=0.05;


        %x-angles
anglecutoff_left=41;
anglecutoff_right=34;
numberofangles=255-anglecutoff_left-anglecutoff_right;
%center_pixel = 156-anglecutoff_left;  
%detector_angle = 40;  
%anglestep= 18.4/111*0.86;
center_pixel = 256/2 - anglecutoff_left + 20;  
%detector_angle = 45;
detector_angle = 43;
anglestep= detector_angle/256;


%y_angle_offset= 7.2;
y_angle_offset=-5;

%normalize to 1 count/pixel
normalizearea=1; 
normalizelinebyline=0;



numberoffiles=size(fileproperties,1);

for fileindex = 1:numberoffiles
        file_name = sprintf('SPE%04d_%d.txt',fileproperties(fileindex,1), fileproperties(fileindex,3));
        data(fileindex,:,:) = dlmread(file_name,'\t',82,0);
        
end 

        %normalize area

if(normalizearea)
        for fileindex = 1:numberoffiles
                s = 0;
                for as = 2: 2+anglecutoff_left: size(data,3)-anglecutoff_right
                        
                        for es = 1:size(data,2)
                                s= s+data(fileindex, es,as);
                        end
                end
                
                s =s/(numberofangles+size(data,2));
                for as = 2: size(data,3)
                        for es = 1:size(data,2)
                                data(fileindex, es,as) = data(fileindex, es,as)/s;
                        end
                end
        end
end



        %normalize line by line

if(normalizelinebyline)
        for fileindex = 1:numberoffiles
                for as = 2+anglecutoff_left: size(data,3)-anglecutoff_right
                        s = 0;
                        for es = 1:size(data,2)
                                s= s+data(fileindex, es,as);
                        end
                        if(s)
                                data(fileindex, :,as) = data(fileindex, :,as)*size(data,2)/s;
                        end
                end
        end
end

numberofenergies=size(data,2);

kx=zeros(numberofangles,numberoffiles,numberofenergies);
ky=zeros(numberofangles,numberoffiles,numberofenergies);
v=zeros(numberofangles,numberoffiles,numberofenergies);

for fileindex = 1:numberoffiles
    angles_y(fileindex)=fileproperties(fileindex,2)-y_angle_offset;
    for ite = 1:numberofenergies;
        for itx = 1:numberofangles
                angles_x(itx)=(itx-center_pixel )* anglestep;
                v(itx, fileindex, ite) = v(itx, fileindex, ite)+data(fileindex, ite,1+itx+anglecutoff_left);
                energies(itx, fileindex, ite) = data(fileindex,ite, 1);
                kx(itx, fileindex, ite) = 0.5123*sqrt(energies(itx, fileindex, ite))*sin( (itx-center_pixel )* anglestep *pi/180);
                %ky(itx, fileindex, ite) = 0.5123*sqrt(energies(itx, fileindex, ite))*cos( (itx-center_pixel )* anglestep *pi/180) *sin( (fileproperties(fileindex,2)-y_angle_offset) *pi/180 );    
                ky(itx, fileindex, ite) = 0.5123*sqrt(energies(itx, fileindex, ite))*cos( (itx-center_pixel )* anglestep *pi/180) *sin( (fileproperties(fileindex,2)-y_angle_offset) *pi/180 );    
        end            
    end
end

clear data;


end