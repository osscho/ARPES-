
% RUN LIST OF FILE GENERATOR FIRST!!!!

close all; 
clear all;
set(0, 'DefaultFigureRenderer', 'zbuffer');
%constants;
%Ep_HeI = 21.22;    %Ep_Laser = 6;  %Ep_Xe = 8.44;  %Ep_Ne = 16.848;    %Ep_HeII = 40.8;

file=263;   %CHANGE

EF = 1.638;                                                 

center_x_pixel = 25;                                     % Zentral Pixel in x-Richtung - definiert die 0 auf der kx-Skala
center_y_pixel = 241;                                    % Zentral Pixel in y-Richtung - definiert die 0 auf der ky-Skala                                            
                                            

                                
colorm = 'temp';                                            % die Colormap
% colorm = 'jet'; 
% colorm = 'hot';
% colorm = flipud(gray);
lower = 1 ;                                                % definiert die untere und obere Grenze der Farbskala (1.1*unterster Wert - 0.9*oberster Wert)
upper = 0.7;

%--------read data: ------------------------------------------------------- 
list_of_files = load(sprintf('list_of_files_%05d.txt',file));                   
number_of_files = size(list_of_files,1); 

name = sprintf('SPE%05d_00000.txt',file); 
[step_name, step]= textread(name,'%s %f',1,'headerlines',52);
[angle_max_name, angle_max]= textread(name,'%s %f',1,'headerlines',53);
[angle_min_name, angle_min]= textread(name,'%s %f',1,'headerlines',54);

[KE_name1,KE_name2,KE_step]= textread(name,'%s %s %f',1,'headerlines',11);
[KE_max_name1,KE_max_name2,KE_max]= textread(name,'%s %s %f',1,'headerlines',12);
[KE_min_name1,KE_min_name2,KE_min]= textread(name,'%s %s %f',1,'headerlines',10);

detector_angle=angle_max-angle_min;
number_of_pixels=int32((angle_max-angle_min)/step)+1;
number_of_energies=int16((KE_max-KE_min)/KE_step);

 for i=1:number_of_energies
         KE(i)=KE_min+(double(i)-1)*KE_step;
 end;

for i=1:number_of_files
    file_name = sprintf('SPE%05d_%05d.txt',file,list_of_files(i,1));           
    data(:,:,i) = dlmread(file_name,'\t',57,0);
%     timestamp = textread(file_name,'%s',5,'headerlines',79);
%     time_start(i) = timestamp(3); 
    for j=1:(number_of_pixels)
        kinetic_energy(:,j,i) = KE(1,:);
        binding_energy(:,j,i) = EF - kinetic_energy(:,j,i);
    end;    display(i);
end;

n_e_angle = list_of_files(center_x_pixel,2) 

siz = size(data);
angle = (list_of_files(:,2)-n_e_angle).*(2*pi)/360; 

angle_step= list_of_files(1,3)  % displays the angle step size

KE_step                         % displays the energy step size



% --------- Normalization: ------------------------------------------------
for i=1:siz(3)
    for j=1:siz(2)
        data_norm_linewise(:,j,i) = smooth(data(:,j,i)./max(data(:,j,i)),5);      %line- / energywise normalization
    end;
    for k=1:siz(1)
        data_norm_columnwise(k,:,i) = smooth(data(k,:,i)./max(data(k,:,i)),5);
    end;
end;

%%%%-------choose which normalization here:--------------------------------
data = data;                                   % no normalization
% data = data_norm_linewise;                   % vertical lines in spectra
% data = data_norm_columnwise;                 % horizontal lines in spectra

%ANGULAR / MOMENTUM SCALE X - deflextion angles
for i = 1:siz(2)
    for j = 1:siz(3)
        k_parallel_x(:,i,j) = 0.5123*sqrt(kinetic_energy(:,1))*sin(angle(j,1));
    end;
end;    
%ANGULAR / MOMENTUM SCALE Y - Detector given angles along the slit
for i = 1:siz(2)
    for j = 1:siz(3)
        k_parallel_y(:,i,j) = 0.5123*sqrt(kinetic_energy(:,1))*sin(((i - center_y_pixel)*detector_angle/siz(2))*2*pi/360);
    end;    
end;





% ------------------------PLOTTING-----------------------------------------


%% -----FIGURE with all spectra taken during the entire scan:--------------

figure(1); clf; hold on;
for i=1:siz(3)
    subplot(ceil(sqrt(siz(3))),ceil(sqrt(siz(3))),i)
        imagesc(data(:,:,i));colormap(colorm);
        title(list_of_files(i,2));
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); set(gca,'ydir','normal'); set(gca, 'CLim', [lower*min(min(data(:,:,i))) upper*max(data(:))]); %INTENSITY
%%% Falls du die einzelnen Bilder als Datenstapel speichern möchtest (um sie
%%% zum Beispiel mit ImageJ anzuschauen... einfach die if-Abfrage
%%% auskommentieren...
%     if i == 1
%         imwrite(uint16(data(:,:,i)),'kxspectra.tif','Compression','None');
%     else
%         imwrite(uint16(data(:,:,i)),'kxspectra.tif','Compression','None','WriteMode','append');
%     end;  
end;


%% --------------certain spectra along kx direction:-----------------------

figure(11); hold on;
pixel = center_y_pixel;
for p = 1:9  
    for i = 1:siz(3)
        spec_x(:,i,p) = mean(data(:,pixel-3:pixel+3,i),2);
        spec_smooth(:,i,p) = smooth(spec_x(:,i,p),5);
    end;    
subplot(3,3,p)
    surf(squeeze(k_parallel_x(:,pixel,:)),squeeze(binding_energy(:,pixel,:)),squeeze(binding_energy(:,pixel,:))*0,...
        real2rgb(spec_smooth(:,:,p),colorm,[lower*min(min(spec_smooth(:,:,p))) upper*max(spec_smooth(:))]));            %INTENSITY
    view(2); shading interp; set(gca,'ydir','reverse');
    em_angle=center_y_pixel - pixel; txt = [num2str(em_angle),' pixel off-normal']; title(txt);
    
    
    
    pixel = pixel + 10;    %CHANGE

    
    
    
    ylim([min(min(min(binding_energy))) max(max(max(binding_energy)))]); 
    xlim([min(min(min(k_parallel_x))) max(max(max(k_parallel_x)))]); 
end;

%% -----------------constant energy cuts ----------------------------------


start_pix = 179;                 %defines the cut at certain energy. averages from start_pix to end_pix 
end_pix = 181; 


for kk = 1:12
    for i=1:siz(3)
        norm_start=end_pix;   norm_end=siz(1);  
        Fermi_surf(kk,:,i) = sum(data(start_pix:end_pix,:,i));
        %Fermi_surf(kk,:,i) = sum(data(start_pix:end_pix,:,i))/sum(sum(data(norm_start:norm_end,:,i)));
        Fermi_surf_smooth(kk,:,i) = smooth(Fermi_surf(kk,:,i),5);
    end;
        
        pix=(start_pix+end_pix)/2;
        
        E_B = round(binding_energy(pix)*1000);
        str = ['E_B = ', num2str(E_B),' meV'];  
        
        
        
% ------------------ alle in einen 3D Plot:--------------------------------            
    figure(111);
    if mod(kk,2) == 0
        surf(squeeze(k_parallel_x(pix,:,:)),squeeze(k_parallel_y(pix,:,:)),squeeze(binding_energy(pix,:,:)),...
            real2rgb(squeeze(Fermi_surf_smooth(kk,:,:)),'temp',[0 upper*max(max(squeeze(Fermi_surf_smooth(kk,:,:))))])); %INTENSITY normalizes every image to its particular maximum intensity
        set(gca,'zdir','reverse'); set(gca,'xdir','reverse'); set(gca,'ydir','reverse');
        shading flat; grid on; view(45,20); 
        axis equal; 
        xlabel('k_{|| x} [1/Å]'); ylabel('k_{|| y} [1/Å]'); zlabel('binding energy [eV]');
        daspect([1 1 0.3]); 
        hold on;
        view(200,30); 
    end;
%     surf(squeeze(k_parallel_x(start_pix,:,:)),squeeze(k_parallel_y(start_pix,:,:)),squeeze(binding_energy(start_pix,:,:)),real2rgb(squeeze(Fermi_surf(1,:,:)),'temp',[0 10000]))
%     
% ------------------ all in one Subplot:-----------------------------------            
    figure(222); 
    subplot(3,4,kk); hold on; title(str);  
        surf(squeeze(k_parallel_x(start_pix,:,:)),squeeze(k_parallel_y(start_pix,:,:)),squeeze(k_parallel_x(start_pix,:,:)*0),...
            real2rgb(squeeze(Fermi_surf(kk,:,:)),colorm,[0 upper*max(max(squeeze(Fermi_surf(kk,:,:))))]));  %INTENSITY normalisiert jedes Bild aufs jeweilige Maximum
        view(2); shading interp; axis equal;
%         ylim([-0.13 0.13]); 
%         xlim([-0.125 0.125]);                                       
                                 
    start_pix = start_pix - 10;   end_pix = end_pix - 10;             % defines energy step between two cuts (number*energy_step)
end;


%% ------------------single constant energy cut ---------------------------



pixels = 1;             %CHANGE




if pixels == 1
% hier in pixelskala: 
figure(11111);                                  
    kk = 2;                                     % choose subplot from figure 111 !
    imagesc(squeeze(Fermi_surf(kk,:,:)));
    

else    
%-------------- with momentum scale:---------------------------------------

figure(11111);
    kk = 2;                                     %choose subplot from figure 111 !
            surf(squeeze(k_parallel_x(kk,:,:)),squeeze(k_parallel_y(kk,:,:)),squeeze(k_parallel_x(kk,:,:)*0),...
            real2rgb(squeeze(Fermi_surf(kk,:,:)),colorm,[lower*min(min(squeeze(Fermi_surf(kk,:,:)))) upper*max(max(squeeze(Fermi_surf(kk,:,:))))]));    %Intensity
    view(2); shading interp; 
    axis equal;
    grid off; box on; 
%     ylim([-0.09 0.09]);
%     xlim([-0.09 0.09]);
    set(gca,'FontSize',18);
    ylabel('k_{|| y} [1/Å]','Fontsize',20); xlabel('k_{|| x} [1/Å]','Fontsize',20);
end;

%% ------------------ 3D Volume plot --------------------------------------
figure(3); clf; hold on; box on;
iso = -1;


% 1.half
ymin=1; ymax=450;                                    % ky, along the slit
xmin=1; xmax=center_x_pixel;                         % kx, perpendicular to the slit
zmin=1; zmax=170;                                    % energy [(start_pix+end_pix)/2 MINUS ONE STEP]
D = permute(data,[2 3 1]);
[x,y,z,D] = subvolume(D,[xmin,xmax,ymin,ymax,zmin,zmax]);
patch(isocaps(permute(k_parallel_x(zmin:zmax,ymin:ymax,xmin:xmax),[2 3 1]),permute(k_parallel_y(zmin:zmax,ymin:ymax,xmin:xmax),[2 3 1]),permute(binding_energy(zmin:zmax,ymin:ymax,xmin:xmax),[2 3 1]),D,iso),'FaceColor', 'interp', 'EdgeColor', 'none');

% + one more quarter:
ymin2=center_y_pixel; ymax2=ymax;                   % ky
xmin2=xmax; xmax2=siz(3);                           % kx 
zmin2=1; zmax2=zmax;                                  % energy
D = permute(data,[2 3 1]);
[x,y,z,D] = subvolume(D,[xmin2,xmax2,ymin2,ymax2,zmin2,zmax2]);
patch(isocaps(permute(k_parallel_x(zmin2:zmax2,ymin2:ymax2,xmin2:xmax2),[2 3 1]),permute(k_parallel_y(zmin2:zmax2,ymin2:ymax2,xmin2:xmax2),[2 3 1]),permute(binding_energy(zmin2:zmax2,ymin2:ymax2,xmin2:xmax2),[2 3 1]),D,iso),'FaceColor', 'interp', 'EdgeColor', 'none');


set(gca,'zdir','reverse');
view(45,25); 
axis equal; daspect([1.5,1.5,1])
colormap(temp); caxis([lower*0 upper*max(data(:))]);
set(gca,'FontSize',16,'FontWeight','bold');
xlabel('k_{x} [1/Å]','FontSize',16,'FontWeight','bold')
ylabel('k_{y} [1/Å]','FontSize',16,'FontWeight','bold')
zlabel('Binding energy [eV]','FontSize',16,'FontWeight','bold')

%-----Lines to highlight the cuts through the 3-dimensional space ---------
xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
zLimits = get(gca,'ZLim');

A=[xLimits(1) k_parallel_x(1,1,xmax) k_parallel_x(1,1,xmax) k_parallel_x(1,1,xmax) k_parallel_x(1,1,xmax) k_parallel_x(1,1,xmax) k_parallel_x(1,1,xmax) xLimits(2) xLimits(2) xLimits(2) xLimits(2)];
B=[yLimits(1) yLimits(1) yLimits(1) yLimits(1) k_parallel_y(1,ymin2,1)  k_parallel_y(1,ymin2,1) k_parallel_y(1,ymin2,1) k_parallel_y(1,ymin2,1) k_parallel_y(1,ymin2,1) k_parallel_y(1,ymin2,1) yLimits(2)];
C=[zLimits(1) zLimits(1) zLimits(2) zLimits(1) zLimits(1) zLimits(2) zLimits(1) zLimits(1) zLimits(2) zLimits(1) zLimits(1)];

plot3(A,B,C,'k-')





















