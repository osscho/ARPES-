
close all;
clear all; 




% -------------files: --------------------------------------------------------
file = 233;         %CHANGE
region=76;          %CHANGE

name = sprintf('SPE%05d_%05d.txt',file,region);        
data = dlmread(name,'\t',57,0);
siz = size(data);   

%--------------------------------------------------------------------------

EF = 1.64;                      %CHANGE

center_y_pixel = 276;           %CHANGE 


%------- Normalization ----------------------------------------------------
for j=2:siz(2)
    data_norm_linewise(:,j) = smooth(data(:,j)./max(data(:,j)),5);      %line- / energywise normalization
end;
for i = 1:siz(1)
    data_norm_columnwise(i,:) = smooth(data(i,:)/max(data(i,:)),5);
end;

data_norm = data_norm_linewise;

for i = 1:siz(2)
    kinetic_energy(:,i) = data(:,1);
    binding_energy(:,i) = EF - kinetic_energy(:,i);
end;

%ANGULAR / MOMENTUM SCALE Y - Detector given range ------------------------
[step_name, step]= textread(name,'%s %f',1,'headerlines',52);
[angle_max_name, angle_max]= textread(name,'%s %f',1,'headerlines',53);
[angle_min_name, angle_min]= textread(name,'%s %f',1,'headerlines',54);

detector_angle=angle_max-angle_min;

for i = 1:siz(2)
    k_parallel_y(:,i) = 0.5123*sqrt(kinetic_energy(:,1))*sin(((i - center_y_pixel)*detector_angle/siz(2))*2*pi/360);   %CHANGE number of pixels to siz(2)
end;
    k_parallel_y_1D = k_parallel_y(siz(1),:);
    delta_y = abs(k_parallel_y(2,1)-k_parallel_y(2,2));
    delta_E = abs(binding_energy(1,1)-binding_energy(2,1));

    
%-------calculation of EDC,MDC, 2ndDerivs and Curvatures-------------------

% erste_diff = zeros(siz(1),siz(2));
% zweite_diff = zeros(siz(1),siz(2));
curvature_MDC_1 = zeros(siz(1),siz(2));
curvature_EDC_1 = zeros(siz(1),siz(2));

strength = 4;
int = 3;
for i = 1:siz(1)
    if i+int<siz(1)     MDC_1(i,:) = mean(data(i:i+int,:),1);
    else                MDC_1(i,:) = mean(data(i-int:i,:),1);  end;

    MDC_norm_1(i,:) = smooth(MDC_1(i,:) ./ max(MDC_1(i,:)),20);
    deriv_MDC_1(i,:) = del2(MDC_norm_1(i,:),delta_y);
    deriv_MDC_norm_1(i,:) = smooth(-1 .* deriv_MDC_1(i,:) ./ max(deriv_MDC_1(i,:)),strength);
    C_1(i,:) = 0.01 * max(abs(gradient(MDC_norm_1(i,:),delta_y))).^2;
    curvature_MDC_1(i,:) = del2(smooth(MDC_norm_1(i,:),strength),delta_y)./(sqrt(C_1(i,:) + gradient(smooth(MDC_norm_1(i,:),strength),delta_y).^2)).^3;
    curvature_MDC_norm_1(i,:) = smooth(-1 .* curvature_MDC_1(i,:) ./ max(curvature_MDC_1(i,:)),strength);
end;
MDC_int = sum(MDC_1(:,:),1);

%------------weiter----- EDC s---------------------------------------------
for j = 2:siz(2)  
    if j+int<siz(2)     EDC_1(:,j) = mean(data(:,j:j+int),2);
    else                EDC_1(:,j) = mean(data(:,j-int:j),2);  end;
    
    EDC_norm_1(:,j) = smooth(EDC_1(:,j) ./ max(EDC_1(:,j)),20);
    deriv_EDC_1(:,j) = del2(EDC_norm_1(:,j),delta_E);
    deriv_EDC_norm_1(:,j) = smooth(-1 .* deriv_EDC_1(:,j) ./ max(deriv_EDC_1(:,j)),strength);
    C2_1(:,j) = 0.01 * max(abs(gradient(EDC_norm_1(:,j),delta_E))).^2;
    curvature_EDC_1(:,j) = del2(smooth(EDC_norm_1(:,j),strength),delta_E) ./ (sqrt(C2_1(:,j) + gradient(smooth(EDC_norm_1(:,j),strength),delta_E).^2)).^3;
    curvature_EDC_norm_1(:,j) = smooth(-1 .* curvature_EDC_1(:,j) ./ max(curvature_EDC_1(:,j)),strength);
end;
EDC_int = sum(EDC_1(:,:),2);
    



%% ----------------plotting -----------------------------------------------

% angle or k-scale: if angles = 1 --> angles;    
angles = 0;                                                 %CHANGE

lower = 1; upper = 0.85;                                     %CHANGE
%colorm = flipud(gray);
% colorm = 'hot';
colorm = 'temp';


figure(1); clf;  

    if angles == 1
        imagesc([1:siz(2)],binding_energy(:,1),data);
        colormap(colorm); caxis([lower*min(min(data)) upper*max(max(data))]);
        xlabel('detector angle [a.u.]','FontSize',15);
        
%         xlim([1 siz(2)]); colorbar; 
    
    else 
        data_surf=real2rgb(data,colorm,[lower*min(min(data)) upper*max(max(data))]); 
        surf(k_parallel_y,binding_energy,0.*data,data_surf);
        view(2);
        shading interp;
        colormap(colorm);
       
%         xlim([k_parallel_y_1D(1,1) k_parallel_y_1D(1,end)]);
        xlim([-0.15 0.15]);                                    %CHOSE / CHANGE
        
        set(gca,'FontSize',25);
        xlabel('k_{||} [1/Å]','FontSize',25);
        s=colorbar('Ticks',[0,1],'TickLabels',{'Low','High'});
        s.Label.String='Photoemission intensity';
        
    end;
        grid off; shading interp; box on;
        ylabel('Binding energy [eV]','FontSize',25);
        
%         ylim([binding_energy(end,1) binding_energy(1,1)]); 
        ylim([-0.1 0.4]);                                     %CHOSE / CHANGE
        
        set(gca,'ydir','reverse');
        pbaspect([2 3 1]);                                      %CHANGE

%% ------------------------------------------------------------------------

% saveas(gcf,sprintf('SPE%05d_%05d',file,region),'png');


  

