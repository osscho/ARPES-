% close all;
clear all; 

%Ep_Xe = 8.44;  %Ep_Ne = 16.848; %Ep_HeI = 21.22; %Ep_HeII = 40.8;

E_photon = 21.22;
EF = 4.43;

% -------------files: --------------------------------------------------------

file_1 = sprintf('SPE3428_11.txt');                



data = dlmread(file_1,'\t',82,0);
siz = size(data);     

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
    binding_energy(:,i) = E_photon - EF - kinetic_energy(:,i);
end;

center_x_pixel = 139;                                               %CHANGE !                                            

%ANGULAR / MOMENTUM SCALE X - Detector given range
detector_angle = 28;
for i = 1:siz(2)
    k_parallel_x(:,i) = 0.5123*sqrt(kinetic_energy(:,1))*sin(((i - center_x_pixel)*detector_angle/256)*2*pi/360);
end;
    k_parallel_x_1D = k_parallel_x(siz(1),:);
    delta_x = abs(k_parallel_x(2,1)-k_parallel_x(2,2));
    delta_E = abs(binding_energy(1,1)-binding_energy(2,1));

       
%-------calculation of EDC,MDC, 2ndDerivs and Curvatures---------------
% erste_diff = zeros(siz(1),siz(2));
% zweite_diff = zeros(siz(1),siz(2));
curvature_MDC_1 = zeros(siz(1),siz(2));
curvature_EDC_1 = zeros(siz(1),siz(2));

strength = 5;
int = 2;

for i = 1:siz(1)
    if i+int<siz(1)     MDC_1(i,:) = mean(data(i:i+int,:),1);
    else                MDC_1(i,:) = mean(data(i-int:i,:),1);  end;

    MDC_norm_1(i,:) = smooth(MDC_1(i,:) ./ max(MDC_1(i,:)),15);
    deriv_MDC_1(i,:) = del2(MDC_norm_1(i,:),delta_x);
    deriv_MDC_norm_1(i,:) = smooth(-1 .* deriv_MDC_1(i,:) ./ max(deriv_MDC_1(i,:)),strength);
    C_1(i,:) = 0.01 * max(abs(gradient(MDC_norm_1(i,:),delta_x))).^2;
    curvature_MDC_1(i,:) = del2(smooth(MDC_norm_1(i,:),strength),delta_x)./(sqrt(C_1(i,:) + gradient(smooth(MDC_norm_1(i,:),strength),delta_x).^2)).^3;
    curvature_MDC_norm_1(i,:) = smooth(-1 .* curvature_MDC_1(i,:) ./ max(curvature_MDC_1(i,:)),strength);
end;
MDC_int = sum(MDC_1(:,:),1);

%------------weiter----- EDC s---------------------------------------------
for j = 2:siz(2)  
    if j+int<siz(2)     EDC_1(:,j) = mean(data(:,j:j+int),2);
    else                EDC_1(:,j) = mean(data(:,j-int:j),2);  end;
    
    EDC_norm_1(:,j) = smooth(EDC_1(:,j) ./ max(EDC_1(:,j)),15);
    deriv_EDC_1(:,j) = del2(EDC_norm_1(:,j),delta_E);
    deriv_EDC_norm_1(:,j) = smooth(-1 .* deriv_EDC_1(:,j) ./ max(deriv_EDC_1(:,j)),strength);
    C2_1(:,j) = 0.01 * max(abs(gradient(EDC_norm_1(:,j),delta_E))).^2;
    curvature_EDC_1(:,j) = del2(smooth(EDC_norm_1(:,j),strength),delta_E) ./ (sqrt(C2_1(:,j) + gradient(smooth(EDC_norm_1(:,j),strength),delta_E).^2)).^3;
    curvature_EDC_norm_1(:,j) = smooth(-1 .* curvature_EDC_1(:,j) ./ max(curvature_EDC_1(:,j)),strength);
end;
EDC_int = sum(EDC_1(:,:),2);

%______________plotting____________________________________________________
xscale = 0.45;


%%%energy range:
en_start=1; en_end=siz(1)

figure(11); clf; fac = 1;
subplot(2,4,[1 5]);
     imagesc(k_parallel_x_1D,binding_energy(en_start:en_end,1),data(en_start:en_end,:));
     colormap(gray); caxis([100 5000]);
     set(gca,'ydir','reverse'); grid off; shading interp; box on;
     xlabel('k_{parallel_x} [1/Å]'); ylabel('Binding energy [eV]');
     xlim([-0.35 0.35]); % ylim([-0.05 0.4]);

subplot(2,4,[2 6]);
     imagesc(k_parallel_x_1D,binding_energy(en_start:en_end,1),data_norm(en_start:en_end,:));
     colormap(temp);
     set(gca,'ydir','reverse'); grid off; shading interp; box on;
     xlabel('k_{parallel_x} [1/Å]'); 
     xlim([-0.35 0.35]); % ylim([-0.05 0.4]);
    
subplot(2,4,[3 7]);
    imagesc(k_parallel_x_1D,binding_energy(en_start:en_end,1),curvature_MDC_norm_1(en_start:en_end,:)*fac);
    colormap(temp); caxis([0 0.9]);
    title('Curvature from MDCs');
    set(gca,'ydir','reverse'); xlabel('k_{parallel_x} [1/Å]');
    xlim([-0.35 0.35]);

subplot(2,4,[4 8]);
    imagesc(k_parallel_x_1D,binding_energy(en_start:en_end,1),curvature_EDC_norm_1(en_start:en_end,:)*fac);
    colormap(temp); caxis([0 0.9]);
    title('Curvature from EDCs');
    set(gca,'ydir','reverse'); xlabel('k_{parallel_x} [1/Å]');
    xlim([-0.35 0.35]);

   
lower = 1;
upper = 0.3;
%colorm = flipud(gray);
colorm = 'hot';

figure(111)
clf; fac = 1;
%subplot(1,2,1);
     surf(k_parallel_x,binding_energy,binding_energy*0,real2rgb(data,colorm,[lower*min(min(data)) upper*max(max(data))]));
     view(2), shading interp;
     set(gca,'ydir','reverse'); grid off; shading interp; box on;
     xlabel('k_{parallel_x} [1/Å]','FontSize',15); ylabel('Binding energy [eV]','FontSize',15);
     xlim([-0.45 0.45]); % ylim([-0.05 0.4]);
     set(gca,'FontSize',15);
% subplot(1,2,2);
%     plot(EDC_int,binding_energy(en_start:en_end,1),'LineWidth',3);
%     set(gca,'ydir','reverse');
%     ylabel('Binding energy [eV]','FontSize',15);
%     xlabel('intensity','FontSize',15);
%     set(gca,'FontSize',15);
   
        
        
  


    
    





