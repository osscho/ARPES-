clear all;
%close all;

%sherman = input('Sherman: ');
sherman = 0.2;
%normally around 0.2

%--------------------------------------------------------------------
% files: 

NE=933;

%Coil1pos = [];
%Coil1neg = [];
Coil2pos = [9638,9641,9642,9645,9646,9649,9650,9653,9654,9657,9658,9661,9662,9665,9666,9669,9670,9673,9674,9677,9678,9681,9682]';
Coil2neg = [9639,9640,9643,9644,9647,9648,9651,9652,9655,9656,9659,9660,9663,9664,9667,9668,9671,9672,9675,9676,9679,9680,9683]';

angles = [927:0.5:938]'-NE;

EF = 15.9;

for i = 1:size(Coil2pos,1)
%    file_1pos = sprintf('SPE%d_07.mot',Coil1pos(1,i));   
%    file_1neg = sprintf('SPE%d_07.mot',Coil1neg(1,i)); 
    file_2pos = sprintf('SPE%d_07.mot',Coil2pos(i,1));   
    file_2neg = sprintf('SPE%d_07.mot',Coil2neg(i,1)); 
    
%    in(:,:,i) = dlmread(file_1pos,'\t',1,0);
%    out(:,:,i) = dlmread(file_1neg,'\t',1,0);
    left(:,:,i) = dlmread(file_2pos,'\t',1,0);
    right(:,:,i) = dlmread(file_2neg,'\t',1,0);
  
energy(:,i) = -(left(:,1) - EF);
end;

%in = squeeze(in(:,2,:));
%out = squeeze(out(:,2,:));
left = squeeze(left(:,2,:));
right = squeeze(right(:,2,:));
siz = size(left); 


for j = 1:siz(2)
    kvector(:,j) = 0.5123*sqrt(-energy(:,j)+EF)*sin(angles(j,1)*pi/180)';
end;

angles = kvector(:,:);

strenght = 5;
% in_smooth = smooth(in,strenght);
% out_smooth = smooth(out,strenght);
left_smooth = smooth(left,strenght);
right_smooth = smooth(right,strenght);

%---------------------------------------------------------------------
% %   normalize   normalize channeltron 2+3 and 4+5 to the biggest number
for i = 1:siz(2)       
%    in_norm(:,i) = in(:,i)./sum(in(:,i));
%    out_norm(:,i) = out(:,i)./sum(out(:,i));
    left_norm(:,i) = left(:,i)./sum(left(:,i),1);
    right_norm(:,i) = right(:,i)./sum(right(:,i),1);
end;

%left = left_norm;
%right = right_norm;


%----------spin polarization calculation------------
   % SPLR = 1/sherman * ( (sqrt(left_up.*right_down) - sqrt(left_down.*right_up))./ (sqrt(left_up.*right_down) + sqrt(left_down.*right_up)));
   DiffLR = left - right; 
   SPLR = 1/sherman.*(left - right)./(left + right);
       
%SPIO=1/sherman*(in - out)./(in + out);    
%SPIO = 1/sherman * ( (sqrt(in_up.*out_down) - sqrt(in_down.*out_up))./ (sqrt(in_up.*out_down) + sqrt(in_down.*out_up)));



delta_x = abs(kvector(2,1)-kvector(1,1));
delta_E = abs(energy(1,1)-energy(2,1));

strength = 4
 int = 1;
%--------------SP - EDC s---------------------------------------------
for j = 2:siz(2)  
    if j+int<siz(2)     EDC_SPLR(:,j) = mean(SPLR(:,j:j+int),2);
    else                EDC_SPLR(:,j) = mean(SPLR(:,j-int:j),2);  end;
    
    EDC_SPLR_norm(:,j) = smooth(EDC_SPLR(:,j) ./ max(EDC_SPLR(:,j)),20);
    deriv_SPLR(:,j) = del2(EDC_SPLR_norm(:,j),delta_E);
    deriv_SPLR_norm(:,j) = smooth(-1 .* deriv_SPLR(:,j) ./ max(deriv_SPLR(:,j)),strength);
    C2_SPLR(:,j) = 0.01 * max(abs(gradient(EDC_SPLR_norm(:,j),delta_E))).^2;
    curvature_SPLR(:,j) = del2(smooth(EDC_SPLR_norm(:,j),strength),delta_E) ./ (sqrt(C2_SPLR(:,j) + gradient(smooth(EDC_SPLR_norm(:,j),strength),delta_E).^2)).^3;
    curvature_SPLR_norm(:,j) = smooth(-1 .* curvature_SPLR(:,j) ./ max(curvature_SPLR(:,j)),strength);
end;
EDC_SPLR_int = sum(EDC_SPLR(:,:),2);
for j = 2:siz(2)  
    if j+int<siz(2)     EDC_DiffLR(:,j) = mean(DiffLR(:,j:j+int),2);
    else                EDC_DiffLR(:,j) = mean(DiffLR(:,j-int:j),2);  end;
    
    EDC_DiffLR_norm(:,j) = smooth(EDC_DiffLR(:,j) ./ max(EDC_DiffLR(:,j)),20);
    deriv_DiffLR(:,j) = del2(EDC_DiffLR_norm(:,j),delta_E);
    deriv_DiffLR_norm(:,j) = smooth(-1 .* deriv_DiffLR(:,j) ./ max(deriv_DiffLR(:,j)),strength);
    C2_DiffLR(:,j) = 0.01 * max(abs(gradient(EDC_DiffLR_norm(:,j),delta_E))).^2;
    curvature_DiffLR(:,j) = del2(smooth(EDC_DiffLR_norm(:,j),strength),delta_E) ./ (sqrt(C2_DiffLR(:,j) + gradient(smooth(EDC_DiffLR_norm(:,j),strength),delta_E).^2)).^3;
    curvature_DiffLR_norm(:,j) = smooth(-1 .* curvature_DiffLR(:,j) ./ max(curvature_DiffLR(:,j)),strength);
end;
EDC_DiffLR_int = sum(EDC_DiffLR(:,:),2);
                       

SPLR_smooth = smooth(SPLR,strenght);
%SPIO_smooth = smooth(SPIO,strenght);

sum_left_right = left + right;
sum_left_right_smooth = smooth(sum_left_right,strenght);
% sum_in_out = in + out;
% sum_in_out_smooth = smooth(sum_in_out,strenght);

%calculate back REAL Spectrum out of measured Asymmetry !!!
LEFT_spectrum = (sum_left_right/2) .* (1+SPLR);
LEFT_spectrum_smooth = smooth(LEFT_spectrum,strenght);
LEFT_spectrum_norm = LEFT_spectrum/max(LEFT_spectrum);
RIGHT_spectrum = (sum_left_right/2) .* (1-SPLR);
RIGHT_spectrum_smooth = smooth(RIGHT_spectrum,strenght);
RIGHT_spectrum_norm = RIGHT_spectrum/max(RIGHT_spectrum);
% IN_spectrum = (sum_in_out/2) .* (1+SPIO);
% IN_spectrum_smooth = smooth(IN_spectrum,strenght);
% IN_spectrum_norm = IN_spectrum/max(IN_spectrum);
% OUT_spectrum = (sum_in_out/2) .* (1-SPIO);
% OUT_spectrum_smooth = smooth(OUT_spectrum,strenght);
% OUT_spectrum_norm = OUT_spectrum/max(OUT_spectrum);



colorm = jet(512);
xscale = 0.3;
Emin = -0.3; Emax = 1.5;

figure(1); clf;
subplot(2,4,1)
    surf(angles,energy,left);
    set(gca,'ydir','reverse','FontSize',14);  xlim([-xscale xscale]); ylim([Emin Emax]);
    view(2); shading interp; box on;
    title('raw data: LEFT');
   
subplot(2,4,2)
    surf(angles,energy,right);
     set(gca,'ydir','reverse','FontSize',14);  xlim([-xscale xscale]); ylim([Emin Emax]);
    view(2); shading interp; box on;
    title('raw data: RIGHT');
subplot(2,4,3)
    surf(angles,energy,DiffLR);
     set(gca,'ydir','reverse','FontSize',14);  xlim([-xscale xscale]); ylim([Emin Emax]);
    view(2); shading interp; box on;
    title('Difference'); 
subplot(2,4,4)
    surf(angles,energy,curvature_DiffLR_norm);
    set(gca,'ydir','reverse','FontSize',14);  xlim([-xscale xscale]); ylim([Emin Emax]);
    view(2); shading interp; box on;
    title('2nd derivative of Diff');
    
subplot(2,4,5)
    surf(angles,energy,LEFT_spectrum);
     set(gca,'ydir','reverse','FontSize',14);  xlim([-xscale xscale]); ylim([Emin Emax]);
    view(2); shading interp; box on;
    title('recalculated LEFT');
subplot(2,4,6)
    surf(angles,energy,RIGHT_spectrum);
     set(gca,'ydir','reverse','FontSize',14);  xlim([-xscale xscale]); ylim([Emin Emax]);
    view(2); shading interp; box on;
    title('recalculated RIGHT');
subplot(2,4,7)
    surf(angles,energy,SPLR);
     set(gca,'ydir','reverse','FontSize',14);  xlim([-xscale xscale]); ylim([Emin Emax]);
    view(2); shading interp; box on;
    title('Spin Polarization'); 
subplot(2,4,8)
    surf(angles,energy,curvature_SPLR_norm);
    set(gca,'ydir','reverse','FontSize',14);  xlim([-xscale xscale]); ylim([Emin Emax]);
    view(2); shading interp; box on;
    title('2nd derivative of spin pol'); 
    
    colormap(colorm);
    
    
create_temperature_colormap_512;

    
figure(56)

subplot(1,6,1)
surf(angles,energy,log(left+right),'edgecolor','none')
view(2); axis('tight');
set(gca,'ydir','reverse');
colorbar;
grid off
box on
title('log(left+right)')

subplot(1,6,2)
colormap(colormap_temp);
surf(angles,energy,SPLR,'edgecolor','none')
shading flat;
view(2); axis('tight'); 
alpha(log(left+right))
grid off
box on
set(gca,'ydir','reverse');
colorbar
caxis([-0.5 0.5])
title('SPRL w log(l+r)')

subplot(1,6,4)
colormap(colormap_temp);
surf(angles,energy,SPLR,'edgecolor','none')
shading flat;
view(2); axis('tight'); 
intensity_matrix = left+right;
intensity_matrix = Saturate_matrix(intensity_matrix,0.3);
alpha(intensity_matrix)
grid off
box on
set(gca,'ydir','reverse');
colorbar
caxis([-0.5 0.5])
title('SPRL w/saturated alpha')

subplot(1,6,3)
surf(angles,energy,intensity_matrix,'edgecolor','none')
view(2); axis('tight');
set(gca,'ydir','reverse');
colorbar;
grid off
box on
title('Left+right saturated')

subplot(1,6,6)
colormap(colormap_temp);
surf(angles,energy,SPLR,'edgecolor','none')
shading flat;
view(2); axis('tight'); 
intensity_matrix = left+right;
intensity_matrix = DOS_normalize_matrix(intensity_matrix);
alpha(intensity_matrix)
grid off
box on
set(gca,'ydir','reverse');
colorbar
caxis([-0.5 0.5])
title('SPRL DOS norm')

subplot(1,6,5)
surf(angles,energy,intensity_matrix,'edgecolor','none')
view(2); axis('tight');
set(gca,'ydir','reverse');
colorbar;
grid off
box on
title('Left+right DOS norm')


%shading interp;

colormap(colormap_temp);
% 
% figure(1);
%     colormap(gray(256));

    
    
    
    
    
    
% 1D plots: 

figure(2); clf; 
imagesc(SPLR), set(gca,'ydir','normal');

cut = 53;          %fermi=112;
%cut2 = 93;
angles = angles(cut,:);

E_int = 4;
left_rim = 4; 
right_rim = 22;
angle_int = 1;

SP_MDC = mean(SPLR(cut-E_int:cut+E_int,:,1));
SP_EDC_left = mean(SPLR(:,left_rim-angle_int:left_rim+angle_int),2);
SP_EDC_right = mean(SPLR(:,right_rim-angle_int:right_rim+angle_int),2);

figure(4); clf; box on; hold on; grid on;
    plot(angles,SP_MDC,'o');
    plot(angles,smooth(SP_MDC,5),'-'); 
    set(gca,'FontSize',14);
    ylabel('Spin polarization [%]'); 
    xlabel('wave vector [1/A]');

figure(5); clf; 
subplot(1,2,1); hold on; box on; grid on;
    plot(SP_EDC_left,energy,'^r');
    plot(smooth(SP_EDC_left,5),energy,'r-');
    set(gca,'ydir','reverse','FontSize',14);
    xlabel('Spin polarization [%]'); 
    ylabel('Binding energy [eV]');
subplot(1,2,2); hold on; box on; grid on;
    plot(SP_EDC_right,energy,'vb');
    plot(smooth(SP_EDC_right,5),energy,'b-');
    set(gca,'ydir','reverse','FontSize',14);
    xlabel('Spin polarization [%]'); 
    ylabel('Binding energy [eV]');

    figure(56)