%   File for non-magnetic samples whose magnetization direction cannot be
%   manipulated and changed... large instrumental asymmetries !!!

%sherman = input('Sherman: ');
clear all;
close all;
sherman = 0.27;
%normally around 0.2

%--------------------------------------------------------------------
%   Load data files: up & down (bzw. left & right) means magnetization
%   direction
%Coilone_pos = dlmread('SPE8784_07.mot','\t',1,0);       % oop, out
%Coilone_neg = dlmread('SPE8785_07.mot','\t',1,0);       % oop, in
Coiltwo_pos = dlmread('SPE%d.mot','\t',1,0);       % ip, up
Coiltwo_neg = dlmread('SPE%d.mot','\t',1,0);       % ip, down
%Coilone_pos = dlmread('SPE8742_07.mot','\t',1,0);       % again ...

%out = Coilone_pos(:,2);
%in = Coilone_neg(:,2);
up = Coiltwo_pos(:,2);
down = Coiltwo_neg(:,2);

left = up;
right = down;

strenght = 5;
up_smooth = smooth(up,strenght);
down_smooth = smooth(down,strenght);
%in_smooth = smooth(in,strenght);
%out_smooth = smooth(out,strenght);

left_smooth = smooth(up,strenght);
right_smooth = smooth(down,strenght);

%Fermi Energy, anpassen damit Bindungsenergie stimmt !!!
EF = 15.8;

%---------------------------------------------------------------------
energy = -(Coiltwo_pos(:,1) - EF);

% %   normalize   normalize channeltron 2+3 and 4+5 to the biggest number
% normalize = input('Normalize ? [1-yes ; 0-no]');
% 
% if normalize == 1
%     file(:,2:3) = file(:,2:3)/sum(sum(file(:,2:3)));
%     file(:,4:5) = file(:,4:5)/sum(sum(file(:,4:5)));
% end;



% further files for calculations...: 
% channel 2 = in    ->red
% channel 3 = out   ->yellow
% channel 4 = left  ->green         /green and blue are exchanged in the Hardware
% channel 5 = right ->blue

% left = file(:,4);
% left_smooth = file_smooth(:,4);
% left_norm = left/max(left);
% right = file(:,5);
% right_smooth = file_smooth(:,5);
% right_norm = right/max(right);
% 
% in = file(:,2);
% in_smooth = file_smooth(:,2);
% in_norm = in/max(in);
% out = file(:,3);
% out_smooth = file_smooth(:,3);
% out_norm = out/max(out);

%----------spin polarization calculation------------
fact = 1;
   % SPLR = 1/sherman * ( (sqrt(left_up.*right_down) - sqrt(left_down.*right_up))./ (sqrt(left_up.*right_down) + sqrt(left_down.*right_up)));
  
   SPLR=1/sherman*(left - right)./(left + right);
%    for i = 1:size(left_up)
%       e_left = fact * sqrt(left);
%       e_right = fact * sqrt(right);
% %             Fehler auf den Z‰hler (von innen nach auﬂen: Fehler auf Produkt,
% %             Fehler auf Wurzel, Fehler auf Summe/Differenz...)
%               dev_zaehler(i) = devsum(sqrt(left_up(i)*right_down(i)),sqrt(left_down(i)*right_up(i)),devpow(left_up(i)*right_down(i),devprod(left_up(i),right_down(i), e_left_up(i) , e_right_down(i)),0.5),devpow(left_down(i)*right_up(i),devprod(left_down(i),right_up(i), e_left_down(i), e_right_up(i)),0.5));
%                         %devtop = devsum(devpow(lue*rde,devprod(lue,rde,sqrt(lue),sqrt(rde)),0.5),devpow(lde*rue,devprod(lde,rue,sqrt(lde),sqrt(rue)),0.5))              
%               devSPLR(i) = (1/sherman) * devdiv(sqrt(left_up(i)*right_down(i))-sqrt(left_down(i)*right_up(i)),sqrt(left_up(i)*right_down(i))+ sqrt(left_down(i)*right_up(i)),dev_zaehler(i),dev_zaehler(i));        
%   end;
            %     % aus Spaltenvektor mach Zeilenvektor:
            %     devSPLR = rot90(devSPLR,3); 
    
  %  topLR = (sqrt(left_up.*right_down) - sqrt(left_down.*right_up));
  %  bottomLR = (sqrt(left_up.*right_down) + sqrt(left_down.*right_up));
  %  error_topLR = sqrt(0.5*(sqrt((right_down.*e_left_up).^2+(left_up.*e_right_down).^2)./sqrt(left_up.*right_down))+0.5*(sqrt((right_up.*e_left_down).^2+(left_down.*e_right_up).^2)./sqrt(left_down.*right_up)));
  %  error_bottomLR = error_topLR;
  %  e_SPLR = 1/sherman .* abs(SPLR) .* sqrt((error_topLR./topLR).^2+(error_bottomLR./bottomLR).^2);
    
%spin polarization IN - OUT out of the 2 states:
%SPIO=1/sherman*(in - out)./(in + out);    
%SPIO = 1/sherman * ( (sqrt(in_up.*out_down) - sqrt(in_down.*out_up))./ (sqrt(in_up.*out_down) + sqrt(in_down.*out_up)));
 %   for i = 1:size(left_up)
  %     e_in_up = fact * sqrt(in_up(i));
   %    e_in_down = fact * sqrt(in_down(i));
    %   e_out_up = fact * sqrt(out_up(i));
     %  e_out_down = fact * sqrt(out_down(i));
               %Fehler auf den Z‰hler (von innen nach auﬂen: Fehler auf Produkt,
                %Fehler auf Wurzel, Fehler auf Summe/Differenz...)
                %        dev_zaehler = devsum(sqrt(e_in_up*e_out_down),sqrt(e_in_down*e_out_up),devpow(e_in_up*e_out_down,devprod(e_in_up,e_out_down,fact*sqrt(e_in_up),fact*sqrt(e_out_down)),0.5),devpow(e_in_down*e_out_up,devprod(e_in_down,e_out_up,fact*sqrt(e_in_down),fact*sqrt(e_out_up)),0.5));
                %             %devtop = devsum(devpow(lue*rde,devprod(lue,rde,sqrt(lue),sqrt(rde)),0.5),devpow(lde*rue,devprod(lde,rue,sqrt(lde),sqrt(rue)),0.5))              
                %        devSPIO(i) = (1/sherman) * devdiv(sqrt(e_in_up*e_out_down)-sqrt(e_in_down*e_out_up),sqrt(e_in_up*e_out_down)+ sqrt(e_in_down*e_out_up),dev_zaehler,dev_zaehler);
    %end
                %     % aus Spaltenvektor mach Zeilenvektor:
                %     devSPIO = rot90(devSPIO,3);  
    
   % topIO = (sqrt(in_up.*out_down) - sqrt(in_down.*out_up));
   % bottomIO = (sqrt(in_up.*out_down) + sqrt(in_down.*out_up));
   % error_topIO = sqrt(0.5*(sqrt((out_down.*e_in_up).^2+(in_up.*e_out_down).^2)./sqrt(in_up.*out_down))+0.5*(sqrt((out_up.*e_in_down).^2+(in_down.*e_out_up).^2)./sqrt(in_down.*out_up)));
   % error_bottomIO = error_topIO;
   % e_SPIO = 1/sherman .* abs(SPIO) .* sqrt((error_topIO./topIO).^2+(error_bottomIO./bottomIO).^2);
   

                             
%evtl SPLR und SPIO noch smoothen,
SPLR_smooth = smooth(SPLR,strenght);
%SPIO_smooth = smooth(SPIO,strenght);

sum_left_right = left + right;
sum_left_right_smooth = smooth(sum_left_right,strenght);
%sum_in_out = in + out;
%sum_in_out_smooth = smooth(sum_in_out,strenght);

%calculate back REAL Spectrum out of measured Asymmetry !!!
LEFT_spectrum = (sum_left_right/2) .* (1+SPLR);
LEFT_spectrum_smooth = smooth(LEFT_spectrum,strenght);
LEFT_spectrum_norm = LEFT_spectrum/max(LEFT_spectrum);
RIGHT_spectrum = (sum_left_right/2) .* (1-SPLR);
RIGHT_spectrum_smooth = smooth(RIGHT_spectrum,strenght);
RIGHT_spectrum_norm = RIGHT_spectrum/max(RIGHT_spectrum);
%{
IN_spectrum = (sum_in_out/2) .* (1+SPIO);
IN_spectrum_smooth = smooth(IN_spectrum,strenght);
IN_spectrum_norm = IN_spectrum/max(IN_spectrum);
OUT_spectrum = (sum_in_out/2) .* (1-SPIO);
OUT_spectrum_smooth = smooth(OUT_spectrum,strenght);
OUT_spectrum_norm = OUT_spectrum/max(OUT_spectrum);
%}

figure(22)
subplot(3,1,1);
box on;hold on; grid on;
title('measured Spectra: green: LEFT, blue: RIGHT, black: SUM');
% dots are measured data, line is given by the smoothed data
plot(energy,left,'g.')
plot(energy,left_smooth,'g','linewidth',2)
plot(energy,right,'b.')
plot(energy,right_smooth,'b','linewidth',2)
plot(energy,sum_left_right,'k.-');
plot(energy,sum_left_right_smooth,'k-','linewidth',2);
%Legende Farben zu Channels
%h = legend('in - Ch 1','out - Ch 2','left - Ch 3','right - Ch 4',2);
%set(h,'Interpreter','none')
%x-Achse falsch rum: 
set(gca,'xdir','reverse');
%axis tight;
%{
subplot(3,1,2);
box on; hold on; grid on;
title('measured spectra "IN-OUT"');
plot(energy,in,'r.')
plot(energy,in_smooth,'r','linewidth',2)
plot(energy,out,'y.')
plot(energy,out_smooth,'y','linewidth',2)
plot(energy,sum_in_out,'k.-');
plot(energy,sum_in_out_smooth,'k-','linewidth',2);
title(texlabel('measured Spectra: red: IN, yellow: OUT, black: SUM'));
%Legende: Farben zu Channels
%h = legend('in - Ch 1','out - Ch 2','left - Ch 3','right - Ch 4',2);
%set(h,'Interpreter','none')
set(gca,'xdir','reverse');
%axis tight;
%}
subplot(3,1,2);
hold on; box on; grid on;
title('recalculated Spin Spectra: green: LEFT (Majority?); blue: RIGHT (Minority)');
plot(energy,LEFT_spectrum,'g.-');
plot(energy,LEFT_spectrum_smooth,'g-','linewidth',3);
plot(energy,RIGHT_spectrum,'b.-');
plot(energy,RIGHT_spectrum_smooth,'b-','linewidth',3);
ylabel('Photoemission signal [counts]','fontsize',12);
set(gca,'xdir','reverse');
%axis tight;

%{
subplot(3,2,4);
hold on; box on; grid on;
title('recalculated Spin Spectra: blue: red (Majority?); yellow: OUT (Minority)');
plot(energy,IN_spectrum,'r.-');
plot(energy,IN_spectrum_smooth,'r-','linewidth',3);
plot(energy,OUT_spectrum,'y.-');
plot(energy,OUT_spectrum_smooth,'y-','linewidth',3);
ylabel('Photoemission signal [counts]','fontsize',12);
set(gca,'xdir','reverse');
%axis tight;
%}
%----wichtigster Plot: SPIN POLARIZATION (+- Std Dev)--------%


subplot(3,1,3); hold on; box on; grid on;
title('SPIN-POLARIZATION:  left-right (red)','fontsize',12);
plot(energy,SPLR,'r.-');
plot(energy,SPLR_smooth,'r-','linewidth',3);
%plot(energy,SPLR_smooth - e_SPLR,'g-','linewidth',3);
%plot(energy,SPLR_smooth + e_SPLR,'g-','linewidth',3);
xlabel('binding energy E_B [eV]','fontsize',12);
ylabel('Spin Polarization','fontsize',12);
set(gca,'xdir','reverse');


%{
subplot(3,2,6);
hold on;
box on;
grid on;
title('SPIN-POLARIZATION:  in-out (black)','fontsize',12);
plot(energy,SPIO,'b.-');
plot(energy,SPIO_smooth,'b-','linewidth',3);
%plot(energy,SPIO_smooth - e_SPIO,'b-','linewidth',3);
%plot(energy,SPIO_smooth + e_SPIO,'b-','linewidth',3);
xlabel('binding energy E_B [eV]','fontsize',12);
ylabel('Spin Polarization','fontsize',12);
set(gca,'xdir','reverse');
%}


















% %nochmal SPIN pol spectra ...
% figure(333)
% hold on;
% box on;
% title('SPIN-POLARIZATION left-right (red) and in-out (black)','fontsize',16);
% plot(energy,SPLR,'r.-');
% plot(energy,SPLR_smooth,'r-','linewidth',3);
% plot(energy,SPLR_smooth - devSPLR,'g-','linewidth',3);
% plot(energy,SPLR_smooth + devSPLR,'g-','linewidth',3);
% plot(energy,SPIO,'k.-');
% plot(energy,SPIO_smooth,'k-','linewidth',3);
% plot(energy,SPIO_smooth - devSPIO,'b-','linewidth',3);
% plot(energy,SPIO_smooth + devSPIO,'b-','linewidth',3);
% xlabel('binding energy E_B [eV]','fontsize',14);
% ylabel('Spin Polarization','fontsize',14);
% set(gca,'xdir','reverse');
% %axis tight;
% 
