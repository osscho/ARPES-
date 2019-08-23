

%sherman = input('Sherman: ');
clear all;
%close all;
sherman = 0.27
%normally around 0.2

%--------------------------------------------------------------------
%   Load data files: up & down (bzw. left & right) means magnetization direction

angles=[925:5:975];
angles_num=numel(angles);
start_file=3141;


file_number=start_file;

for i=1:angles_num
Coilone_neg_files(:,i)=file_number;
Coilone_pos_files(:,i)=file_number+1;
Coiltwo_neg_files(:,i)=file_number+2;
Coiltwo_pos_files(:,i)=file_number+3;
file_number=file_number+4;
end;






for i=1:angles_num
    f_1_n=Coilone_neg_files(i);
    f_1_p=Coilone_pos_files(i);
    f_2_n=Coiltwo_neg_files(i);
    f_2_p=Coiltwo_pos_files(i);
    
file_name_1_n = sprintf('SPE%d_07.mot',f_1_n);
file_name_1_p = sprintf('SPE%d_07.mot',f_1_p);
file_name_2_n = sprintf('SPE%d_07.mot',f_2_n);
file_name_2_p = sprintf('SPE%d_07.mot',f_2_p);

Coilone_neg_single = dlmread(file_name_1_n,'\t',2,0);
Coilone_pos_single = dlmread(file_name_1_p,'\t',2,0);
Coiltwo_neg_single = dlmread(file_name_2_n,'\t',2,0);
Coiltwo_pos_single = dlmread(file_name_2_p,'\t',2,0);

Coilone_neg(:,i)=Coilone_neg_single(:,2);
Coilone_pos(:,i)=Coilone_pos_single(:,2);
Coiltwo_neg(:,i)=Coiltwo_neg_single(:,2);
Coiltwo_pos(:,i)=Coiltwo_pos_single(:,2);
end


out = Coilone_pos(:,:);
in = Coilone_neg(:,:);
up = Coiltwo_pos(:,:);
down = Coiltwo_neg(:,:);

%----Normalize:-------------
% out = out/sum(out);
% in = in/sum(in);
% up = up/sum(up);
% down = down/sum(down);
%------- above Fermi level
% out = out/sum(out(end-5:end));
% in = in/sum(in(end-5:end));
% up = up/sum(up(end-5:end));
% down = down/sum(down(end-5:end));



strenght = 5;

for i=1:angles_num
up_smooth(:,i) = smooth(up(:,i),strenght);
down_smooth(:,i) = smooth(down(:,i),strenght);

in_smooth(:,i) = smooth(in(:,i),strenght);
out_smooth(:,i) = smooth(out(:,i),strenght);

up_smooth(:,i) = smooth(up(:,i),strenght);
down_smooth(:,i) = smooth(down(:,i),strenght);

end;
%Fermi Energy, anpassen damit Bindungsenergie stimmt !!!
EF = 16.3;

%---------------------------------------------------------------------
energy = -(Coilone_pos_single(:,1) - EF);


%----------spin polarization calculation------------
fact = 1;
for i=1:angles_num
SPIP(:,i)=1/sherman*(up(:,i) - down(:,i))./(up(:,i) + down(:,i));
end;


    
for i=1:angles_num
SPOOP(:,i)=1/sherman*(in(:,i) - out(:,i))./(in(:,i) + out(:,i));    
end;


for i=1:angles_num
                            
SPIP_smooth(:,i) = smooth(SPIP(:,i),strenght);
SPOOP_smooth(:,i) = smooth(SPOOP(:,i),strenght);

end;

for i=1:angles_num
sum_up_down(:,i) = up(:,i) + down(:,i);
sum_up_down_smooth(:,i) = smooth(sum_up_down(:,i),strenght);
sum_in_out(:,i) = in(:,i) + out(:,i);
sum_in_out_smooth(:,i) = smooth(sum_in_out(:,i),strenght);
end;



%calculate back REAL Spectrum out of measured Asymmetry !!!

for i=1:angles_num

up_spectrum(:,i) = (sum_up_down(:,i)/2) .* (1+SPIP(:,i));
up_spectrum_smooth(:,i) = smooth(up_spectrum(:,i),strenght);
up_spectrum_norm(:,i) = up_spectrum(:,i)/max(up_spectrum(:,i));
down_spectrum(:,i) = (sum_up_down(:,i)/2) .* (1-SPIP(:,i));
down_spectrum_smooth(:,i) = smooth(down_spectrum(:,i),strenght);
down_spectrum_norm(:,i) = down_spectrum(:,i)/max(down_spectrum(:,i));
IN_spectrum(:,i) = (sum_in_out(:,i)/2) .* (1+SPOOP(:,i));
IN_spectrum_smooth(:,i) = smooth(IN_spectrum(:,i),strenght);
IN_spectrum_norm(:,i) = IN_spectrum(:,i)/max(IN_spectrum(:,i));
OUT_spectrum(:,i) = (sum_in_out(:,i)/2) .* (1-SPOOP(:,i));
OUT_spectrum_smooth(:,i) = smooth(OUT_spectrum(:,i),strenght);
OUT_spectrum_norm(:,i) = OUT_spectrum(:,i)/max(OUT_spectrum(:,i));

end;

for i=1:angles_num
figure(i); clf;
subplot(3,2,1);
box on;hold on; grid on;
title('raw data: green: UP, blue: DOWN');
%errorbar(energy,up,e_up,'g.');
plot(energy,up(:,i),'g.')
plot(energy,up_smooth(:,i),'g','linewidth',2)
%errorbar(energy,down,e_down,'b.')
plot(energy,down(:,i),'b.')
plot(energy,down_smooth(:,i),'b','linewidth',2)
%plot(energy,sum_up_down,'k.-');
%plot(energy,sum_up_down_smooth,'k-','linewidth',2);
%Legende Farben zu Channels
%h = legend('in - Ch 1','out - Ch 2','up - Ch 3','down - Ch 4',2);
%set(h,'Interpreter','none')
%x-Achse falsch rum: 
set(gca,'xdir','reverse');
%axis tight;

subplot(3,2,2);
box on; hold on; grid on;
title(texlabel('raw data: red: UP, black: DOWN'));
%errorbar(energy,in,e_in,'r.')
plot(energy,in(:,i),'r.')
plot(energy,in_smooth(:,i),'r','linewidth',2)
%errorbar(energy,out,e_out,'k.')
plot(energy,out(:,i),'k.')
plot(energy,out_smooth(:,i),'k','linewidth',2)
%plot(energy,sum_in_out,'k.-');
%plot(energy,sum_in_out_smooth,'k-','linewidth',2);

%Legende: Farben zu Channels
%h = legend('in - Ch 1','out - Ch 2','up - Ch 3','down - Ch 4',2);
%set(h,'Interpreter','none')
set(gca,'xdir','reverse');
%axis tight;

subplot(3,2,3);
hold on; box on; grid on;
title('recalculated Spin Spectra: green: UP; blue: DOWN');
plot(energy,up_spectrum(:,i),'g.-');
plot(energy,up_spectrum_smooth(:,i),'g-','linewidth',3);
plot(energy,down_spectrum(:,i),'b.-');
plot(energy,down_spectrum_smooth(:,i),'b-','linewidth',3);
ylabel('Photoemission signal [counts]','fontsize',12);
set(gca,'xdir','reverse');
%axis tight;

subplot(3,2,4);
hold on; box on; grid on;
title('recalculated Spin Spectra: red: UP; black: DOWN');
plot(energy,IN_spectrum(:,i),'r.-');
plot(energy,IN_spectrum_smooth(:,i),'r-','linewidth',3);
plot(energy,OUT_spectrum(:,i),'k.-');
plot(energy,OUT_spectrum_smooth(:,i),'k-','linewidth',3);
ylabel('Photoemission signal [counts]','fontsize',12);
set(gca,'xdir','reverse');
%axis tight;

%----wichtigster Plot: SPIN POLARIZATION (+- Std Dev)--------%
subplot(3,2,5); hold on; box on; grid on;
title('SPIN-POLARIZATION:  in plane Theta=938','fontsize',12);
plot(energy,SPIP(:,i),'r.-');
plot(energy,SPIP_smooth(:,i),'r-','linewidth',3);
%plot(energy,SPIP_smooth - e_SPIP,'g-','linewidth',3);
%plot(energy,SPIP_smooth + e_SPIP,'g-','linewidth',3);
xlabel('binding energy E_B [eV]','fontsize',12);
ylabel('Spin Polarization','fontsize',12);
set(gca,'xdir','reverse');
ylim([-1 1]);

subplot(3,2,6);
hold on;
box on;
grid on;
title('SPIN-POLARIZATION:   OOP Theta=938','fontsize',12);
plot(energy,SPOOP(:,i),'b.-');
plot(energy,SPOOP_smooth(:,i),'b-','linewidth',3);
%plot(energy,SPOOP_smooth - e_SPOOP,'k-','linewidth',3);
%plot(energy,SPOOP_smooth + e_SPOOP,'k-','linewidth',3);
xlabel('binding energy E_B [eV]','fontsize',12);
ylabel('Spin Polarization','fontsize',12);
set(gca,'xdir','reverse');
ylim([-1 1]);

end;








