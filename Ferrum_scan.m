clear all;
close all;


strenght = 5;
sherman=0.27;
C = {'r','k','b','g','y','m','c',[.5 .6 .7],[.8 .2 .6],[.1 .2 .6],[.8 .2 .1],[.8 .8 .6],[.5 .5 .6],[.1 .3 .6],[.0 .0 .6],[.2,.2,.2]} % Cell array of colors.
files_up=[8784:1:8799];

num_of_files_up=numel(files_up);
for i=1:num_of_files_up
    f=files_up(i);
    file_name = sprintf('SPE%d_07.mot',f); 
    data_spin_up(:,:,i)=dlmread(file_name,'\t',1,0);
end;





for i=1:num_of_files_up
plot(smooth(data_spin_up(:,1,i),strenght),smooth(data_spin_up(:,2,i),strenght),'color',C{i},'marker','o', 'LineWidth',2)
set(gca,'FontSize',10);
%axis([15.0,17.5,0,0.6e4]);
grid on;
%legend('3563','3573');
title('coil2 positive','FontSize',18 ); xlabel('Kinetic energy [eV]','FontSize',14); ylabel('Counts [arb.u]','FontSize',14)
hold on;
end;

files_down=[8800:1:8815];
num_of_files_down=numel(files_down);
for i=1:num_of_files_down
    f=files_down(i);
    file_name = sprintf('SPE%d_07.mot',f); 
    data_spin_down(:,:,i)=dlmread(file_name,'\t',1,0);
end;


for i=1:num_of_files_up
plot(smooth(data_spin_down(:,1,i),strenght),smooth(data_spin_down(:,2,i),strenght),'color',C{i},'marker','+', 'LineWidth',2)
set(gca,'FontSize',10);
%axis([15.0,17.5,0,0.6e4]);
grid on;
%legend('3563','3573');
title('up&down','FontSize',18 ); xlabel('Kinetic energy [eV]','FontSize',14); ylabel('Counts [arb.u]','FontSize',14)
hold on;
end;



for i=1:num_of_files_up
    figure();
plot(smooth(data_spin_down(:,1,i),strenght),smooth(data_spin_down(:,2,i),strenght),'color',C{i},'marker','+', 'LineWidth',2)
set(gca,'FontSize',10);
%axis([15.0,17.5,0,0.6e4]);
grid on;
%legend('3563','3573');
title('up&down','FontSize',18 ); xlabel('Kinetic energy [eV]','FontSize',14); ylabel('Counts [arb.u]','FontSize',14)
hold on;
plot(smooth(data_spin_up(:,1,i),strenght),smooth(data_spin_up(:,2,i),strenght),'color',C{i},'marker','o', 'LineWidth',2)
end;



figure();


for i=1:num_of_files_up

    SPUD(:,i)=1/sherman*(smooth(data_spin_down(:,2,i),strenght)- smooth(data_spin_up(:,2,i),strenght))./(smooth(data_spin_down(:,2,i),strenght)+smooth(data_spin_up(:,2,i),strenght));

plot(data_spin_down(:,1,i),smooth(SPUD(:,i),strenght),'color',C{i},'LineWidth',2);
hold on; grid on;
end;
figure();

for i=1:num_of_files_up
UP_spectrum(:,i) = ((smooth(data_spin_down(:,2,i),strenght)+smooth(data_spin_up(:,2,i),strenght))/2) .* (1+SPUD(:,i));
DOWN_spectrum(:,i) = ((smooth(data_spin_down(:,2,i),strenght)+smooth(data_spin_up(:,2,i),strenght))/2) .* (1-SPUD(:,i));
plot(data_spin_down(:,1,i),UP_spectrum(:,i),'color',C{i},'LineWidth',2);
hold on;
end;