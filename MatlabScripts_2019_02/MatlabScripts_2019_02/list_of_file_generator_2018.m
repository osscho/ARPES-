clear all;
close all;
%List_of_files generator:
%---------------------------INPUT-------------------------------------

name_number = 263;   %CHANGE




name = sprintf('SPE%05d_00000.txt',name_number);        
% data = dlmread(name,'\t',57,0);



[step_name, step]= textread(name,'%s %f',1,'headerlines',48);
[angle_max_name, angle_max]= textread(name,'%s %f',1,'headerlines',49);
[angle_min_name, angle_min]= textread(name,'%s %f',1,'headerlines',50);

detector_angle=angle_max-angle_min;



start_number =0;    
start_angle = angle_min;      %CHANGE, can also be ML23 for Ferrum calibration
last_angle = angle_max;       %CHANGE, can also be ML23 for Ferrum Calibration 
angle_step = step;        %CHANGE, can also be delta ML23 for Ferrum calibration

%---------------------------------------------------------------------

size = abs((start_angle-last_angle)/angle_step);
for i = 1:(size+1)
    numbers(i,1) = start_number;
    start_number = start_number + 1;
    
    angles(i,1) = start_angle;
    start_angle = start_angle + angle_step;
    
    list(i,1) = numbers(i,1);
    list(i,2) = angles(i,1);
    list(i,3) = step;        % shows the angle step size in mapping mode
end;

name2 = sprintf('list_of_files_%05d.txt',name_number);
dlmwrite(name2,list,'\t');