clear all;
close all;

file = 283;      %CHANGE


name = sprintf('SPE%05d_00000.txt',file);        
data = dlmread(name,'\t',57,0);
siz=size(data);

EF=16.8;        %CHANGE

[x_step_name, x_step]= textread(name,'%s %f',1,'headerlines',48);
[x_max_name, x_max]= textread(name,'%s %f',1,'headerlines',49);
[x_min_name, x_min]= textread(name,'%s %f',1,'headerlines',50);

[y_step_name, y_step]= textread(name,'%s %f',1,'headerlines',52);
[y_max_name, y_max]= textread(name,'%s %f',1,'headerlines',53);
[y_min_name, y_min]= textread(name,'%s %f',1,'headerlines',54);


DX=0.07567;       %CHANGE shift x
DY=-0.2684;       %CHANGE shift y

x=[x_min:x_step:x_max];
y=[y_min:y_step:y_max];

for i=1:siz(2)
    kx(i) = 0.5123*sqrt(EF)*sin(x(i)*2*pi/360)-DX;
end

for j=1:siz(1)
      ky(j)= 0.5123*sqrt(EF)*sin(y(j)*2*pi/360)-DY; 
end;

[Kx,Ky]=meshgrid(kx,ky);

surf(Kx,Ky,Ky*0,data);

xlabel('k_{|| x} [1/Å]'); ylabel('k_{|| y} [1/Å]')
shading interp; view(2);
axis equal;

% xlim([-0.25 0.25]);
% ylim([-0.3 0.3]);


%% saving as image

% name = [file,'_Fermisurface_2D'];
% name = sprintf(name); 
% saveas(gcf,name,'png');