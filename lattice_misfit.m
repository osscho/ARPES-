a_Au=2.883;% [A]
a_Fe=2.866;% [A]
a_Au_m=a_Au*10^-4;
a_Fe_m=a_Fe*10^-4;
length=10;% [microns]
cell_count=length/a_Fe_m; %[microns]

diff=cell_count*a_Au_m-cell_count*a_Fe_m; %[microns]
fprintf('Difference at the distance of 10 microns [in microns] %f', diff);

