%creates matrix F with Gaussian profile,
%centered around the centre of the matrix
%and filters matrix M with the created filter matrix F
%results in matrix N with fewer elements than M
%clear all;
%close all;
%---------------


%file_1 = sprintf('SPE0881_21.txt');                    %BST510 - 6QL

%data = dlmread(file_1,'\t',82,0);


M = data(:,:,1);

% M=NF;
%---------------

xsizeF=10;              %INPUT number of matrix elements in x direction
ysizeF=10;              %INPUT number of matrix elements in y direction
FWHMx=45;               %INPUT full width at half maximum of Gauss profile in x direction
FWHMy=45;               %INPUT full width at half maximum of Gauss profile in y direction

sigx=FWHMx/2.3548;      %determination of sigma for Gauss distribution
sigy=FWHMy/2.3548;      %determination of sigma for Gauss distribution
F=zeros(ysizeF,xsizeF); %creation of matrix with appropriate size
for i=1:1:ysizeF;
    for j=1:1:xsizeF;
        F(i,j)=exp(-0.5*(((i-0.5*ysizeF-0.5)/sigy)^2+((j-0.5*xsizeF-0.5)/sigx)^2));
    end
end
S=sum(sum(F));          %sum of all elements of F
F=F/S;                  %renormalization => sum of all elements is now = 1
N=conv2(M,F,'valid');   %filter of matrix M with matrix F, 'valid': edges are cut away
colormap(bone);
imagesc(N);figure(gcf);
%dlmwrite('filt.txt',N);