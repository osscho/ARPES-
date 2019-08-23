%determines the 2D curvature of an array input: matrix N which should be low-pass filtered,
%(e.g., with ml_create_matrix_gauss_and_filter.m) computes numerical derivatives up to second order
%and from this the curvatures according to Zhang, RSI 82,043712 (2011)
%at the end, filtering with same matrix F
%
%first derivative computed by using 5x5 matrix
%and averaging over the two non-cancelling possible contributions:
%D1=1/2*[1 1 1 1 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; -1 -1 -1 -1 -1]/20
%+1/2*[0 0 0 0 0; 1 1 1 1 1; 0 0 0 0 0; -1 -1 -1 -1 -1; 0 0 0 0 0]/10
D1=[1 1 1 1 1;2 2 2 2 2; 0 0 0 0 0; -2 -2 -2 -2 -2; -1 -1 -1 -1 -1]/40;
%second derivative, similar averaging as above
%D2=1/2*[1 1 1 1 1; 0 0 0 0 0; -2 -2 -2 -2 -2; 0 0 0 0 0; 1 1 1 1 1]/20
%+1/2*[0 0 0 0 0; 1 1 1 1 1; -2 -2 -2 -2 -2; 1 1 1 1 1; 0 0 0 0 0]/5
D2=[1 1 1 1 1; 4 4 4 4 4; -10 -10 -10 -10 -10;4 4 4 4 4; 1 1 1 1 1]/40;
d1=[1 1 1; 0 0 0; -1 -1 -1]/6;%short version for mixed derivatives
Dy=D1;%creation of matrices for derivatives for different orders and directions
D2y=D2;
Dx=D1';
D2x=D2';
dy=d1;
dx=d1';
%computing the derivatives, "valid" means that the boundaries are cut off
%result is 4 pixels smaller in each direction (important for scaling)
Fx=conv2(N,Dx,'valid');
F2x=conv2(N,D2x,'valid');
Fy=conv2(N,Dy,'valid');
F2y=conv2(N,D2y,'valid');
fx=conv2(N,dx,'valid');
fy=conv2(N,dy,'valid');
Fxy=conv2(fx,dy,'valid');
Fyx=conv2(fy,dx,'valid');
F2xy=(Fxy+Fyx)/2;
[m,n]=size(Fx);
%computation of first estimate for Cx, Cy (s. Zhang)
%Cx*df/dx should be on the order of 1
vx=diag(Fx'*Fx)/m;
vy=diag(Fy'*Fy)/m;

scale=0.05;      %INPUT relative weight of Laplace (ca. 0.01) vs. highlighting the zero slope parts (ca. 10)
weight=1;        %INPUT relative weight of derivatives along x (rows) or y (columns), set this according to spectrum details and noise levels

Cx=scale*weight*sqrt(n/(vx'*vx));
Cy=(scale/weight)*sqrt(n/(vy'*vy));
O=ones(m,n);
Laplace=-(Cx*F2x+Cy*F2y);%Laplace operator for comparison
zaehler=(O+Cx*Fx.*Fx).*F2y*Cy-2*Cx*Cy*Fx.*Fy.*F2xy+(O+Cy*Fy.*Fy).*F2x*Cx;
nenner=(O+Cx*Fx.*Fx+Cy*Fy.*Fy).^(3/2);
curv=-zaehler./nenner;%this is the 2D curvature, according to Zhang
colormap(bone);
%imagesc(Laplace);figure(gcf);%just for comparison, if you want
NF=conv2(curv,F,'valid');%low-pass filter with matrix F
%imagesc(curv);figure(gcf);%display of matrix "curv",
%good for judging the best setting of "scale" and "weight"


figure(1);
subplot(1,3,3);
    imagesc(real2rgb(NF,colorm,[0 upper*max(max(NF))]));figure(gcf); %display of 2D curvature, again low-pass filtered
    set(gca,'ydir','normal');




% % ---------- MDC and EDC ---------- %
% 
% siz = size(data_1); channels = siz(1);     %Energy channels: e.g. 1..480
% 
% 
% figure(2); hold on;
% energy_step = 4;        % <- Nr of integrated energies
% for i = 1:siz(1)
%            
%     if mod(i,energy_step) == 0 && i < siz(1)-30  
%     
%         title('Momentum Distribution Curves MDCs');
%              plot(NF(i,:) + i/100 );
%     end;
% end;