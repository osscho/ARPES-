close all;

hold on;

files=[3563:3573]
%files=[2386:1:2387]
file_no=numel(files);
C = {'r','k','b','g','y','m','c',[.5 .6 .7],[.8 .2 .6],[.1 .2 .6],[.8 .2 .1],[.8 .8 .6]} % Cell array of colors.
cmap = colormap(jet); 
for k=1:file_no
    f=files(k);

file_name = sprintf('SPE%d_10.txt',f);
data_EDC(:,:) = dlmread(file_name,'\t',82,0); %1000X256
size_data=size(data_EDC);
EDC_sum=zeros(1,size_data(1));
for j=1:size_data(1)
 for i=2:247
 EDC_sum(1,j)=EDC_sum(1,j)+data_EDC(j,i);
 end;
end;


g=(i-1914)*0.0001;
if f==3061
    EDC_sum(1,:)=EDC_sum(1,:)*0.55;
elseif f==2321
    EDC_sum(1,:)=EDC_sum(1,:)*20;
else
    EDC_sum(1,:)=EDC_sum(1,:);
end
plot(data_EDC(:,1),EDC_sum(1,:),'color',C{k},'marker','o', 'LineWidth',2)
set(gca,'FontSize',10);
%axis([15.0,17.5,0,0.6e4]);
grid on;
legend('3563','3573');
title('EM13','FontSize',18 ); xlabel('Kinetic energy [eV]','FontSize',14); ylabel('Counts [arb.u]','FontSize',14)
hold on;
end;

%saveas(gcf, 'EM13 3055 3061', 'png')

