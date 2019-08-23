Data=importdata('p02b-DaneAll_1-28.dat');

n=1;        %nr petli
offset=0;%18;
dopliku=[];
k=1;
for i=1:(size(Data)*[1;0]) %536
    if Data(i,3)==n
        dopliku(k, 1)=Data(i, 1);
        dopliku(k, 2)=Data(i, 2);
        k=k+1;
    elseif isnan(Data(i,3))
        if (n+offset)<10
            liczba=['0' mat2str(n+offset)];
        else
            liczba= mat2str(n+offset);
        end            
        save(['podzielone\p02-0' liczba '.dat'], 'dopliku','-ASCII')
        n=n+1;
        k=1;
        dopliku=[];
    end
end

if (n+offset)<10
    liczba=['0' mat2str(n+offset)];
else
    liczba= mat2str(n+offset);
end            
save(['podzielone\p02-0' liczba '.dat'], 'dopliku','-ASCII')
n=n+1;
k=1;
dopliku=[];
