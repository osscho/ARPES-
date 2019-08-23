%Do analizy pêtli ze slope'm
%
%podaæ pthprfx np dla pliku "Disk:Dir/Subdir/Wn-xx.dat" bêdzie
% 'Disk:Dir/Subdir/Wn-'   xx - oznacza kolejny nr pliku
%
%bias dla nasycenia podaæ w procentach np: 95 !bez znaku %!
%
%wskazaæ parê pkt-ów (lewy, lewy) dla 1 prostej
%wskazaæ parê pkt-ów (lewy, prawy) dla 2 prostej
%
%wskazaæ doln¹ asymptotê (lewy) nastêpnie górn¹ (prawy) !dla pêtli normalnych!
%wskazaæ górn¹ asymptotê (lewy) nastêpnie doln¹ (prawy) !dla pêtli odwrotnych!
%
%wskazaæ dolne przeciêcie (lewym) z markerem Rem
%wskazaæ górne przeciêcie (prawym) z markerem Rem
%
%wskazaæ dolne przeciêcie (lewym) z markerem Sat
%wskazaæ dolne przeciêcie (prawym) z markerem Sat

clear all
warning off all
pthprfx = '24-oprac\field_cooling\24fc-'; %'MgO(111)_2-oprac\MgO(111)-';%'24-oprac\24-'; %input('Podaj przednazwê pliku: ','s');
sffx = '.dat';%input('Podaj rozszerzenie nazwy pliku: ','s');
start = input('Podaj numer pocz¹tku serii: ');
koniec = input('Podaj numer koñca serii: ');
%LastSave = input('Podaj nazwe ostatniego zapisanego pliku: ');
%FileExist=exist('LastSave.mat')
%if FileExist==2
%    LastSave=importdata('LastSave.mat');
%else
LastSave=[];
%end
bias = 0.99; %input('Podaj poziom bias: ')/100;
clf reset;
%pause;
try
for i=start:koniec;
    if i==start;
        if i<10
            liczba=['00' mat2str(i)];
        elseif i<100
            liczba=['0' mat2str(i)];
        else 
            liczba=mat2str(i);
        end
        nazwa=char([pthprfx liczba sffx]);
        Dane=importdata(nazwa);
%%{
        clf reset;
        set(gcf, 'Position', [180 47 937 665]);
        subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(Dane(:,1),Dane(:,2),'c-');
        axis([min(Dane(:,1)) max(Dane(:,1)), min(Dane(:,2)) max(Dane(:,2))]);
        %axis([min(Dane(:,2)) max(Dane(:,2)) -12e-5 -7e-5]);
        title(char(['Nazwa pliku: ' nazwa]));
        Slp = [];    % Initially, the list of points is empty.
        n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'ro')
            n = n+1;
            Slp(n,1:2) = [xi, yi];
            hold off
            end
        d=[Slp(2,1)-Slp(1,1), Slp(2,2)-Slp(1,2); Slp(4,1)-Slp(3,1), Slp(4,2)-Slp(3,2)]; %delta
        a(i,1)=(d(1,2)/d(1,1)*sqrt(d(1,2)^2+d(1,1)^2)+d(2,2)/d(2,1)*sqrt(d(2,2)^2+d(2,1)^2))/(sqrt(d(1,2)^2+d(1,1)^2)+sqrt(d(2,2)^2+d(2,1)^2));
        Dane(:,2)=Dane(:,2)-a(i,1)*Dane(:,1);
%%}
        clf reset;
        subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(Dane(:,1),Dane(:,2),'bx-');
        axis([min(Dane(:,1)) max(Dane(:,1)) min(Dane(:,2)) max(Dane(:,2))]);
        %axis([min(Dane(:,2)) max(Dane(:,2)) -12e-5 -7e-5]);
        title(mat2str(['Nazwa pliku: ' nazwa]));
        Norm = [];    % Initially, the list of points is empty.
        n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'ro')
            n = n+1;
            Norm(i,n) = [yi];
            hold off
            end
        clf reset;
        subplot(3,4,1); plot(Norm);
        title('Amplituda')
        bias = [-10.5, bias; 10.5, bias; 10.5, -bias; -10.5, -bias];
        Dane(:,2)=2*(Dane(:,2)-(Norm(i,2)+Norm(i,1))/2)/(Norm(i,2)-Norm(i,1));
        subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(Dane(:,1),Dane(:,2),'k-',[0;0],[-1;1],'m:');
        title(char(['Nazwa pliku: ' nazwa]));
        axis([.1*min(Dane(:,1)) .1*max(Dane(:,1)) -1.02 1.02]);
        Rem = [];    % Initially, the list of points is empty.
        n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'co')
            n = n+1;
            Rem(i,n) = [yi];
            hold off
            end
        subplot(3,4,5); plot(Rem);
        title('Remanencja')
        Sat = [];    % Initially, the list of points is empty.
        subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(Dane(:,1), Dane(:,2),'k-', bias(:,1), bias(:,2), 'b:');
        title(char(['Nazwa pliku: ' nazwa]));
        axis([min(Dane(:,1)) max(Dane(:,1)) -1.02 1.02]);
        n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'mo')
            n = n+1;
            Sat(i,n) = [xi];
            Saty(i,n)= [yi];
            hold off
            end
        subplot(3,4,9); plot(Sat);
        title('Pole nasycenia')
        j=1;
        Dane(:,3)=j;
        %WnAll=[];
        DaneAll=[Dane];
        DaneOld=Dane;
        pause(1);
%{        
        Sat2 = [];    % Initially, the list of points is empty.
        subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(Dane(:,1), Dane(:,2),'k-', bias(:,1), bias(:,2), 'b:');
        title(char(['Nazwa pliku: ' nazwa]));
        axis([min(Dane(:,1)) max(Dane(:,1)) -1.02 1.02]);
        n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'mo')
            n = n+1;
            Sat2(i,n) = [xi];
            hold off
            end
        subplot(3,4,9); plot(Sat2);
        title('Pole nasycenia 2')
        j=1;
        Dane(:,3)=j;
        %WnAll=[];
        DaneAll=[Dane];
        DaneOld=Dane;
        pause(1); 
%} 
    else
        if i<10
            liczba=['00' mat2str(i)];
        elseif i<100
            liczba=['0' mat2str(i)];
        else 
            liczba=mat2str(i);
        end
        nazwa=char([pthprfx liczba sffx]);
        while exist(nazwa)==0
            disp('brak pliku')
            i=i+1;
            if i<10
                liczba=['00' mat2str(i)];
            elseif i<100
                liczba=['0' mat2str(i)];
            else 
                liczba=mat2str(i);
            end            
            nazwa=char([pthprfx liczba sffx])
        end
        %disp('koniec while')
        Dane=importdata(nazwa);
%%{
        clf reset;
        subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(Dane(:,1),Dane(:,2),'c-');
        axis([min(Dane(:,1)) max(Dane(:,1)) min(Dane(:,2)) max(Dane(:,2))]);
        %axis([min(Dane(:,2)) max(Dane(:,2)) -3e-5 0]);
        title(char(['Nazwa pliku: ' nazwa]));
          n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'ro')
            n = n+1;
            Slp(n,1:2) = [xi, yi];
            hold off
            end
        d=[Slp(2,1)-Slp(1,1), Slp(2,2)-Slp(1,2); Slp(4,1)-Slp(3,1), Slp(4,2)-Slp(3,2)]; %delta
        a(i,1)=(d(1,2)/d(1,1)*sqrt(d(1,2)^2+d(1,1)^2)+d(2,2)/d(2,1)*sqrt(d(2,2)^2+d(2,1)^2))/(sqrt(d(1,2)^2+d(1,1)^2)+sqrt(d(2,2)^2+d(2,1)^2));
        Dane(:,2)=Dane(:,2)-a(i,1)*Dane(:,1);
%%}
        clf reset;
        subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(Dane(:,1),Dane(:,2),'bx-');
        axis([min(Dane(:,1)) max(Dane(:,1)) min(Dane(:,2)) max(Dane(:,2))]);
        title(char(['Nazwa pliku: ' nazwa]));
        n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'ro')
            n = n+1;
            Norm(i,n) = [yi];
            hold off
            end
        clf reset;
        subplot(3,4,1); plot(Norm);
        title('Amplituda');
        %axis([start koniec  ]);
        Dane(:,2)=2*(Dane(:,2)-(Norm(i,2)+Norm(i,1))/2)/(Norm(i,2)-Norm(i,1));
        subplot(3,4,[2 3 4 6 7 8 10 11 12]);plot(Dane(:,1),Dane(:,2),'k-',DaneOld(:,1),DaneOld(:,2),'g:',[0;0],[-1;1],'m:');
        axis([.1*min(Dane(:,1)) .1*max(Dane(:,1)) -1.02 1.02]);
        title(char(['Nazwa pliku: ' nazwa]));
        n = 0;     % Loop, picking up the points.
        but = 1;
            while (but == 1)%&&(n > 0);
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'co')
            n = n+1;
            Rem(i,n) = [yi];
            hold off
            end
        subplot(3,4,5); plot(Rem);
        title('Remanencja');
        axis([start-.5 koniec+.5 -1 1]);
        subplot(3,4,[2 3 4 6 7 8 10 11 12]);
            plot(Dane(:,1),Dane(:,2),'k-',DaneOld(:,1),DaneOld(:,2),'g:', bias(:,1), bias(:,2), 'b:', ...
                Sat((i-1), 1), -1.01, '+r', Sat((i-1), 2), 1.01, '+r');
        axis([min(Dane(:,1)) max(Dane(:,1)) -1.02 1.02]);
        title(char(['Nazwa pliku: ' nazwa]));
        n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'mo');
            n = n+1;
            Sat(i,n) = [xi];
            hold off;
            end
        subplot(3,4,9); plot(Sat);
        title('Pole nasycenia');
        %axis([start-.5 koniec+.5 min(Dane(:,2)) max(Dane(:,2))]);
        axis([start-.5 koniec+.5 1.1*min(min(Sat)) 1.1*max(max(Sat))]);
        j=j+1;
        Dane(:,3)=j;
        DaneAll=[DaneAll; NaN NaN NaN; Dane]; % vector [NaN NaN NaN] is a separator
        DaneOld=Dane;
        pause(1);
    end
    subplot(3,1,1); plot(Rem,'o:');
    axis([start-.5 koniec+.5 -1 1]);
    title('Remanencja');
    subplot(3,1,2:3); plot(Sat,'o:'); %hold on; plot(Sat2,'o:k');
    axis([start-.5 koniec+.5 min(Dane(:,1)) max(Dane(:,1))]);
    title('Pole nasycenia');
    %pause;
end
NormSr=(Norm(:,2)-Norm(:,1))/2;
RemSr=(Rem(:,2)-Rem(:,1))/2;
%SatSr=(Sat(:,2)-Sat(:,1))/2;

%RemSr2=(Rem(:,4)-Rem(:,3))/2;
%SatSr2=(Sat(:,4)-Sat(:,3))/2;
%RemSr3=(Rem(:,6)-Rem(:,5))/2;
%SatSr3=(Sat(:,6)-Sat(:,5))/2;

Remplt=[[1:1:1+koniec-start]' 0*[1:1:1+koniec-start]' RemSr(start:koniec,1)];
Satplt1=[[1:1:1+koniec-start]' Sat(start:koniec,1) -1+0*[1:1:1+koniec-start]'];
Satplt2=[[1:1:1+koniec-start]' Sat(start:koniec,2) 1+0*[1:1:1+koniec-start]'];
%Satplt2=[[1:1:1+koniec-start]' SatSr2(start:koniec,1) 1+0*[1:1:1+koniec-start]'];
%Satplt3=[[1:1:1+koniec-start]' SatSr3(start:koniec,1) 1+0*[1:1:1+koniec-start]'];

subplot(3,1,1); plot(NormSr,'o:');
axis([start-.001 koniec min(NormSr) max(NormSr)]);
title('Amplituda');
subplot(3,1,2); plot(RemSr,'o:');
axis([start-.001 koniec 0 1.05*max(RemSr)]);
title('Remanencja');
subplot(3,1,3); plot(Sat(:,1),'o:r'); hold on; plot(Sat(:,2),'o:k');
axis([start-.001 koniec 1.05*max(Sat(:,1)) 1.05*max(Sat(:,2))]);
title('Pole nasycenia');
%{
pause
%i=1
n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,1,2); plot(xi,yi,'rx');
            axis([start-.001 koniec 0 1.05*max(RemSr)]);
            title('Remanencja');%subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'bo');
            n = n+1;
            RemanencjaSr(1,n) = [yi];
            hold off;
            end
 %i=1           
 n = 0;     % Loop, picking up the points.
        but = 1;
            while but == 1
            [xi,yi,but] = ginput(1);
            hold on
            subplot(3,1,3); plot(xi,yi,'gx');%subplot(3,4,[2 3 4 6 7 8 10 11 12]); plot(xi,yi,'go');
            n = n+1;
            SaturacjaSr(1,n) = [yi];
            hold off;
            end

LastSave(TmpDeg+1,:,1)=RemanencjaSr;
LastSave(TmpDeg+1,:,2)=SaturacjaSr;
%}
save(char([pthprfx 'LastSave.mat']),'LastSave','-mat');
disp('+')


            
subplot(1,1,1) ; plot3(DaneAll(:,3),DaneAll(:,1),DaneAll(:,2),'k-');
pause
subplot(1,1,1) ; plot3(DaneAll(:,3),DaneAll(:,1),DaneAll(:,2),'k-')
hold on
plot3(Remplt(:,1), Remplt(:,2), Remplt(:,3),'r-', 'LineWidth', 2)
plot3(Satplt1(:,1), Satplt1(:,2), Satplt1(:,3),'b-', 'LineWidth', 2);%,Satplt2(:,1), Satplt2(:,2), Satplt2(:,3),'g-', Satplt3(:,1), Satplt3(:,2), Satplt3(:,3),'m-');
plot3(Satplt2(:,1), Satplt2(:,2), Satplt2(:,3),'g-', 'LineWidth', 2);%,Satplt2(:,1), Satplt2(:,2), Satplt2(:,3),'g-', Satplt3(:,1), Satplt3(:,2), Satplt3(:,3),'m-');
axis([.5 koniec-start+1+.5 min(DaneAll(:,2)) max(DaneAll(:,2)) -1.2 1.2]);
%ylim([-0.4 0.4]);
%zlim([-1.2 1.2]);
%xlim([.5 koniec-start+1+.5]);
xlabel('kat')
ylabel('pole [T]')
%set(gca,'XTick', 1:18) %0:5:90)
%set(gca,'XTickLabel', 0:20:340)
hold off
%Results=[RemSr(start:koniec,:) SatSr(start:koniec,:) NormSr(start:koniec,:) Rem(start:koniec,:) Sat(start:koniec,:) Norm(start:koniec,:)];
Results=[RemSr(start:koniec,:) NormSr(start:koniec,:) Rem(start:koniec,:) Sat(start:koniec,:) Norm(start:koniec,:)];
disp('Zapis danych... (Enter = zapis ;  Ctrl + c = anuluj)')
pause;
save(char([pthprfx '_res_' mat2str(start) '-' mat2str(koniec) sffx]),'Results','-ASCII');
disp('+')
%save(char([pthprfx '_slope_' mat2str(start) '-' mat2str(koniec) sffx]),'a','-ASCII');
%disp('+')
save(char([pthprfx '_Splt1_' mat2str(start) '-' mat2str(koniec) sffx]),'Satplt1','-ASCII');
disp('+')
save(char([pthprfx '_Splt2_' mat2str(start) '-' mat2str(koniec) sffx]),'Satplt2','-ASCII');
disp('+')
save(char([pthprfx '_Rplt_' mat2str(start) '-' mat2str(koniec) sffx]),'Remplt','-ASCII');
disp('+')
save(char([pthprfx 'DaneAll_' mat2str(start) '-' mat2str(koniec) sffx]),'DaneAll','-ASCII');
disp('OK.')
catch me
end
warning on
