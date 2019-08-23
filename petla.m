clear all;
clf reset;
wsp_abs=0.480;%tyle procent dolnej warstwy widzimy w MOKE
K1=-4.9e5;%erg/cm3
K2=-1.5e5;%erg/cm3
Ku1=-0.00e5;%erg/cm3
Ku2=-0.0e5;%erg/cm3
J1=-0.02; %erg/cm2 
J2=-0.0 %erg/cm2<-------
n=180%input('liczba krokow - na ile dzielimy k¹ty: ') %30; 
krokpola=200%input('liczba krokow pola: ')
d=[300 100]*1e-8;
omega=45;%kat anizotropii jednoosiowej wzgledem kierunku latwego
%stale
%d=[300 100]*1e-8;                 %grubosci poszczegolnych warstw[cm]
%% interface,...,surface with relation to [001] direction
dtotal=d(1)+d(2);
fi_down=[0 0];   %wektor dolnych granic katow poszczegolnych warstw - wymiar musi byc taki jak w wektorze d[] bo to jest liczba warstw w obu przypadkach
fi_up=[360 360];   %wektor dolnych granic katow poszczegolnych warstw - wymiar musi byc taki jak w wektorze d[] bo to jest liczba warstw w obu przypadkach
s=size(d)*[0 1]';
%n=180%input('liczba krokow - na ile dzielimy k¹ty: ') %30;   %%tu wpisz liczbe krokow dla kazdego z parametrow
%krok pola
%krokpola=70%input('liczba krokow pola: ')
%%teraz przedzialy zmiennosci parametrow:
data=importdata('petla103.dat');
s=size(data);
max(max(data))
for l=1:s(1)
    h_exp(l)=data(l,2);
    m_exp(l)=data(l,1);
    l=l+1;
end
alfa=90;
m1=1700;
m2=1700;
hmaxexp=max(h_exp);
Hmax=hmaxexp; % 
deltaH=Hmax/krokpola; %
%stale
%%%
h(1)=Hmax;
Fi1(1)=alfa;
Fi2(1)=alfa;
namagnesowanie(1)=1;
for i=1:2*krokpola
   h(i+1)=Hmax-i*deltaH; %wyliczenie pola w danym kroku
E=zeros([n n]);
                for m=0:n-1
                fi_1=fi_down(1)+m*(fi_up(1)-fi_down(1))/(n-1); %interface angle
                m=m+1;
                    for g=0:n-1
                    fi_2=fi_down(2)+g*(fi_up(2)-fi_down(2))/(n-1);  %surface angle
                    g=g+1;
                    fi_1out(m,g)=fi_1;  %katy warstwy dolnej
                    fi_2out(m,g)=fi_2;  %katy warstwy powierzchniowej
                    %liczenie energii!!
                        E(m,g)=E(m,g)-K1/4*d(1)*sin((pi/180)*2*fi_1out(m,g))^2-K2/4*d(2)*sin((pi/180)*2*fi_2out(m,g))^2+Ku1*d(1)*sin(pi/180*(fi_1out(m,g)-omega))^2+Ku2*d(2)*sin(pi/180*(fi_2out(m,g)-omega))^2-J1*cos(pi/180*(fi_1out(m,g)-fi_2out(m,g)))-J2*cos(pi/180*(fi_1out(m,g)-fi_2out(m,g)))^2-m1*d(1)*h(i)*cos((pi/180)*(fi_1out(m,g)-alfa))-m2*d(2)*h(i)*cos((pi/180)*(fi_2out(m,g)-alfa));
                    end
                end
            E;
           minimum_energii=(min(min(E))); %%petla do szukania indeksu w tablicy energi dla ktorego mamy minimum energii
              for m=1:n
                  for g=1:n
                    if ~(E(m,g)-minimum_energii)
                    x1=m;
                    x2=g;
                    end
                    g=g+1;
                  end
              m=m+1;
              end
fi_1out;
fi_2out;
Energia=E;
%plot3(fi_1out,fi_2out,E)
m;
g;
x1;
x2;
E=E(x1,x2);
kat_dolnej_warstwy=fi_1out(x1,x2);
kat_na_powierzchni=fi_2out(x1,x2);
MinE(i+1)=minimum_energii;
%figurepalette;
%surfc (fi_1out, fi_2out, Energia ,'DisplayName', 'fi_1out, fi_2out, Energia'); figure(gcf);
%pcolor(fi_1out, fi_2out, Energia)
%plot3(fi_1out, fi_2out, Energia)
%text(kat_dolnej_warstwy,kat_na_powierzchni,'\leftarrow minE','VerticalAlignment','Middle','FontSize',38) 
%close
%pause%(0.1)
namagnesowanie1(i+1)=wsp_abs*d(1)*cos(pi/180*(kat_dolnej_warstwy-alfa))/(wsp_abs*d(1)+d(2));
namagnesowanie2(i+1)=d(2)*cos(pi/180*(kat_na_powierzchni-alfa))/(wsp_abs*d(1)+d(2));
namagnesowanie(i+1)=(wsp_abs*d(1)*cos(pi/180*(kat_dolnej_warstwy-alfa))+d(2)*cos(pi/180*(kat_na_powierzchni-alfa)))/(wsp_abs*d(1)+d(2));
Fi1(i+1)=kat_dolnej_warstwy;
Fi2(i+1)=kat_na_powierzchni;
X1(i+1)=300*cos(pi/180*Fi1(i+1));
X2(i+1)=100*cos(pi/180*Fi2(i+1));
Y1(i+1)=300*sin(pi/180*Fi1(i+1));
Y2(i+1)=100*sin(pi/180*Fi2(i+1));
%[X1(i+1),Y1(i+1)]=pol2cart(Fi1(i+1),1);
%[X2(i+1),Y2(i+1)]=pol2cart(Fi2(i+1),1);
end
%%teraz zaczynamy to samo co wyzej - drugi kierunek pola!
%m
%g
%wskaznik=0;
%for i=2*krokpola+1:4*krokpola
%    wskaznik=wskaznik+1;
%   h(i+1)=-Hmax+wskaznik*deltaH; %wyliczenie pola w danym kroku
%   E=zeros([n n]);
%                for m=0:n-1
%                fi_1=fi_down(1)+m*(fi_up(1)-fi_down(1))/(n-1); %interface angle
%                m=m+1;
%                    for g=0:n-1
%                    fi_2=fi_down(2)+g*(fi_up(2)-fi_down(2))/(n-1);  %surface angle
%                    g=g+1;
%                    fi_1out(m,g)=fi_1;  %katy warstwy dolnej
%                    fi_2out(m,g)=fi_2;  %katy warstwy powierzchniowej
%                    %liczenie energii!!
%                        E(m,g)=E(m,g)-K1/4*d(1)*sin((pi/180)*2*fi_1out(m,g))^2+K2/4*d(2)*sin((pi/180)*2*fi_2out(m,g))^2+Ku1*d(1)*sin(pi/180*fi_1out(m,g))^2+Ku2*d(2)*sin(pi/180*fi_2out(m,g))^2-J1*cos(pi/180*(fi_1out(m,g)-fi_2out(m,g)))-J2*cos(pi/180*(fi_1out(m,g)-fi_2out(m,g)))^2-m1*d(1)*h(i)*cos((pi/180)*(fi_1out(m,g)-alfa))-m2*d(2)*h(i)*cos((pi/180)*(fi_2out(m,g)-alfa));
%                    end
%                end
%            E;
%           minimum_energii=(min(min(E))); %%petla do szukania indeksu w tablicy energi dla ktorego mamy minimum energii
%              for m=1:n
%                  for g=1:n
%                    if ~(E(m,g)-minimum_energii)
%                    x1=m;
%                    x2=g;
%                    end
%                    g=g+1;
%                  end
%              m=m+1;
%              end
%fi_1out;
%fi_2out;
%Energia=E;
%plot3(fi_1out,fi_2out,E)
%m;
%g;
%x1;
%x2;
%E=E(x1,x2);
%kat_dolnej_warstwy=fi_1out(x1,x2);
%kat_na_powierzchni=fi_2out(x1,x2);
%minimum_energii;
%namagnesowanie(i+1)=(wsp_abs*d(1)*cos(pi/180*(kat_dolnej_warstwy-alfa))+d(2)*cos(pi/180*(kat_na_powierzchni-alfa)))/(wsp_abs*d(1)+d(2));
%end
%subplot(2,1,1)
clc
figure(1)
plot(h,namagnesowanie,'b.-',h_exp,m_exp,'r.','MarkerSize',5)
daneplot=[h',namagnesowanie'];
roznica=Fi1-Fi2;
danekaty=[h',roznica'];

TMR=[h',1-cosd(roznica)']
save('roznicakatow.dat','danekaty','-ASCII')
save('TMR.dat','TMR','-ASCII')
save('dane.dat','daneplot','-ASCII')
set(gcf,'Units','normalized','Position',[0 0.08 1 .8])
axis([-1.1*Hmax 1.1*Hmax,min(namagnesowanie)-.05*(max(namagnesowanie)-min(namagnesowanie)) max(namagnesowanie)+.05*(max(namagnesowanie)-min(namagnesowanie))])
title('Namagnesowanie ukladu')
grid on

hold all

t(2,1)={num2str(wsp_abs)};
t(4,1)={num2str(K1,'%1.1e')};
t(6,1)={num2str(K2,'%1.1e')};
t(2,2)={num2str(J1)};
t(4,2)={num2str(J2)};
t(6,2)={num2str(krokpola)};
t(1,1)={'wspolczynnik absorpcji'};
t(3,1)={'K1'};
t(5,1)={'K2'};
t(1,2)={'J1'};
t(3,2)={'J2'};
t(5,2)={'Krok pola'};


s=text(.9*min(h),.5*max(namagnesowanie),t);
set(s,'FontSize',12,'BackgroundColor',[1 1 1]);

saveas(gcf,'103.bmp')

figure(2)
plot(h,1-cosd(roznica),'k.-')
axis([min(h) max(h) -.2 2.2])
grid on
saveas(gcf,'TMR.bmp')
%pause
%for i=1:2*krokpola
%compass(X1(i),Y1(i))
%hold all
%compass(X1(i),Y1(i),'r-')
%hold on
%compass(X2(i),Y2(i),'b-')
%hold off
%pause(0.4)
%end