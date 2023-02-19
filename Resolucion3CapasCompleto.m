syms v1 v2 v3 E1 E2 E3 a b c d P Pc1 Pc2 SR1 SR2 SR3 ST1 ST2 ST3 SZ1 SZ2 SZ3 C1 C2 C3 C4 C5 C6 r

%Z direction equations
eqn1=SZ1*pi*(b^2-a^2)+SZ2*pi*(c^2-b^2)+SZ3*pi*(d^2-c^2)==P*pi*a^2; %Sumatorio fuerzas Z
eqn2=(1/E1)*(SZ1-v1*(SR1+ST1))==(1/E2)*(SZ2-v2*(SR2+ST2)); %Deformacion en Z 1=2
eqn3=(1/E2)*(SZ2-v2*(SR2+ST2))==(1/E3)*(SZ3-v3*(SR3+ST3)); %Deformacion en Z 2=3

%R direction
eqn4=C1*b+C2/b==C3*b+C4/b; %Desplazamiento u(r=b) en 1=2 
eqn5=C3*c+C4/c==C5*c+C6/c; %Desplazamiento u(r=c) en 2=3

%Calculo componentes tensiones R y theta
%Capa 1
eqn6=SR1==(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/r^2)+v1*SZ1/(1-v1);
eqn7=ST1==(E1/(1-v1^2))*(C1*(1+v1)+C2*(1-v1)/r^2)+v1*SZ1/(1-v1);
%Capa 2
eqn8=SR2==(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/r^2)+v2*SZ2/(1-v2);
eqn9=ST2==(E2/(1-v2^2))*(C3*(1+v2)+C4*(1-v2)/r^2)+v2*SZ2/(1-v2);
%Capa 3
eqn10=SR3==(E3/(1-v3^2))*(C5*(1+v3)-C6*(1-v3)/r^2)+v3*SZ3/(1-v3);
eqn11=ST3==(E3/(1-v3^2))*(C5*(1+v3)+C6*(1-v3)/r^2)+v3*SZ3/(1-v3);

%Igualdad de tension en R
%En r=b SR1=SR2
eqn12=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/b^2)+v1*SZ1/(1-v1)==(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/b^2)+v2*SZ2/(1-v2);
%En r=c SR2=SR3
eqn13=(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/c^2)+v2*SZ2/(1-v2)==(E3/(1-v3^2))*(C5*(1+v3)-C6*(1-v3)/c^2)+v3*SZ3/(1-v3);

%Presiones conocidas
%En r=a SR1=-P
eqn14=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/a^2)+v1*SZ1/(1-v1)==-P;
%En r=b SR1=-Pc1
eqn15=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/b^2)+v1*SZ1/(1-v1)==-Pc1;
%En r=c SR2=-Pc2
eqn16=(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/c^2)+v2*SZ2/(1-v2)==-Pc2;
%En r=d P=Patm=0, SR3=0
eqn17=(E3/(1-v3^2))*(C5*(1+v3)-C6*(1-v3)/d^2)+v3*SZ3/(1-v3)==0;

eqns=[eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11,eqn12,eqn13,eqn14,eqn15,eqn16,eqn17];

Sol=solve(eqns,[Pc1,Pc2,SR1,SR2,SR3,ST1,ST2,ST3,SZ1,SZ2,SZ3,C1,C2,C3,C4,C5,C6],'Real',true);


%Valores numericos y resolucion

%Prueba para todas las capas 1 mm

%3 layers
E1=131.7e9;
E2=12682000000;
E3=131.7e9;
v1=0.274;
v2=0.4;
v3=0.274;
a=0.003; %inner radius
b=0.004; %interfase 1 radius
c=0.005; %interfase 2 radius
d=0.006;
P=20e6; %internal pressure

SZ1=double(subs(Sol.SZ1));
SZ2=double(subs(Sol.SZ2));
SZ3=double(subs(Sol.SZ3));

%Calculo tensiones Tramo 1
r1=[a:0.00001:b];
SR1=zeros(1,size(r1,2));
ST1=zeros(1,size(r1,2));
SZ1R=zeros(1,size(r1,2));
SEQ1=zeros(1,size(r1,2));
for i=1:size(r1,2)
    r=r1(i);
    SR1(i)=double(subs(Sol.SR1));
    ST1(i)=double(subs(Sol.ST1));
    SZ1R(i)=SZ1;
    SEQ1(i)=sqrt(((SR1(i)-ST1(i))^2+(ST1(i)-SZ1)^2+(SZ1-SR1(i))^2)/2);
end
%Calculo tensiones Tramo 2
r2=[b:0.00001:c];
SR2=zeros(1,size(r2,2));
ST2=zeros(1,size(r2,2));
SZ2R=zeros(1,size(r2,2));
SEQ2=zeros(1,size(r2,2));
for i=1:size(r2,2)
    r=r2(i);
    SR2(i)=double(subs(Sol.SR2));
    ST2(i)=double(subs(Sol.ST2));
    SZ2R(i)=SZ2;
    SEQ2(i)=sqrt(((SR2(i)-ST2(i))^2+(ST2(i)-SZ2)^2+(SZ2-SR2(i))^2)/2);
end
%Calculo tensiones Tramo 3
r3=[c:0.00001:d];
SR3=zeros(1,size(r3,2));
ST3=zeros(1,size(r3,2));
SZ3R=zeros(1,size(r3,2));
SEQ3=zeros(1,size(r3,2));
for i=1:size(r3,2)
    r=r3(i);
    SR3(i)=double(subs(Sol.SR3));
    ST3(i)=double(subs(Sol.ST3));
    SZ3R(i)=SZ3;
    SEQ3(i)=sqrt(((SR3(i)-ST3(i))^2+(ST3(i)-SZ3)^2+(SZ3-SR3(i))^2)/2);
end

%Ploteo solucion

%R
figure(1)
plot(r1,SR1,'b-')
hold on
plot(r2,SR2,'b-')
plot(r3,SR3,'b-')

%Theta
figure(2)
plot(r1,ST1,'b-')
hold on
v1=[b,b];
s1=[ST1(size(r1,2)),ST2(1)];
plot(v1,s1,'b-')
plot(r2,ST2,'b-')
v2=[c,c];
s2=[ST2(size(r2,2)),ST3(1)];
plot(v2,s2,'b-')
plot(r3,ST3,'b')

%Z
figure(3)
plot(r1,SZ1R,'b-')
v1=[b,b];
s1=[SZ1,SZ2];
hold on
plot(r2,SZ2R,'b-')
v2=[c,c];
s2=[SZ2,SZ3];
plot(r3,SZ3R,'b-')
plot(v1,s1,'b-')
plot(v2,s2,'b-')

%Equivalent Stress
figure(4)
plot(r1,SEQ1,'b-')
v1=[b,b];
s1=[SEQ1(size(r1,2)),SEQ2(1)];
hold on
plot(r2,SEQ2,'b-')
v2=[c,c];
s2=[SEQ2(size(r2,2)),SEQ3(1)];
plot(r3,SEQ3,'b-')
plot(v1,s1,'b-')
plot(v2,s2,'b-')


%Resultados Abaqus
SR=xlsread('3layers_StressComponentsTest_1mm','SR');
SR(:,1)=SR(:,1)+a;
figure(1)
plot(SR(:,1),SR(:,2),'r*')

ST=xlsread('3layers_StressComponentsTest_1mm','ST');
ST(:,1)=ST(:,1)+a;
figure(2)
plot(ST(:,1),ST(:,2),'r*')

SZ=xlsread('3layers_StressComponentsTest_1mm','SZ');
SZ(:,1)=SZ(:,1)+a;
figure(3)
plot(SZ(:,1),SZ(:,2),'r*')

VMS=xlsread('3layers_StressComponentsTest_1mm','VMS');
VMS(:,1)=VMS(:,1)+a;
figure(4)
plot(VMS(:,1),VMS(:,2),'r*')

%Prueba para capa cobre 2 mm y las otras de 0.5 mm

%3 layers
E1=131.7e9;
E2=12682000000;
E3=131.7e9;
v1=0.274;
v2=0.4;
v3=0.274;
a=0.003; %inner radius
b=0.0035; %interfase 1 radius
c=0.0055; %interfase 2 radius
d=0.006;
P=20e6; %internal pressure

SZ1=double(subs(Sol.SZ1));
SZ2=double(subs(Sol.SZ2));
SZ3=double(subs(Sol.SZ3));

%Calculo tensiones Tramo 1
r1=[a:0.00001:b];
SR1=zeros(1,size(r1,2));
ST1=zeros(1,size(r1,2));
SZ1R=zeros(1,size(r1,2));
SEQ1=zeros(1,size(r1,2));
for i=1:size(r1,2)
    r=r1(i);
    SR1(i)=double(subs(Sol.SR1));
    ST1(i)=double(subs(Sol.ST1));
    SZ1R(i)=SZ1;
    SEQ1(i)=sqrt(((SR1(i)-ST1(i))^2+(ST1(i)-SZ1)^2+(SZ1-SR1(i))^2)/2);
end
%Calculo tensiones Tramo 2
r2=[b:0.00001:c];
SR2=zeros(1,size(r2,2));
ST2=zeros(1,size(r2,2));
SZ2R=zeros(1,size(r2,2));
SEQ2=zeros(1,size(r2,2));
for i=1:size(r2,2)
    r=r2(i);
    SR2(i)=double(subs(Sol.SR2));
    ST2(i)=double(subs(Sol.ST2));
    SZ2R(i)=SZ2;
    SEQ2(i)=sqrt(((SR2(i)-ST2(i))^2+(ST2(i)-SZ2)^2+(SZ2-SR2(i))^2)/2);
end
%Calculo tensiones Tramo 3
r3=[c:0.00001:d];
SR3=zeros(1,size(r3,2));
ST3=zeros(1,size(r3,2));
SZ3R=zeros(1,size(r3,2));
SEQ3=zeros(1,size(r3,2));
for i=1:size(r3,2)
    r=r3(i);
    SR3(i)=double(subs(Sol.SR3));
    ST3(i)=double(subs(Sol.ST3));
    SZ3R(i)=SZ3;
    SEQ3(i)=sqrt(((SR3(i)-ST3(i))^2+(ST3(i)-SZ3)^2+(SZ3-SR3(i))^2)/2);
end

%Ploteo solucion
set(0,'defaulttextinterpreter','Latex')
%Propiedades de las graficas
width = 5;     % Width in inches
height = 3.09;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 1.5;      % LineWidth
msz =5;       % MarkerSize
%R
figure(5)
box on
hold on
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
figure_size=get(gcf,'position');
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
set(gca,'TickLabelInterpreter','latex')
plot(r1*1000,SR1/10^6,'r-','LineWidth',lw,'MarkerSize',msz)
plot(r2*1000,SR2/10^6,'r-','LineWidth',lw,'MarkerSize',msz)
plot(r3*1000,SR3/10^6,'r-','LineWidth',lw,'MarkerSize',msz)

%Theta
%figure(6)
plot(r1*1000,ST1/10^6,'g-','LineWidth',lw,'MarkerSize',msz)
hold on
v1=[b,b];
s1=[ST1(size(r1,2)),ST2(1)];
plot(v1*1000,s1/10^6,'g-','LineWidth',lw,'MarkerSize',msz)
plot(r2*1000,ST2/10^6,'g-','LineWidth',lw,'MarkerSize',msz)
v2=[c,c];
s2=[ST2(size(r2,2)),ST3(1)];
plot(v2*1000,s2/10^6,'g-','LineWidth',lw,'MarkerSize',msz)
plot(r3*1000,ST3/10^6,'g-','LineWidth',lw,'MarkerSize',msz)

%Z
%figure(7)
plot(r1*1000,SZ1R/10^6,'b-','LineWidth',lw,'MarkerSize',msz)
v1=[b,b];
s1=[SZ1,SZ2];
hold on
plot(r2*1000,SZ2R/10^6,'b-','LineWidth',lw,'MarkerSize',msz)
v2=[c,c];
s2=[SZ2,SZ3];
plot(r3*1000,SZ3R/10^6,'b-','LineWidth',lw,'MarkerSize',msz)
plot(v1*1000,s1/10^6,'b-','LineWidth',lw,'MarkerSize',msz)
plot(v2*1000,s2/10^6,'b-','LineWidth',lw,'MarkerSize',msz)

%Equivalent Stress
%figure(8)
plot(r1*1000,SEQ1/10^6,'k-','LineWidth',lw,'MarkerSize',msz)
v1=[b,b];
s1=[SEQ1(size(r1,2)),SEQ2(1)];
hold on
plot(r2*1000,SEQ2/10^6,'k-','LineWidth',lw,'MarkerSize',msz)
v2=[c,c];
s2=[SEQ2(size(r2,2)),SEQ3(1)];
plot(r3*1000,SEQ3/10^6,'k-','LineWidth',lw,'MarkerSize',msz)
plot(v1*1000,s1/10^6,'k-','LineWidth',lw,'MarkerSize',msz)
plot(v2*1000,s2/10^6,'k-','LineWidth',lw,'MarkerSize',msz)


%Resultados Abaqus
SR=xlsread('3layers_StressComponentsTest_2mm','SR');
SR(:,1)=SR(:,1)+a;
%figure(5)
plot(SR(:,1)*1000,SR(:,2)/10^6,'ro','LineWidth',lw,'MarkerSize',msz)

ST=xlsread('3layers_StressComponentsTest_2mm','ST');
ST(:,1)=ST(:,1)+a;
%figure(6)
plot(ST(:,1)*1000,ST(:,2)/10^6,'go','LineWidth',lw,'MarkerSize',msz)

SZ=xlsread('3layers_StressComponentsTest_2mm','SZ');
SZ(:,1)=SZ(:,1)+a;
%figure(7)
plot(SZ(:,1)*1000,SZ(:,2)/10^6,'bo','LineWidth',lw,'MarkerSize',msz)

VMS=xlsread('3layers_StressComponentsTest_2mm','VMS');
VMS(:,1)=VMS(:,1)+a;
%figure(8)
plot(VMS(:,1)*1000,VMS(:,2)/10^6,'ko','LineWidth',lw,'MarkerSize',msz)
ylabel('$\sigma$ [MPa]')
xlabel('r [mm]')
h=legend('$\sigma_r$','$\sigma_\theta$','$\sigma_z$','$\sigma_{eq}$','Interpreter','latex')
h.FontSize=14;
ylim([-20 100])