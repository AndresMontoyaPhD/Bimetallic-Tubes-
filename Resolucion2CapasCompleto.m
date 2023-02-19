syms v1 v2 E1 E2 a b c d P Pc SR1 SR2 ST1 ST2 SZ1 SZ2 C1 C2 C3 C4 r
set(0,'defaulttextinterpreter','Latex')
%Propiedades de las graficas
width = 5;     % Width in inches
height = 3.09;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 1.5;      % LineWidth
msz =5;       % MarkerSize
%Z direction equations
eqn1=SZ1*pi*(b^2-a^2)+SZ2*pi*(c^2-b^2)==P*pi*a^2; %Sumatorio fuerzas Z
eqn2=(1/E1)*(SZ1-v1*(SR1+ST1))==(1/E2)*(SZ2-v2*(SR2+ST2)); %Deformacion en Z 1=2

%R direction
eqn3=C1*b+C2/b==C3*b+C4/b; %Desplazamiento u(r=b) en 1=2 

%Calculo componentes tensiones R y theta
%Capa 1
eqn4=SR1==(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/r^2)+v1*SZ1/(1-v1);
eqn5=ST1==(E1/(1-v1^2))*(C1*(1+v1)+C2*(1-v1)/r^2)+v1*SZ1/(1-v1);
%Capa 2
eqn6=SR2==(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/r^2)+v2*SZ2/(1-v2);
eqn7=ST2==(E2/(1-v2^2))*(C3*(1+v2)+C4*(1-v2)/r^2)+v2*SZ2/(1-v2);

%Igualdad de tension en R
%En r=b SR1=SR2
eqn8=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/b^2)+v1*SZ1/(1-v1)==(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/b^2)+v2*SZ2/(1-v2);

%Presiones conocidas
%En r=a SR1=-P
eqn9=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/a^2)+v1*SZ1/(1-v1)==-P;
%En r=b SR1=-Pc1
eqn10=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/b^2)+v1*SZ1/(1-v1)==-Pc;
%En r=c SR2=0
eqn11=(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/c^2)+v2*SZ2/(1-v2)==0;

eqns=[eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11];

Sol=solve(eqns,[Pc,SR1,SR2,ST1,ST2,SZ1,SZ2,C1,C2,C3,C4],'Real',true);

%Valores numericos y resolucion

%Prueba para todas las capas 1 mm

%2 layers
E1=131.7e9;
E2=12682000000;
v1=0.274;
v2=0.4;
a=0.003; %inner radius
b=0.0045; %interfase 1 radius
c=0.006; %Outer radius radius
P=20e6; %internal pressure

SZ1=double(subs(Sol.SZ1));
SZ2=double(subs(Sol.SZ2));

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

%Ploteo solucion

%R
figure(1)
box on
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
figure_size=get(gcf,'position');
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
set(gca,'TickLabelInterpreter','latex')
hold on
plot(r1*1000,SR1/10^6,'r-','LineWidth',lw,'MarkerSize',msz)
hold on
plot(r2*1000,SR2/10^6,'r-','LineWidth',lw,'MarkerSize',msz)

%Theta
%figure(2)
plot(r1*1000,ST1/10^6,'g-','LineWidth',lw,'MarkerSize',msz)
hold on
t1=[b,b];
s1=[ST1(size(r1,2)),ST2(1)];
plot(t1*1000,s1/10^6,'g-','LineWidth',lw,'MarkerSize',msz)
plot(r2*1000,ST2/10^6,'g-','LineWidth',lw,'MarkerSize',msz)

%Z
%figure(3)
plot(r1*1000,SZ1R/10^6,'b-','LineWidth',lw,'MarkerSize',msz)
t11=[b,b];
s1=[SZ1,SZ2];
hold on
plot(r2*1000,SZ2R/10^6,'b-','LineWidth',lw,'MarkerSize',msz)
plot(t11*1000,s1/10^6,'b-','LineWidth',lw,'MarkerSize',msz)

%Equivalent Stress
%figure(4)
plot(r1*1000,SEQ1/10^6,'k-','LineWidth',lw,'MarkerSize',msz)
t11=[b,b];
s1=[SEQ1(size(r1,2)),SEQ2(1)];
hold on
plot(r2*1000,SEQ2/10^6,'k-','LineWidth',lw,'MarkerSize',msz)
plot(t11*1000,s1/10^6,'k-','LineWidth',lw,'MarkerSize',msz)

%Resultados Abaqus
SR=xlsread('2layers_StressComponentsTest_15mm','SR');
SR(:,1)=SR(:,1)+a;
%figure(1)
plot(SR(:,1)*1000,SR(:,2)/10^6,'ro','LineWidth',lw,'MarkerSize',msz)

ST=xlsread('2layers_StressComponentsTest_15mm','ST');
ST(:,1)=ST(:,1)+a;
%figure(2)
plot(ST(:,1)*1000,ST(:,2)/10^6,'go','LineWidth',lw,'MarkerSize',msz)

SZ=xlsread('2layers_StressComponentsTest_15mm','SZ');
SZ(:,1)=SZ(:,1)+a;
%figure(3)
plot(SZ(:,1)*1000,SZ(:,2)/10^6,'bo','LineWidth',lw,'MarkerSize',msz)

VMS=xlsread('2layers_StressComponentsTest_15mm','VMS');
VMS(:,1)=VMS(:,1)+a;
%figure(4)
plot(VMS(:,1)*1000,VMS(:,2)/10^6,'ko','LineWidth',lw,'MarkerSize',msz)
ylabel('$\sigma$ [MPa]')
xlabel('r [mm]')
h=legend('$\sigma_r$','$\sigma_\theta$','$\sigma_z$','$\sigma_{eq}$','Interpreter','latex')
h.FontSize=14;