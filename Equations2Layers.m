syms v1 v2 E1 E2 a b c d P Pc SR1 SR2 ST1 ST2 SZ1 SZ2 C1 C2 C3 C4 r
set(0,'defaulttextinterpreter','Latex')
%Figure Properties
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

%Equations to calculate R and theta stress components
%Layer 1
eqn4=SR1==(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/r^2)+v1*SZ1/(1-v1);
eqn5=ST1==(E1/(1-v1^2))*(C1*(1+v1)+C2*(1-v1)/r^2)+v1*SZ1/(1-v1);
%Layer 2
eqn6=SR2==(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/r^2)+v2*SZ2/(1-v2);
eqn7=ST2==(E2/(1-v2^2))*(C3*(1+v2)+C4*(1-v2)/r^2)+v2*SZ2/(1-v2);

%Stress equilibrium in R direction between layers
%In r=b SR1=SR2
eqn8=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/b^2)+v1*SZ1/(1-v1)==(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/b^2)+v2*SZ2/(1-v2);

%Known pressure
%In r=a SR1=-P
eqn9=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/a^2)+v1*SZ1/(1-v1)==-P;
%In r=b SR1=-Pc1
eqn10=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/b^2)+v1*SZ1/(1-v1)==-Pc;
%In r=c SR2=0
eqn11=(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/c^2)+v2*SZ2/(1-v2)==0;

%Solving equations
eqns=[eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11];

Sol=solve(eqns,[Pc,SR1,SR2,ST1,ST2,SZ1,SZ2,C1,C2,C3,C4],'Real',true);

%Numerical values and plotting

%Test all layers 1 mm thickness

%2 layers
E1=131.7e9; %Young modulus material 1
E2=12682000000; %Young modulus material 2
v1=0.274; %Poisson's ratio material 1
v2=0.4; %Poisson's ratio material 2
a=0.003; %inner radius
b=0.0045; %interfase 1 radius
c=0.006; %Outer radius radius
P=20e6; %internal pressure

%Axial stress constant through layers
SZ1=double(subs(Sol.SZ1));
SZ2=double(subs(Sol.SZ2));

%Stress R theta layer eq 1
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
%Stress R theta eq layer 2
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

%Plotting solution

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
figure(2)
plot(r1*1000,ST1/10^6,'g-','LineWidth',lw,'MarkerSize',msz)
hold on
t1=[b,b];
s1=[ST1(size(r1,2)),ST2(1)];
plot(t1*1000,s1/10^6,'g-','LineWidth',lw,'MarkerSize',msz)
plot(r2*1000,ST2/10^6,'g-','LineWidth',lw,'MarkerSize',msz)

%Z
figure(3)
plot(r1*1000,SZ1R/10^6,'b-','LineWidth',lw,'MarkerSize',msz)
t11=[b,b];
s1=[SZ1,SZ2];
hold on
plot(r2*1000,SZ2R/10^6,'b-','LineWidth',lw,'MarkerSize',msz)
plot(t11*1000,s1/10^6,'b-','LineWidth',lw,'MarkerSize',msz)

%Equivalent Stress
figure(4)
plot(r1*1000,SEQ1/10^6,'k-','LineWidth',lw,'MarkerSize',msz)
t11=[b,b];
s1=[SEQ1(size(r1,2)),SEQ2(1)];
hold on
plot(r2*1000,SEQ2/10^6,'k-','LineWidth',lw,'MarkerSize',msz)
plot(t11*1000,s1/10^6,'k-','LineWidth',lw,'MarkerSize',msz)