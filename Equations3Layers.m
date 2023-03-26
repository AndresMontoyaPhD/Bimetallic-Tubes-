syms v1 v2 v3 E1 E2 E3 a b c d P Pc1 Pc2 SR1 SR2 SR3 ST1 ST2 ST3 SZ1 SZ2 SZ3 C1 C2 C3 C4 C5 C6 r

%Z direction equations
eqn1=SZ1*pi*(b^2-a^2)+SZ2*pi*(c^2-b^2)+SZ3*pi*(d^2-c^2)==P*pi*a^2; %Sumatorio fuerzas Z
eqn2=(1/E1)*(SZ1-v1*(SR1+ST1))==(1/E2)*(SZ2-v2*(SR2+ST2)); %Deformacion en Z 1=2
eqn3=(1/E2)*(SZ2-v2*(SR2+ST2))==(1/E3)*(SZ3-v3*(SR3+ST3)); %Deformacion en Z 2=3

%R direction
eqn4=C1*b+C2/b==C3*b+C4/b; %Desplazamiento u(r=b) en 1=2 
eqn5=C3*c+C4/c==C5*c+C6/c; %Desplazamiento u(r=c) en 2=3

%Equations to calculate R and Theta stress components
%Layer 1
eqn6=SR1==(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/r^2)+v1*SZ1/(1-v1);
eqn7=ST1==(E1/(1-v1^2))*(C1*(1+v1)+C2*(1-v1)/r^2)+v1*SZ1/(1-v1);
%Layer 2
eqn8=SR2==(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/r^2)+v2*SZ2/(1-v2);
eqn9=ST2==(E2/(1-v2^2))*(C3*(1+v2)+C4*(1-v2)/r^2)+v2*SZ2/(1-v2);
%Layer 3
eqn10=SR3==(E3/(1-v3^2))*(C5*(1+v3)-C6*(1-v3)/r^2)+v3*SZ3/(1-v3);
eqn11=ST3==(E3/(1-v3^2))*(C5*(1+v3)+C6*(1-v3)/r^2)+v3*SZ3/(1-v3);

%Stress equilibrium in R direction between layers
%In r=b SR1=SR2
eqn12=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/b^2)+v1*SZ1/(1-v1)==(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/b^2)+v2*SZ2/(1-v2);
%In r=c SR2=SR3
eqn13=(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/c^2)+v2*SZ2/(1-v2)==(E3/(1-v3^2))*(C5*(1+v3)-C6*(1-v3)/c^2)+v3*SZ3/(1-v3);

%Known pressure
%In r=a SR1=-P
eqn14=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/a^2)+v1*SZ1/(1-v1)==-P;
%In r=b SR1=-Pc1
eqn15=(E1/(1-v1^2))*(C1*(1+v1)-C2*(1-v1)/b^2)+v1*SZ1/(1-v1)==-Pc1;
%In r=c SR2=-Pc2
eqn16=(E2/(1-v2^2))*(C3*(1+v2)-C4*(1-v2)/c^2)+v2*SZ2/(1-v2)==-Pc2;
%In r=d P=Patm=0, SR3=0
eqn17=(E3/(1-v3^2))*(C5*(1+v3)-C6*(1-v3)/d^2)+v3*SZ3/(1-v3)==0;

%Solving system of equations
eqns=[eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11,eqn12,eqn13,eqn14,eqn15,eqn16,eqn17];

Sol=solve(eqns,[Pc1,Pc2,SR1,SR2,SR3,ST1,ST2,ST3,SZ1,SZ2,SZ3,C1,C2,C3,C4,C5,C6],'Real',true);


%Numerical values and resolution

%Test for layers of 1 mm, layer 1 and 3 same material

%3 layers
E1=131.7e9; %Young modulus inner layer
E2=12682000000; %Young modulus middle layer
E3=131.7e9; %Young modulus outer layer
v1=0.274; %Poisson's ratio inner layer
v2=0.4; %Poisson's ratio middle layer
v3=0.274; %Poisson's ratio outer layer
a=0.003; %inner radius
b=0.004; %interfase 1 radius
c=0.005; %interfase 2 radius
d=0.006; %Outer radius
P=20e6; %internal pressure

%Axial stresses (constat through layers)
SZ1=double(subs(Sol.SZ1));
SZ2=double(subs(Sol.SZ2));
SZ3=double(subs(Sol.SZ3));

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
%Stress R theta eq layer 3
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

%Plotting solution

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
