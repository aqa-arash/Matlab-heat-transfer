clc
clear%this program was written by Arash Moradian 
disp('This program is written to solve unsteady 1D heat transfer in a rod');
disp('with boundry condition of radiation on one side and constant flux on the other');
choice=menu('Do you want to crate a custom problem or see an example solved?','create custom problem','solve an example'); 
if (choice==1)
t0=input('Initial temperature of the rod: ');
e=input('Emisivity: ');
k=input('Conductivity: ');
a=input('Thermal diffusivity: ');
s=input('Stefan-Boltzmann constant: ');
g0=input('Heat generation: ');
tinf=input('Amient temperature: ');
q0=input('Constant heat flux: ');
l=input('Wall thickness: ');
n=input('Number of nodes: ');
p=input('p = dt/dx^2: ');
nt=input('Number of steps in time: ');
np=input('Interval of ploted steps in time: ');
else
disp('initial temperature has been set to 303K');
t0=273+30;
disp('Emisivity = 0.76 ; Conductivity = 386 ; Thermal diffusivity = 11*10^-5');
e=0.76;
k=386;
a=11*10^-5;
disp('Stefan-Boltzmann constant = 5.668*10^-8 ; Heat generation = 0');
s=5.668*10^-8;
g0=20;
disp('Ambient temperature = 373 ; Constant heat flux = 480');
tinf=273+100;
q0=480;disp('Length of rod = 0.1m ; Number of nodes = 50 ; steps in time = 1500 ');
l=0.1;
n=50;
p=4000;
nt=5000;
np=100;
end
dx=l/(n-1);
dt=p*dx*dx;
x=0:dx:l;
To=zeros(nt,1);
T1=t0.* ones(n,1);
T2=zeros(n,1);
A=zeros(n,n);
b=zeros(n,1);
for i=1:n
    if i==1
        A(1,1)=-(2*p*(a+(e*s*a*dx/k)* (T1(1)^3))-1);
        A(1,2)=2*a*p; b(i)=g0*a*p*dx*dx/k + 2*e*s*a*p*dx*(tinf^4)/k;
    elseif i==n
           A(i,i)=-(2*a*p -1) ;
           A(i,i-1)=2*a*p;
           b(i)=g0*a*p*dx*dx/k - 2*q0*a*p*dx/k;
        else
            A(i,i)=-(2*a*p-1);
            A(i,i-1)=a*p; A(i,i+1)=a*p;
            b(i)=g0*a*p*dx*dx/k;
        end
end
figure(1)
plot(x,T1)
hold on
% Marching in time
for it=1:nt
To(it)=T1(1);
T2=A*T1+b;
T1=T2;
if mod(it,np)==0
    figure(1)
plot(x,T1)
xlabel({'X'});
ylabel({'Temperature'});
hold on
end
A(1,1)=-(2*p*(a+(e*s*a*dx/k)* (T1(1)^3))-1);
end
time=0:1:nt-1;
figure(2)
plot(time,To);
xlabel({'time'});
ylabel({'T at x=0'});