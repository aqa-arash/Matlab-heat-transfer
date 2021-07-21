clc
clear
%this program was written by Arash Moradian 
disp('This program is written to solve 2D heat transfer on a square');
choice=menu('Do you want to crate a custom problem or see an example solved?','create custom problem','solve an example'); 
if (choice==1)
disp('Boundry condition for all walls has been set to convective heat transfer');
X= input('Length of axies: ');
N=input('Number of nodes along each axies: ');
G0=input('Heat generation inside the square: ');
hu=input('Convective heat coefficient on upper wall: ');
hd=input('Convective heat coefficient on bottom wall: ');
hr=input('Convective heat coefficient on right wall: ');
hl=input('Convective heat coefficient on left wall: ');
Tu=input('Ambient temperature on upper side: ');
Td=input('Ambient temperature on lower side: ');
Tr=input('Ambient temperature on right side: ');
Tl=input('Ambient temperature on left side: ');
kc=input('Conductivity: ');
else
    disp('considering a square with following conditions')
    disp('axies lenght is 8cm');
    X= 8;
    disp('number of nodes along each axies is 80');
N=80;
disp('considering no heat generation')
G0=0;
disp('upper and lower wall have convective heat coefficient of 1400 W/m^2 K')
hu=1400;
hd=1400;
disp('ambient temperature is 300K')
Tu=300;
Td=300;
disp('right wall is isolated')
hr=0;
Tr=1;
disp('left wall has a constant temperature of 423K');
hl=hu*10;
Tl=423;
disp('and Conductivity is 20 W/mK')
kc=20;
end
k=0;
NN=N*N;
dx=X/(N-1);
A=zeros(NN,NN);
disp('... calculating coefficients matrix ...')
for j=1:N
    for i=1:N
        k=k+1;
        if (i==1 & j==1)
A(k,k)=-2-(hd+hl)*dx/kc; A(k,k+1)=1;  A(k,k+N)=1; B(k)=-(dx*(hl*Tl+hd*Td)/kc) -(G0*dx*dx/(2*kc));
        elseif (j==1 & i>=2 & i<=(N-1))
A(k,k)=-2*(2+hd*dx/kc); A(k,k+1)=1; A(k,k-1)=1; A(k,k+N)=2; B(k)=-(2*hd*dx*Td/kc) -(G0*dx*dx/kc);
 elseif (i==N & j==1)
A(k,k)=-2-(hd+hr)*dx/kc; A(k,k-1)=1;  A(k,k+N)=1; B(k)=-(dx*(hr*Tr+hd*Td)/kc) -(G0*dx*dx/(2*kc));
        elseif (i==1 & j>=2 & j<=(N-1))
A(k,k)=-2*(2+hl*dx/kc); A(k,k+1)=2; A(k,k-N)=1; A(k,k+N)=1;
B(k)=-(2*hl*dx*Tl/kc) -(G0*dx*dx/kc);
        elseif (i==N & j>=2 & j<=(N-1))
A(k,k)=-2*(2+hr*dx/kc); A(k,k+N)=1; A(k,k-N)=1; A(k,k-1)=2; B(k)=-(2*hr*dx*Tr/kc) -(G0*dx*dx/kc);
 elseif (i==1 & j==N)
A(k,k)=-2-(hu+hl)*dx/kc; A(k,k+1)=1;  A(k,k-N)=1; B(k)=-(dx*(hl*Tl+hu*Tu)/kc) -(G0*dx*dx/(2*kc));
 elseif (j==N & i>=2 & i<=(N-1))
A(k,k)=-2*(2+hu*dx/kc); A(k,k+1)=1; A(k,k-1)=1; A(k,k-N)=2; B(k)=-(2*hu*dx*Tu/kc) -(G0*dx*dx/kc);
 elseif (i==N & j==N)
A(k,k)=-2 -(hu+hr)*dx/kc; A(k,k-1)=1;  A(k,k-N)=1; B(k)=-(dx*(hr*Tr+hu*Tu)/kc) -(G0*dx*dx/(2*kc));
 else
A(k,k)=-4; A(k,k+1)=1; A(k,k-1)=1; A(k,k-N)=1; A(k,k+N)=1;
B(k)=-G0*dx*dx/kc;
end
    end
end
disp('...Done ... now inversing coefficients matrix (this might take a while) ...')
A=inv(A);
disp('... Done ... now calculating the results ...')
tt=A*B';
kk=0;
for i=1:N
    for j=1:N
        kk=kk+1; T(i,j)=tt(kk);
    end
end
x=0:dx:X; y=x;
disp('... Done ... Ploting contour.')
cs=contourf(x,y,T,20);