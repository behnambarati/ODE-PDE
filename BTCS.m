%*****************************************
t=4000;       %final time of Run
dt=1;         %time step
%*****************************************
% constants from table 1
cm=0.01;
cp=2500;
k=0.65;
delta= 2;
ro= 370;
dm=2.2e-8;
Tb= 10;        %initial temperature(t=0)
Ti= 60;        %boundary condition(x=0)
T0= 110;       %boundary condition(x=l)
l=0.024;      %length
hlv= 2.5e-9;  
mb= 86;        %initial moisture(t=0)
mi= 45;        %boundary condition(x=0)
m0= 4;         %boundary condition(x=l)
gama= 0;
eps=0.3;
%*****************************************
L=k/(ro*cp);
d=k*dm/(ro*cm*(k+dm*delta*(eps*hlv+gama)));
landa=cp*dm*delta/(cm*(k+dm*delta*(eps*hlv+gama)));
v=cm*(eps*hlv+gama)/cp;
ti=0:dt:t;
n=t/dt;
dx=.001;
h=l/dx;
xi=0:dx:l;
h=h+1;
n=n+1;
nx=24;
T=zeros(h,n);
m=zeros(h,n);
T(1:h,1)=Tb;
m(1:h,1)=mb;
for k=1:h-2
         x(2*k-1)=Tb;
         x(2*k)=mb;
end
teta=1;
aa=L*dt/(dx^2);
bb=d*dt/(dx^2);
T(1,2:n)=Ti;
T(h,2:n)=T0;
m(1,2:n)=mi;
m(h,2:n)=m0;
 for j=2:n
     
A(1,1)=1+2*aa*teta;
A(1,3)=-aa*teta;
A(1,2)=-v;
B(1)=aa*(1-teta)*(x(3)+Ti)+(1-2*aa*(1-teta))*x(1)-v*x(2)+aa*teta*Ti;
%*******************************************************************
A(2,2)=1+2*bb*teta;
A(2,4)=-bb*teta;
A(2,1)=-landa;
B(2)=bb*(1-teta)*(x(4)+mi)+(1-2*bb*(1-teta))*x(2)-landa*x(1)+bb*teta*mi;
%*******************************************************************
A(2*h-5,2*h-5)=1+2*aa*teta;
A(2*h-5,2*h-7)=-aa*teta;
A(2*h-5,2*h-4)=-v;
B(2*h-5)=aa*(1-teta)*(x(2*h-7)+T0)+(1-2*aa*(1-teta))*x(2*h-5)-v*x(2*h-4)+aa*teta*T0;
%*******************************************************************
A(2*h-4,2*h-4)=1+2*bb*teta;
A(2*h-4,2*h-6)=-bb*teta;
A(2*h-4,2*h-5)=-landa;
B(2*h-4)=bb*(1-teta)*(x(2*h-6)+m0)+(1-2*bb*(1-teta))*x(2*h-4)-landa*x(2*h-5)+bb*teta*m0;
%********************************************************************
     for i=3:2:2*h-6
         A(i,i)=1+2*aa*teta;
         A(i,i+2)=-aa*teta;
         A(i,i-2)=-aa*teta;
         A(i,i+1)=-v;
         B(i)=aa*(1-teta)*(x(i+2)+x(i-2))+(1-2*aa*(1-teta))*x(i)-v*x(i+1);
         %*****************************************************************
         A(i+1,i+1)=1+2*bb*teta;
         A(i+1,i+3)=-bb*teta;
         A(i+1,i-1)=-bb*teta;
         A(i+1,i)=-landa;
         B(i+1)=bb*(1-teta)*(x(i+3)+x(i-1))+(1-2*bb*(1-teta))*x(i+1)-landa*x(i);
         
     end
     x=pinv(A)*(B)';
     for k=2:h-1
         T(k,j)=x(2*k-3);
         m(k,j)=x(2*k-2);
       
     end
 end
 %**********************************************************************
 
%figure

subplot(2,2,1)
G=plot(xi,T(:,[40,60,120,240,480,960]))
set(G,'LineWidth',2,'MarkerSize',3)
xlabel('X(m)'),ylabel('Temperature(C)'),title('Temperature Vs Length')

subplot(2,2,2)
H= plot(ti,T(nx/2,:),'y')
set(H,'LineWidth',2,'MarkerSize',3)
xlabel('time(sec)'),ylabel('Temperature(C)'),title('Temperature Vs time')

subplot(2,2,3)
M=plot(ti,m(nx/2,:))
set(M,'LineWidth',2,'MarkerSize',3)
xlabel('time(sec)'),ylabel('Moisture(M)'),title('Moisture Vs time')

subplot(2,2,4)
mesh(ti,xi,T)
xlabel('time(sec)'),ylabel('X(m)'),zlabel('Temperature(C)'),title('Temperature Vs Length and time')

%********************************************************************
%t=90000;
%F= plot(xi,m(:,[100,3000,6000,10000,60000,90000]))
 %set(F,'LineWidth',2,'MarkerSize',3)
 %xlabel('X(m)'),ylabel('Moisture(M)'),title('Moisture Vs Length')

   
   
         
