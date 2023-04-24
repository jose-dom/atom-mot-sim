alp=19.48/180*pi;
Gcoef=2/0.0132^2;
lambda=780E-9;
gamma=6.06E6;
hbar=1.05457182E-34;
k=1/lambda;
Is=hbar * (2.42E15)^3 * 2*pi*6.06E6 /(12*pi*9E16);
I=128.138;
num=100;
kb=1.38E-23;
m=146E-28;
T=0.72;
sigma=sqrt(kb*T/m);

[t,r]=ode15s(@regular,[0,0.3], [0.01, 5, 0.01, 0, 0.003, 5]);

c = m.*(r(:,2).^2 + r(:,4).^2 + r(:,6).^2)./(2*kb); 

figure

subplot(3,3,1) 

plot(t,r(:,1)) 

xlabel("t [s]") 

ylabel("x [m]") 

subplot(3,3,2) 

plot(t,r(:,2)) 

xlabel("t [s]") 

ylabel("v_x [m/s]") 

subplot(3,3,4) 

plot(t,r(:,3)) 

xlabel("t [s]") 

ylabel("y [m]") 

subplot(3,3,5) 

plot(t,r(:,4)) 

xlabel("t [s]") 

ylabel("v_y [m/s]") 

subplot(3,3,7) 

plot(t,r(:,5)) 

xlabel("t [s]") 

ylabel("z [m]") 

subplot(3,3,8) 

plot(t,r(:,6)) 

xlabel("t [s]") 

ylabel("v_z [m/s]") 

subplot(3,3,9) 

plot(t, c )

xlabel("t [s]") 

ylabel("R [m]") 

subplot(3,3,[3 6]) 
 
scatter3(r(:,1), r(:,3), r(:,5),5, c,"filled")

colorbar

xlabel("x [m]") 

ylabel("y [m]") 

zlabel("z [m]") 

%axis([-0.0132 0.0132 -0.0132 0.0132 -0.0132 0.0132 ])