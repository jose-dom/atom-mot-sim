%Set const%

num=100;
kb=1.38E-23;
m=146E-28;
T=0.1;
sigma=sqrt(kb*T/m);

%M-B distribution in Cartesian Coordinate%

Vx = normrnd(0,sigma,[1,num]) ;
Vy = normrnd(0,sigma,[1,num]) ;
Vz = normrnd(0,sigma,[1,num]) ;
v=sqrt(Vx.^2+Vy.^2+Vz.^2);
figure
histogram(v,25)

%Generate atom position distribution in Cartesian Coordinate%

X=0.01*rand(1,num)-0.005;
Y=0.01*rand(1,num)-0.005;
Z=0.01*rand(1,num)-0.005;

%Gaussian beam generate%

I=128.138; %W/m^2%
Rmot=-0.0132:0.0001:0.0132;
Gcoef=2/(Rmot(end).^2);
GI=I*exp(-Gcoef*Rmot.^2);
%figure
%plot(Rmot,GI)

%Constants setting for Doppler cooling%

lambda=780E-9;
gamma=6.06E6;
hbar=1.05457182E-34;
k=1/lambda;

%Saturation Intensity%

Is=hbar * (2.42E15)^3 * 2*pi*6.06E6 /(12*pi*9E16);
finalplot=[0,0,0,0,0,0];

flagtrap=0;

for i=1:1:num
    
    tbegin=0;
    tend=52E-7;            
    %parameter for test%
%     X(i)=0.01;
%     Vx(i)=10;
%     Y(i)=0.01;
%     Vy(i)=-10;
%     Z(i)=0.01;
%     Vz(i)=5;    
    rplot=[X(i), Vx(i), Y(i), Vy(i), Z(i), Vz(i)];
    
    for j=1:1:4E3

        %Generate spontaneous emissions%
        
        thk=pi*rand(1,1);
        phk=2*pi*rand(1,1);
        hbark=hbar*k*2*pi;
        
        %Run RK4 to solve ODEs% 
        
        [t,r]=ode45(@regular,[tbegin, tend], [X(i), Vx(i), Y(i), Vy(i), Z(i), Vz(i)]);
        
        Vx(i)=r(end,2)+hbark.*sin(thk).*cos(phk)./m;
        Vy(i)=r(end,4)+hbark.*sin(thk).*sin(phk)./m;
        Vz(i)=r(end,6)+hbark.*cos(thk)./m;
        X(i)=r(end,1);
        Y(i)=r(end,3);
        Z(i)=r(end,5);
        tbegin=tbegin+52E-7;
        tend=tend+52E-7;      
        rplot=[rplot;r];
        

%         if sqrt(r(end,2).^2+r(end,6).^2+r(end,4).^2)<1.8904     %1000 uK%
%             flagtrap=flagtrap+1;
%         end

         if sqrt(r(end,2).^2+r(end,6).^2+r(end,4).^2)<0.03     %1000 uK%
            break;
        end
        
    end
    sz(i)=size(rplot,1);
    finalplot=[finalplot;rplot];
end
%Plot the trajectory%

for i=1:1:num
    plotbegin=2+sum(sz(1:i))-sz(i);
    plotend=sum(sz(1:i))+1;
    c= m.*(finalplot(plotbegin:plotend,2).^2 + finalplot(plotbegin:plotend,4).^2 + finalplot(plotbegin:plotend,6).^2)./(2*kb);
    Tend(i)=c(sz(i));
    scatter3(finalplot(plotbegin:plotend,1), finalplot(plotbegin:plotend,3), finalplot(plotbegin:plotend,5),5,c,"filled")
    hold on
    xlabel("x/m",'LineWidth',5)
    ylabel("y/m",'LineWidth',5)
    zlabel("z/m",'LineWidth',5)
    axis([-0.05 0.05 -0.05 0.05 -0.05 0.05 ])
end
colorbar

count=0;
Ttotal=0;

for i=1:1:num
    if Tend(i)<0.001
       Ttotal=Ttotal+Tend(i);
       count=count+1;
    end
end

count
flagtrap
Ttotal/count
