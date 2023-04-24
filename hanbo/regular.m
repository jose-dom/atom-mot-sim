function acc=regular(t,r)
    alp=19.48/180*pi;
    Gcoef=2/0.0132^2;
    lambda=780E-9;
    gamma=6.06E6;
    hbar=1.05457182E-34;
    k=1/lambda;
    Is=hbar * (2.42E15)^3 * 2*pi*6.06E6 /(12*pi*9E16);
    I=128.138;
    m=146E-28;
    Aco=I*hbar*k*gamma/(2*m*Is);
    Bco=I/Is;
    
        r1=r(2);
          
        r2=+Aco .* expf(r(3),r(5)) ./ (1 + Bco .* expf(r(3),r(5)) + 4*(kvf(r(2))+1+Zf(r(1)).*Bf(r(1),r(3),r(5))./gamma).^2) ...
           -Aco .* expf(r(3),r(5)) ./ (1 + Bco .* expf(r(3),r(5)) + 4*(kvf(r(2))-1+Zf(r(1)).*Bf(r(1),r(3),r(5))./gamma).^2);
       
        r3=r(4);
          
        r4=+Aco .* expf(r(1),r(5)) ./ (1 + Bco .* expf(r(1),r(5)) + 4*(kvf(r(4))+1+Zf(r(3)).*Bf(r(1),r(3),r(5))./gamma).^2) ...
           -Aco .* expf(r(1),r(5)) ./ (1 + Bco .* expf(r(1),r(5)) + 4*(kvf(r(4))-1+Zf(r(3)).*Bf(r(1),r(3),r(5))./gamma).^2);
       
        r5=r(6);
          
        r6=+Aco .* expf(r(1),r(3)) ./ (1 + Bco .* expf(r(1),r(3)) + 4*(kvf(r(5))+1+Zf(r(5)).*Bf(r(1),r(3),r(5))./gamma).^2) ...
           -Aco .* expf(r(1),r(3)) ./ (1 + Bco .* expf(r(1),r(3)) + 4*(kvf(r(5))-1+Zf(r(5)).*Bf(r(1),r(3),r(5))./gamma).^2);
           
       aof=[r1,r2,r3,r4,r5,r6];
       acc=aof';
       
       
    function exp0=expf(r1,r3)
        exp0=exp(-Gcoef.*(r1.^2 + r3.^2));
    end    

    function kv=kvf(r2)
        kv=(k.*r2)./gamma;
    end
    
    function Bfield=Bf(r1,r3,r5)
        
            Bfield=2.3E8 .* sqrt(r1.^2 + r3.^2 + r5.^2);
    end

    function Zfine=Zf(r5)
        if r5>0
            Zfine=1;
        elseif r5==0
            Zfine=0;
        else 
            Zfine=-1;
        end
    end


end