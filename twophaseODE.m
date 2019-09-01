function [ dydz ] = twophaseODE( z,y,A )
%TWOPHASEODE ode system of equations for two phase system
eos = eosf(A.delF);
p = y(1); phi = y(2); du = y(3);
[rhog, ~, um] = eos.calcvars(A,phi,p);
ug = du + um;

dydz = zeros(3,1);
g = A.g;

if (A.delF)
    %below fragmentation depth, no gas loss
    Fmg = interphase();
    G = 0;
    Fmw = meltwallfriction();
    Fgw1 = 0; %CORRECT
    Fgw2 = 0;
    
else
    %above fragmentation depth
    Fmg = interphase();
    G = gasloss();
    Fmw = 0; %CORRECT
    Fgw1 = gaswallfriction1();
    Fgw2 = gaswallfriction2();
  
end

lambda = (1/(um*(1-phi)) + 1/(ug*phi*rhog*A.C.delta));
%Set RHS of equations
%dydz(1) = -G/(phi*rhog);
dydz(1) = 0;
dydz(2) = -1/A.C.Fr^2*(phi*rhog*A.C.delta + (1-phi)) - Fmw - Fgw1;
dydz(3) = -1/A.C.Fr^2*(1/ug-1/um) + Fmw/(1-phi) - Fgw2/phi + lambda*Fmg;

    function G = gasloss()
        % put gas loss function here
        d = abs(z); %depth
        ppore = A.rhow*A.g*d+A.Patm_; % hyrdostatic pore pressure
        %ppore = A.Patm_ + 10000*d;
        %kw = A.kw0*exp(-ppore/A.Pstar);        
        kw = 10.^(-14 - 3.2*log10(d/1000+.1));
        if (A.delF)
            % below fragmentation, add in effect of magma permeability
            km = A.kc*(phi)^3;
            k = 2*(1/kw+1/km)^-1;
        else
            % above fragmentation, no magma permeability effect
            k = kw;
        end
        
        G = 2*rhog*phi*k*(p - ppore)/(A.mug*A.r^2);
        
        G(G<0) = 0; % Set no incoming gas
        G(:) = 0;
    end

    function Fmg = interphase()
     
            rb = (3*phi/(4*pi*A.nb))^(1/3); %bubble radius
            
            if (A.useForchheimer)
                k1 = (phi/A.phi0)^(2/3);
                k2 = (phi/A.phi0)^((1+3*A.m)/2+1/3);
               
                gamma = (1/(um*(1-phi)) + 1/(ug*phi*rhog*A.C.delta));
                Fmg1 = -gamma/A.C.St*(1+A.C.Fo*k1/k2*abs(du)*rhog)*phi*(1-phi)/k1*du;
                Fmg1 = sign(Fmg1)*min(abs(Fmg1),1e5);
                Fmg1 = -du*1e6;
                

            else
                %Fmg1 = -(9/2*phi*(1-phi)*A.mu(phi,p)*(ug-um)/rb^2);
            end
            
            
            Fmg2 = -3/8*phi*(1-phi)*A.dragC*A.C.rc/A.Rash*rhog*abs(du)*(du)*A.C.delta;
          
            pf = A.phiforce;
            
            if (A.delF)
                Fmg = Fmg1;
               
            else
                Fmg = Fmg2;
            end
            
            
%             if phi<A.phi0
%                 Fmg = Fmg1;
%             elseif (phi>=A.phi0) && (phi<pf)
%                 t = (phi-A.phi0)/(pf-A.phi0);
%                 Fmg = -1*(abs(Fmg1))^(1-t)*(abs(Fmg2))^(t)*sign(ug-um);
%             elseif phi>=pf
%                 Fmg = Fmg2;
%             else
%                 Fmg=0;
%             end   
%                    
            

    end

    function Fmw = meltwallfriction()
      
        Fmw = 8*A.mu(phi,p*A.C.p0)/A.C.mu0*um/A.C.Re;
    end

    function Fgw = gaswallfriction1()
        % put gas wall friction here
        Fgw = A.C.delta*A.f0*rhog*ug.^2;
    end

    function Fgw = gaswallfriction2()
        % put gas wall friction here
        Fgw = A.f0*ug/phi;
    end
end

