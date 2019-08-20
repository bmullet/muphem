function [ dydz ] = twophaseODE( z,y,A )
%TWOPHASEODE ode system of equations for two phase system
eos = eosf(A.delF);
p = y(1); phi = y(2); du = y(3);
[rhog, ~, um] = eos.calcvars(A,phi,p);
ug = du + um;
rhom = A.rhom0;
Qg = rhog*ug*phi;
Qm = rhom*um*(1-phi);
Qt = A.v_chamber_i*rhom;

dydz = zeros(3,1);
g = A.g;

if (A.delF)
    %below fragmentation depth, no gas loss
    Fmg = interphase();
    G = 0;
    Fmw = meltwallfriction();
    Fgw = 0; %CORRECT
    delF = 1;
else
    %above fragmentation depth
    Fmg = interphase();
    G = gasloss();
    Fmw = 0; %CORRECT
    Fgw = gaswallfriction();
    delF = 0;
end

beta = (2 + Qm/Qg + (1-phi)/phi)/(1/ug + (1-phi)/phi * (1/um));

%Set RHS of equations
dydz(1) = -G/(phi*rhog);
dydz(2) = -(phi * rhog + (1-phi) * rhom)*g + beta*G - (1-delF)*Fgw - delF*Fmw;
dydz(3) = -(1/(ug) - 1/(um))*g + (1/Qm + 1/Qg)*Fmg  + 1/Qg*G - 1/Qg*Fgw + 1/Qm*Fmw;

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
                k1 = (A.ftb*rb)^2/8 * phi^A.m;
                k2 = (A.ftb*rb)/A.Ff0 * phi^(1+3*A.m)/2;
                
                k1 = max(1e-15,k1);
                k2 = max(1e-15,k2);
                     
                Fmg1 = -(A.mug/k1 + rhog/k2*abs(ug-um))*(ug-um)*phi*(1-phi);

            else
                Fmg1 = -(9/2*phi*(1-phi)*A.mu(phi,p)*(ug-um)/rb^2);
            end
            
            
            Fmg2 = -3/8*phi*(1-phi)*A.dragC/A.Rash*rhog*abs(ug-um)*(ug-um);
            pf = A.phiforce;
            
            if (A.delF)
                Fmg = Fmg1;
            else
                Fmg = Fmg2;
            end
            
%             
%             if phi<A.phi0
%                 Fmg = Fmg1;
%             elseif (phi>=A.phi0) && (phi<pf)
%                 t = (phi-A.phi0)/(pf-A.phi0);
%                 Fmg = -1*(abs(Fmg1))^(1-t)*(abs(Fmg2))^(t)*sign(ug-um);
%            
%             elseif phi>=pf
%                 Fmg = Fmg2;
%             else
%                 Fmg=0;
%             end   


            

    end

    function Fmw = meltwallfriction()
      
        Fmw = 8*A.mu(phi,p)*um/A.r^2;
    end

    function Fgw = gaswallfriction()
        % put gas wall friction here
        tau = A.f0.*0.5*rhog.*ug.^2;
        Fgw = 2*tau./A.r;
    end
end

