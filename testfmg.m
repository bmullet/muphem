function [ Fmg, phi, Fmg1, Fmg2 ] = testfmg( z,y,A )
%TWOPHASEODE ode system of equations for two phase system
eos = eosf(A.delF);
p = y(1); ug = y(2); um = y(3);
[rhog, ~, phi] = eos.calcvars(A,um,p);
rhom = A.rhom0;
Qg = rhog*ug*phi;
Qm = rhom*um*(1-phi);
Qt = A.v_chamber_i*rhom;
% 
% if (abs(1- (Qm + Qg)/Qt) > .01) %throws error if continuity is off by 1% or more
%     disp('[Qm, Qg, Qt]')
%     disp([Qm, Qg, Qt])
%     error('Continuity violated!')
% end

dydz = zeros(3,1);
g = A.g;


if (A.delF)
    %below fragmentation depth, no gas loss
    Fmg = interphase();
    G = gasloss();
    Fmw = meltwallfriction();
    %Fgw = gaswallfriction();
    Fgw = 0; %CORRECT
    delF = 1;
else
    %above fragmentation depth
    Fmg = interphase();
    G = gasloss();
    Fmw = 0; %CORRECT
    %Fmw = meltwallfriction();
    Fgw = gaswallfriction();
    delF = 0;
end

%%% CORRECT EQNS
% Set RHS of equations
dydz(1) = -p*G/Qg;
dydz(2) = -phi*rhog*g + Fmg + G*ug - (1-delF)*Fgw;
dydz(3) = -(1-phi)*rhom*g - Fmg - delF*Fmw;


% %%%% INCORRECT EQNS
% % Set RHS of equations
% dydz(1) = -p*G/Qg;
% dydz(2) = -phi*rhog*g + Fmg + G*ug - phi*Fgw;
% dydz(3) = -(1-phi)*rhom*g - Fmg - (1-phi)*Fmw;
% %%%% INCOORECT

    function G = gasloss()
        % put gas loss function here
        d = A.depth - z; %depth
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
        G=0;
    end

    function Fmg = interphase()
%         if phi == 0
%             Fmg =0;
%         elseif A.delF
%             % Below fragmentation depth
%             rb = (3*phi/(4*pi*A.nb))^(1/3); %bubble radius
%             Fmg = -(9/2*phi*(1-phi)*A.mu*(ug-um)/rb^2);
%             if phi>0.5
%                 t = (phi-0.5)/(A.phi0-0.5);
%                 Fmg2 = -3/8*phi*(1-phi)*A.dragC/A.Rash*rhog*abs(ug-um)*(ug-um);
%                 Fmg = -1*(abs(Fmg))^(1-t)*(abs(Fmg2))^(t)*sign(ug-um);
%             end
%         else
%             % Above fragmentation depth
%             Fmg = -3/8*phi*(1-phi)*A.dragC/A.Rash*rhog*abs(ug-um)*(ug-um);
% 
%         end        

      
            rb = (3*phi/(4*pi*A.nb))^(1/3); %bubble radius
            Fmg1 = -(9/2*phi*(1-phi)*A.mu*(ug-um)/rb^2);
            Fmg2 = -3/8*phi*(1-phi)*A.dragC/A.Rash*rhog*abs(ug-um)*(ug-um);
 
            pf = A.phiforce;
            if phi<pf
                Fmg = Fmg1;
            elseif (phi>=pf) && (phi<A.phi0)
                t = (phi-pf)/(A.phi0-pf);
                Fmg = -1*(abs(Fmg1))^(1-t)*(abs(Fmg2))^(t)*sign(ug-um);
                %Fmg = Fmg1*(1-t) + Fmg2*(t);
            elseif phi>=A.phi0
                Fmg = Fmg2;
            end
    end

    function Fmw = meltwallfriction()
        Fmw = 8*A.mu*um/A.r^2;
        
        
    end

    function Fgw = gaswallfriction()
        % put gas wall friction here
        tau = A.f0.*0.5*rhog.*ug.^2;
        Fgw = 2*tau./A.r;
    end
end

