function [zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A)
%INCOODES: (In)tegrate (Co)nduit (ODEs)
%   Integrates conduit ODEs from base of conduit (z=-A.depth) to top of
%   conduit (z=0)

%% First integrate until we reach p critical (where gas first exsolves)
A.delF = 1; % Turns on/off mass transfer
eos = eosf(A.delF);
pcrit = A.Pcrit*.99; % allow pressure to drop slightly below exsolution so that Qg ~= 0 (overpressure develops)
zspan = [-A.depth 0];
options = odeset('Events',@ExsolutionDepth,'NormControl','on','RelTol',2.5e-14,'AbsTol',1e-17);
y0 = [A.Pchamber];
u0 = A.v_chamber_i;

if A.Pchamber < pcrit
    error('pcrit < A.pchamber!')
end

[z1,y1,ze,ze,ze] = ode45(@(z,y) singlephaseODE(z,y,u0,A),zspan,y0,options);
p1 = y1(:,1); um1=u0*ones(size(z1)); ug1 = um1;
phivec = zeros(size(z1));
rhogvec = zeros(size(z1));
chidvec = zeros(size(z1));
Qmvec = A.rhom0*u0*ones(size(z1));
Qgvec = zeros(size(z1));

%% Now integrate to fragmentation depth using two phase model
delF = 1; % Turns on/off mass transfer
eos = eosf(A.delF);
nz = length(z1);
zstart = z1(nz); A.exdepth = zstart;
zspan = [zstart 0];
p0 = y1(nz); % new p0 = (should be pcrit)
if abs((pcrit-p0)/pcrit) > 0.01 
    disp([pcrit p0])
    error('Pcrit and p0 do not match!!')
end

phi0 = fzero(@(phi) exslvphi(phi,u0,p0), 0); % Find phi0

y0 = [p0 phi0 0]; % format is [p phi delta0];

options = odeset('Events',@FragmentationDepth,'Mass',@mass, 'MStateDependence', 'strong', 'NormControl','on','RelTol',2.5e-9,'AbsTol',1e-9);
sol = ode15s(@(z,y) twophaseODE(z,y,A), zspan, y0, options);

zfrag = sol.x';
pfrag = sol.y(1,:)'; phifrag = sol.y(2,:)'; dufrag = sol.y(3,:)'; 
[ rhogfrag, chidfrag, umfrag] = eos.calcvars(A,phifrag,pfrag);
ugfrag = umfrag + dufrag;

Qmfrag = (1-phifrag).*umfrag.*A.rhom0;
Qgfrag = phifrag.*ugfrag.*rhogfrag;
phivec = [phivec; phifrag];
Qmvec = [Qmvec; Qmfrag];
Qgvec = [Qgvec; Qgfrag];
rhogvec = [rhogvec; rhogfrag];
chidvec = [chidvec; chidfrag];

nz = length(zfrag);
zstart = zfrag(nz);
A.fragdepth = zstart;
zspan = [zstart 0];

A.umf = umfrag(nz); % to be used for new phi calculation




% Do fragmentation depth to surface integration, if needed
if abs(max(zfrag)) <  1
    % We reached the surface with no fragmentation!
    
    A.fragdepth = 0;
    zvec  = [z1; zfrag];
    pvec = [p1; pfrag];
    ugvec = [ug1;  ugfrag];
    umvec = [um1; umfrag];
    
    
else
    % Did not reach surface, so keep going!
    A.delF = 0;
    delF = 0;
    eos = eosf(A.delF);
    
    options = odeset('Events',@RegimeChangeDepth,'Mass',@mass2, 'MStateDependence', 'strong', 'NormControl','on','RelTol',2.5e-9,'AbsTol',1e-9);

    solext = ode15s(@(z,y) twophaseODE(z,y,A), zspan, sol.y(:,end), options);
    
    p2e = solext.y(1,:)'; phi2e = solext.y(2,:)'; du2e = solext.y(3,:)'; 
    z2e = solext.x';
    [ rhog2e, chi_d2e, um2e ] = eos.calcvars(A,phi2e,p2e);
    ug2e = du2e + um2e;

    Qm2e = (1-phi2e).*um2e.*A.rhom0;
    Qg2e = phi2e.*ug2e.*rhog2e;
    phivec = [phivec; phi2e];
    Qmvec = [Qmvec; Qm2e];
    Qgvec = [Qgvec; Qg2e];
    rhogvec = [rhogvec; rhog2e];
    chidvec = [chidvec; chi_d2e];
    
    nz = length(z2e);
    zstart = z2e(nz);
    
    zspan = [zstart 0];

    y0 = [p2e(nz) phi2e(nz) du2e(nz)];
    
    options = odeset('Events',@BlowUp, 'Mass',@mass2, 'MStateDependence', 'strong', 'NormControl','on','RelTol',2.5e-9,'AbsTol',1e-9,'InitialStep',1e-6);
    warning off MATLAB:ode15s:IntegrationTolNotMet
    [z3,y3] = ode15s(@(z,y) twophaseODE(z,y,A), zspan, y0, options);
    
    p3 = y3(:,1); phi3 = y3(:,2); du3 = y3(:,3);
    
    [ rhog3, chi_d3, um3 ] = eos.calcvars(A,phi3,p3);
    ug3 = du3 + um3;
    
    Qm3 = (1-phi3).*um3.*A.rhom0;
    Qg3 = phi3.*ug3.*rhog3;
    phivec = [phivec; phi3];
    Qmvec = [Qmvec; Qm3];
    Qgvec = [Qgvec; Qg3];
    rhogvec = [rhogvec; rhog3];
    chidvec = [chidvec; chi_d3];
    zvec  = [z1; zfrag; z2e; z3];
    pvec = [p1; pfrag; p2e; p3];
    ugvec = [ug1; ugfrag; ug2e; ug3];
    umvec = [um1; umfrag; um2e; um3];
    
end

    function resid = exslvphi(phitest,u0,p)
        [ rhog, ~, u ] = eos.calcvars(A,phitest,p);
        resid = A.rhom0*u0 - (A.rhom0*(1-phitest)*u + rhog*phitest*u); % delta u is zero
    end

    function M = mass(~,y)
        p = y(1); phi = y(2); du = y(3); 
        chid = eos.chidofp(A,p); 
        um = eos.calcum(A,phi,chid);
        ug = du + um;
        rhog = eos.rhogofp(A,p);
        rhom = A.rhom0;
        Qm = (1-phi)*um*rhom; Qg = phi*rhog*ug;
        I = -(Qm)/(1-chid)*A.hs*A.hb*p^(A.hb-1);
        
        alpha = - (Qm+Qg)/(p*(1/ug + (1-phi)/(phi*um)));
        Gamma = (Qm/Qg)*(rhog*ug/(rhom*um) - 1)*(A.hb*chid/(1-chid));
        
        gammat = 1 + alpha*(1-Gamma) - I*du;
        md = (rhom*um^2 - rhog*ug^2)/(um/(1-phi) + ug/phi);
        
        gamma3 = (1/(rhog*ug) - 1/(rhom*um) - 1/p*(Qm/(rhog*phi) - um));
        
        M = zeros(3,3);
        M(1,:) = [(ug/p - (1/(phi*rhog) + 1/((1-phi)*rhom))*I), (ug/phi + um/(1-phi)), 1]; %mass balance
        M(2,:) = [gammat, 0, -md]; % Add momentum balance
        M(3,:) = [gamma3,0, 1]; % Subtract momentum balance    
        
    end

    function M = mass2(~,y)
        p = y(1); phi = y(2); du = y(3); 
        chid = eos.chidofp(A,p); 
        um = eos.calcum(A,phi,chid);
        ug = du + um;
        rhog = eos.rhogofp(A,p);
        rhom = A.rhom0;
        Qm = (1-phi)*um*rhom; Qg = phi*rhog*ug;
        
        alpha = - (Qm+Qg)/(p*(1/ug + (1-phi)/(phi*um)));
        
        gammat = 1 + alpha;
        md = (rhom*um^2 - rhog*ug^2)/(um/(1-phi) + ug/phi);
        
        M = zeros(3,3);
        M(1,:) = [ug/p, (ug/phi + um/(1-phi)), 1]; %mass balance
        M(2,:) = [gammat, 0, -md]; % Add momentum balance
        M(3,:) = [1/(rhog*ug) - 1/(rhom*um),0, 1]; % Subtract momentum balance    
    end

    function [value,isterminal,direction] = ExsolutionDepth(~,y)
        % Exsolution depth function
        value = y(1)-pcrit;     % The value that we want to be zero (p = pcrit)
        isterminal = 1;         % Halt integration
        direction = 0;          % The zero can be approached from either direction
        
    end
     
    function [value,isterminal,direction] = RegimeChangeDepth(~,y)
        % Regime change depth function
        phi = y(2);
        value = A.phiforce-phi;     % The value that we want to be zero (p = pcrit)
        isterminal = 1;         % Halt integration
        direction = 0;          % The zero can be approached from either direction
        
    end

    function [value,isterminal,direction] = FragmentationDepth(~,y)
        % Exsolution depth function
        phi = y(2);
        value = A.phi0-phi;     % The value that we want to be zero (p = pcrit)
        isterminal = 1;         % Halt integration
        direction = 0;          % The zero can be approached from either direction
        
    end

    function [value,isterminal,direction] = BlowUp(~,y)
        % Exsolution depth function
        isterminal = 1;
        value = 1;
        direction = 0;
        if (any(isnan(y)))
            value = 0;
        end
        
    end
end


