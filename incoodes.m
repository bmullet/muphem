function [zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A)
%INCOODES: (In)tegrate (Co)nduit (ODEs)
%   Integrates conduit ODEs from base of conduit (z=0) to top of
%   conduit (z=A.depth)
zprint =[];
Qmprint =[];
Qgprint = [];
% First integrate until we reach p critical (where gas first exsolves)
A.delF = 1; % Turns on/off mass transfer
eos = eosf(A.delF);
pcrit = A.Pcrit*.99; % allow pressure to drop slightly below exsolution so that Qg ~= 0 (overpressure develops)
zspan = [0 A.depth];
options = odeset('Events',@ExsolutionDepth,'NormControl','on','RelTol',2.5e-14,'AbsTol',1e-17);
y0 = [A.Pchamber];
u0 = A.v_chamber_i;

[z1,y1,ze,ze,ze] = ode45(@(z,y) singlephaseODE(z,y,u0,A),zspan,y0,options);
p1 = y1(:,1); um1=u0*ones(size(z1)); ug1 = um1;
phivec = zeros(size(z1));
rhogvec = zeros(size(z1));
chidvec = zeros(size(z1));
Qmvec = A.rhom0*u0*ones(size(z1));
Qgvec = zeros(size(z1));

% Now integrate to fragmentation depth using two phase model
delF = 1; % Turns on/off mass transfer
eos = eosf(A.delF);
nz = length(z1);
zstart = z1(nz); A.exdepth = zstart;
zspan = [zstart A.depth];
p0 = y1(nz); % new p0 = (should be pcrit)
if abs((pcrit-p0)/pcrit) > 0.01 
    disp([pcrit p0])
    error('Pcrit and p0 do not match!!')
end
%[ rhog0, ~, phi0 ] = eos.calcvars(A,u0,p0);
% ug0 = u0;
% um0 = (A.rhom0*u0 - phi0*rhog0*u0)/((1-phi0)*A.rhom0); % to conserve continuity melt velocity maybe needs to drop
ug0 = fzero(@(u) exslvv(u,u0,p0),u0); % Find v0 for exsolve phase

y0 = [y1(nz) ug0 ug0]; % format is [p ug um];
%warning on MATLAB:ode15s:IntegrationTolNotMet
options = odeset('Events',@RegimeChangeDepth,'Mass',@mass, 'MStateDependence', 'strong', 'NormControl','on','RelTol',2.5e-9,'AbsTol',1e-10);
sol = ode15s(@(z,y) twophaseODE(z,y,A), zspan, y0, options);

%p2 = y2(:,1);ug2 = y2(:,2); um2 = y2(:,3);
% p2 = sol.y(1,:)'; ug2 = sol.y(2,:)'; um2 = sol.y(3,:)'; 
% z2 = sol.x';
% [ rhog2, chi_d2, phi2 ] = eos.calcvars(A,um2,p2);
% Qm2 = (1-phi2).*um2.*A.rhom0;
% Qg2 = phi2.*ug2.*rhog2;
% phivec = [phivec; phi2];
% Qmvec = [Qmvec; Qm2];
% Qgvec = [Qgvec; Qg2];
% rhogvec = [rhogvec; rhog2];
% chidvec = [chidvec; chi_d2];

options = odeset('Events',@FragmentationDepth,'Mass',@mass, 'MStateDependence', 'strong', 'NormControl','on','RelTol',2.5e-9,'AbsTol',1e-10);

solext = odextend(sol,[],A.depth,sol.y(:,end),options);
p2e = solext.y(1,:)'; ug2e = solext.y(2,:)'; um2e = solext.y(3,:)'; 
z2e = solext.x';
[ rhog2e, chi_d2e, phi2e ] = eos.calcvars(A,um2e,p2e);
Qm2e = (1-phi2e).*um2e.*A.rhom0;
Qg2e = phi2e.*ug2e.*rhog2e;
phivec = [phivec; phi2e];
Qmvec = [Qmvec; Qm2e];
Qgvec = [Qgvec; Qg2e];
rhogvec = [rhogvec; rhog2e];
chidvec = [chidvec; chi_d2e];

if max(z2e) == A.depth
    % We reached the surface with no fragmentation!
    
%     A.fragdepth = 10000;
%     zvec  = [z1; z2; z2e];
%     pvec = [p1; p2; p2e];
%     ugvec = [ug1; ug2; ug2e];
%     umvec = [um1; um2; um2e];
    
    A.fragdepth = 0;
    zvec  = [z1; z2e];
    pvec = [p1; p2e];
    ugvec = [ug1;  ug2e];
    umvec = [um1; um2e];
    
    
else
    %delF = 0; % Turns on/off mass transfer
    A.delF = 0;
    delF = 0;
    eos = eosf(A.delF);
    nz = length(z2e);
    zstart = z2e(nz);
    A.fragdepth = zstart;
    zspan = [zstart A.depth];
    y0 = [p2e(nz) ug2e(nz) um2e(nz)];
    
    A.umf = um2e(nz); % to be used for new phi calculation
    options = odeset('Events',@BlowUp, 'Mass',@mass2, 'MStateDependence', 'strong', 'NormControl','on','RelTol',2.5e-14,'AbsTol',1e-10,'InitialStep',1e-6);
    warning off MATLAB:ode15s:IntegrationTolNotMet
    [z3,y3] = ode15s(@(z,y) twophaseODE(z,y,A), zspan, y0, options);
  
    p3 = y3(:,1); ug3 = y3(:,2); um3 = y3(:,3);
    [ rhog3, chi_d3, phi3 ] = eos.calcvars(A,um3,p3);
    Qm3 = (1-phi3).*um3.*A.rhom0;
    Qg3 = phi3.*ug3.*rhog3;
    phivec = [phivec; phi3];
    Qmvec = [Qmvec; Qm3];
    Qgvec = [Qgvec; Qg3];
    rhogvec = [rhogvec; rhog3];
    chidvec = [chidvec; chi_d3];
    zvec  = [z1; z2e; z3];
    pvec = [p1; p2e; p3];
    ugvec = [ug1; ug2e; ug3];
    umvec = [um1; um2e; um3];
    
end

    function resid = exslvv(utest,u0,p)
        [ rhog, ~, phi ] = eos.calcvars(A,utest,p);
        resid = A.rhom0*u0 - (A.rhom0*(1-phi)*utest + rhog*phi*utest);
    end

    function M = mass(z,y)
        p = y(1); ug = y(2); um = y(3); 
        chid = eos.chidofp(A,p); 
        phi = eos.calcphi(A,um,chid);
        rhog = eos.rhogofp(A,p);
        rhom = A.rhom0;
        Qm = (1-phi)*um*rhom; Qg = phi*rhog*ug;
        I = -(Qm)/(1-chid)*A.hs*A.hb*p^(A.hb-1);
        
        M = zeros(3,3);
%         M(1,:) = [1+I*p/Qg*(rhog*ug/(rhom*um)-1), p/ug, (1-phi)/phi*p/um]; %dpdz equation
%         M(2,:) = [phi+I*ug, Qg, 0]; % Gas momentum balance
%         M(3,:) = [1-phi-I*um,0, Qm]; % Melt momentum balance    
        M(1,:) = [1+I*p/Qg*(rhog*ug/(rhom*um)-1), p/ug, (1-phi)/phi*p/um]; %dpdz equation
        M(2,:) = [phi, 0, 0]; % Gas momentum balance
        M(3,:) = [1-phi,0, Qm]; % Melt momentum balance  
        zprint =[zprint; z];
        Qmprint =[Qmprint; Qm];
        Qgprint = [Qgprint; Qg];
    end

    function M = mass2(~,y)
        p = y(1); ug = y(2); um = y(3); 
        chid = eos.chidofp(A,p); 
        phi = eos.calcphi(A,um,chid);
        rhog = eos.rhogofp(A,p);
        rhom = A.rhom0;
        Qm = (1-phi)*um*rhom; Qg = phi*rhog*ug;
        M = zeros(3,3);
        M(1,:) = [1, p/ug, (1-phi)/phi*p/um]; %dpdz equation
        M(2,:) = [phi, Qg, 0]; % Gas momentum balance
        M(3,:) = [1-phi,0, Qm]; % Melt momentum balance    
    end
    function [value,isterminal,direction] = ExsolutionDepth(~,y)
        % Exsolution depth function
        value = y(1)-pcrit;     % The value that we want to be zero (p = pcrit)
        isterminal = 1;         % Halt integration
        direction = 0;          % The zero can be approached from either direction
        
    end
     
    function [value,isterminal,direction] = RegimeChangeDepth(~,y)
        % Regime change depth function
        p = y(1); ug = y(2); um = y(3); 
        chid = eos.chidofp(A,p); 
        phi = eos.calcphi(A,um,chid);
        value = A.phiforce-phi;     % The value that we want to be zero (p = pcrit)
        isterminal = 1;         % Halt integration
        direction = 0;          % The zero can be approached from either direction
        
    end

    function [value,isterminal,direction] = FragmentationDepth(~,y)
        % Exsolution depth function
        p = y(1); ug = y(2); um = y(3); 
        chid = eos.chidofp(A,p); 
        phi = eos.calcphi(A,um,chid);
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


