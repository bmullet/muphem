function [zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A)
%INCOODES: (In)tegrate (Co)nduit (ODEs)
%   Integrates conduit ODEs from base of conduit (z=-A.depth) to top of
%   conduit (z=0)
%

% Number of terms on the right hand side of the governing equations (see twophaseODE)
NUM_TERMS_EQ1 = 1;
NUM_TERMS_EQ2 = 3;
NUM_TERMS_EQ3 = 4;

% Integration error tolerances
rtol = 1e-14;
atol = 1e-14;

warning off MATLAB:ode15s:IntegrationTolNotMet

debug = false;

%% First integrate until we reach p critical (where gas first exsolves)
A.delF = 1; % Turns on/off mass transfer
A.endTransition = false; % Flag to mark end of transition phase
eos = eosf(A.delF);
pcrit = A.Pcrit*0.99; % allow pressure to drop slightly below exsolution so that Qg ~= 0 (overpressure develops)
zspan = [-A.depth 0];
options = odeset('Events',@ExsolutionDepth,'NormControl','on','RelTol',2.5e-14,'AbsTol',1e-17);
y0 = [A.Pchamber];
u0 = A.v_chamber_i;

LHS = [];
RHS = [];

if A.Pchamber > pcrit

    [z1,y1,ze,ze,ze] = ode45(@(z,y) singlephaseODE(z,y,u0,A),zspan,y0,options);
    p1 = y1(:,1); um1=u0*ones(size(z1)); ug1 = um1;
    phivec = zeros(size(z1));
    rhogvec = zeros(size(z1));
    chidvec = zeros(size(z1));
    Qmvec = A.rhom0*u0*ones(size(z1));
    Qgvec = zeros(size(z1));
    nz = length(z1);
    zstart = z1(nz)/A.r; A.exdepth = zstart*A.r;
    zspan = [zstart 0];
    p0 = y1(nz)/A.Pchamber; % new p0 = (should be pcrit)
 

else

    phivec = [];
    rhogvec = [];
    chidvec = [];
    Qmvec = [];
    Qgvec = [];
    um1=[];
    ug1=[];
    p1=[];

    zstart = -A.depth/A.r; A.exdepth = zstart;
    zspan = [zstart 0];
    p0 = 1; 
    z1 = [];
    
end
%% Now integrate to fragmentation depth using two phase model
% Non-dimensionalize


C.p0 = A.Pchamber;
C.rc = A.r;
C.c0 = A.hg;
C.rhom = A.rhom0;
C.mu0 = A.mu0l;
C.U0 = sqrt(C.p0/C.rhom);

C.rhog0 = C.p0/(A.Rw*(A.T));
C.Re = C.rc*C.rhom*C.U0/C.mu0;

C.Fr = C.U0/sqrt(C.rc*9.8);
C.k10 = A.phi0^A.m*(A.ftb*A.rb0)^2/8;
C.k20 = (A.ftb*A.rb0)*A.phi0^((1+3*A.m)/2)/A.Ff0;
C.St = A.rhom0*C.k10*C.U0/(A.mug*C.rc);
C.Fo = C.k10*C.U0*C.rhom/(C.k20*A.mug);
C.delta = C.rhog0/C.rhom;

A.C = C;

u0 = u0/C.U0;

delF = 1; % Turns on/off mass transfer
eos = eosf(A.delF);

phi0 = fzero(@(phi) exslvphi(phi,u0,p0), 0); % Find phi0

C.rb0 = (3*phi0/(4*pi*A.nb))^(1/3); %characteristic bubble radius
C.Reb = (C.rb0*C.rhom*C.U0)/C.mu0;

A.C = C;

[ rhogtest, ~, utest ] = eos.calcvars(A,phi0,p0);

if (abs((utest*(phi0*rhogtest*C.rhog0 + (1-phi0)*C.rhom)) - u0*C.rhom)>1e-3)
    error('Mass is not conserved!')
end


y0 = [p0 phi0 0]; % format is [p phi delta0];

% Set integration if doing critical strain rate
if (strcmp(A.fragcond, 'strain rate'))
    A.phi0 = 0.95;
end

warning('');
options = odeset('Events',@FragmentationDepth,'Mass',@mass, 'MStateDependence','strong', 'Stats', 'off', 'NormControl','off','RelTol',rtol*1e4,'AbsTol',[atol, atol, atol]*1e4);
sol = ode15s(@(z,y) twophaseODE(z,y,A, false), zspan, y0, options);
[warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg)
        disp("Warning received")
        disp(warnMsg)
        disp(A.v_chamber_i);
        disp(A.lambda);
        disp(A.r);
    end

% modify output if doing critical strain rate
if (strcmp(A.fragcond, 'strain rate'))
   % calculate strain rate
   zover = sol.x'*C.rc;
   pover = sol.y(1,:)'; phiover = sol.y(2,:)'; duover = sol.y(3,:)'; 
   [ ~, ~, umover] = eos.calcvars(A,phiover,pover);
   ez = DGradient(umover*C.U0,zover);
   
   % calculate critical strain rate
   crit = A.ezz.k*A.ezz.Ginf./A.mu(phiover, pover*C.p0);
   fi = (ez < crit)';
   sol.x = sol.x(fi);
   sol.y = sol.y(:,fi);
   A.fragphi = max(phiover(fi));
   A.phi0 = A.fragphi;
   A.phiforce = A.fragphi + 0.05;
end

zfrag = sol.x'*C.rc;

pfrag = sol.y(1,:)'; phifrag = sol.y(2,:)'; dufrag = sol.y(3,:)'; 
[ rhogfrag, chidfrag, umfrag] = eos.calcvars(A,phifrag,pfrag);

LHS = zeros(length(pfrag), 3, 3);

RHS.eq1 = zeros(length(pfrag),NUM_TERMS_EQ1);
RHS.eq2 = zeros(length(pfrag),NUM_TERMS_EQ2);
RHS.eq3 = zeros(length(pfrag),NUM_TERMS_EQ3);

grads = DGradient(sol.y, sol.x, 2,  '2ndorder');

for i = 1:length(pfrag)
    M = mass(nan, [pfrag(i), phifrag(i), dufrag(i)]);
    LHS(i,:,:) = M.*repmat(grads(:,i)',3,1);
    RHStemp = twophaseODE(nan, [pfrag(i), phifrag(i), dufrag(i)], A, true);
    RHS.eq1(i,:) = RHStemp.RHS1(:);
    RHS.eq2(i,:) = RHStemp.RHS2(:);
    RHS.eq3(i,:) = RHStemp.RHS3(:);
    
end

A.LHS = LHS;
A.RHS = RHS;


ugfrag = umfrag + dufrag;

pfrag = pfrag*C.p0;
rhogfrag = rhogfrag*C.rhog0;
umfrag = umfrag*C.U0;
ugfrag = ugfrag*C.U0;

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
    A.umf = umfrag(end)/C.U0; % to be used for new phi calculation
    A.phi0 = sol.y(2,end);
    A.du0 = sol.y(:,end);
    
    A.delF = 0;
    delF = 0;
    eos = eosf(A.delF);
    
    lastwarn('');
    
    % Set flag for marking we are done with transition
    A.endTransition = true;
    
    A.delF = 0;
    delF = 0;
    eos = eosf(A.delF);
    
    zspan = [zstart 0]/C.rc;
    
    y0 = sol.y(:,end);

    options = odeset('Events',@BlowUp, 'Mass',@mass2, 'MStateDependence', 'strong', 'NormControl','off','RelTol',rtol*1e4,'AbsTol',atol*1e4,'InitialStep',1e-6);
    %warning off MATLAB:ode15s:IntegrationTolNotMet

    sol = ode15s(@(z,y) twophaseODE(z,y,A, false), zspan, y0, options);

    z3 = sol.x'*C.rc;
    p3 = sol.y(1,:)'; phi3 = sol.y(2,:)'; du3 = sol.y(3,:)';
    
    [ rhog3, chi_d3, um3 ] = eos.calcvars(A,phi3,p3);
    
    RHS.eq1 = zeros(length(p3),NUM_TERMS_EQ1);
    RHS.eq2 = zeros(length(p3),NUM_TERMS_EQ2);
    RHS.eq3 = zeros(length(p3),NUM_TERMS_EQ3);
    
    grads = DGradient(sol.y, sol.x, 2,  '2ndorder');

    for i = 1:length(p3)
        M = mass2(nan, [p3(i), phi3(i), du3(i)]);
        LHS(i,:,:) = M.*repmat(grads(:,i)',3,1);
        
        RHStemp = twophaseODE(nan, [p3(i), phi3(i), du3(i)], A, true);
        RHS.eq1(i,:) = RHStemp.RHS1(:);
        RHS.eq2(i,:) = RHStemp.RHS2(:);
        RHS.eq3(i,:) = RHStemp.RHS3(:);
    end
    
    A.LHS = [A.LHS; LHS];
    A.RHS.eq1 = [A.RHS.eq1; RHS.eq1];
    A.RHS.eq2 = [A.RHS.eq2; RHS.eq2];
    A.RHS.eq3 = [A.RHS.eq3; RHS.eq3];
    p3 = p3*C.p0;
    du3 = du3*C.U0;
    rhog3 = rhog3*C.rhog0;
    um3 = um3*C.U0;
    
    ug3 = du3 + um3;
    
    Qm3 = (1-phi3).*um3.*A.rhom0;
    Qg3 = phi3.*ug3.*rhog3;
    phivec = [phivec; phi3];
    Qmvec = [Qmvec; Qm3];
    Qgvec = [Qgvec; Qg3];
    rhogvec = [rhogvec; rhog3];
    chidvec = [chidvec; chi_d3];
    zvec  = [z1; zfrag;  z3];
    pvec = [p1; pfrag; p3];
    ugvec = [ug1; ugfrag;  ug3];
    umvec = [um1; umfrag; um3];
    end
    
if (debug)
       fprintf('Ex depth = %.2f\n', A.exdepth);    
       fprintf('frag depth = %.2f\n', A.fragdepth); 
       
end

    function resid = exslvphi(phitest,u0,p)
        [ rhog, ~, u ] = eos.calcvars(A,phitest,p);
        resid = A.rhom0*u0 - (A.rhom0*(1-phitest)*u + rhog*C.rhog0*phitest*u); % delta u is zero
    end

    function M = mass(~,y)
        p = y(1); phi = y(2); du = y(3); 
        chid = eos.chidofp(A,p); 
        um = eos.calcum(A,phi,chid,p);
        ug = du + um;
        rhog = eos.rhogofp(A,p);
        rhom = A.rhom0;
        beta = eos.calcbeta(A,p);
        

        alpha = -1./p.*((1-phi).*um + phi.*ug.*rhog.*C.delta)./(1./ug + (1-phi)./phi.*1./um);
        
        gam2 = (1./p.*ug + um.*((1-phi).*rhom./(rhog.*phi.*(C.rhog0)) + 1).*(beta));
        
        Gamma = (1-phi)./phi.*(1-um./(ug.*rhog.*C.delta)).*(beta);
        
        gammat = 1 + alpha.*(1-p.*Gamma) + (1-phi).*um.*beta.*du;
        
        md = -(alpha.*p./ug + ug.*rhog.*phi.*C.delta);
        
        qm = (1-phi).*um;
        
        gamma3 = (1./(rhog.*ug).*1./C.delta - 1./(um) - (um.*(1-phi)./(rhog.*phi).*1./C.delta + um).*beta);
        
        M = zeros(3,3);
        M(1,:) = [gam2, (ug/phi + um/(1-phi)), 1]; %mass balance
        M(2,:) = [gammat, 0, -md]; % Add momentum balance
        M(3,:) = qm*[gamma3,0, 1]; % Subtract momentum balance    
        
    end

    function M = mass2(~,y)
        p = y(1); phi = y(2); du = y(3); 
       
        um = eos.calcum(A,phi);
        ug = du + um;
        rhog = eos.rhogofp(A,p);

        alpha = -1/p*((1-phi)*um + phi*ug*rhog*C.delta)/(1/ug + (1-phi)/phi*1/um);
     
        gammat = 1 + alpha;
        
        md = -(alpha*p/ug + ug*rhog*phi*C.delta);
        
        gamma3 = (1/(rhog*ug)*1/C.delta - 1/(um) );
        
        qm = (1-phi)*um;
        
        M = zeros(3,3);
        M(1,:) = [ug/p, (ug/phi + um/(1-phi)), 1]; %mass balance
        M(2,:) = [gammat, 0, -md]; % Add momentum balance
        M(3,:) = qm*[gamma3,0, 1]; % Subtract momentum balance    
        
        
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
        value = (A.phi0-phi)-atol*2;     % The value that we want to be zero (p = pcrit)
        isterminal = 1;         % Halt integration
        direction = 0;          % The zero can be approached from either direction
        
    end

    function [value,isterminal,direction] = BlowUp(~,y)
       
        isterminal = 1;
        value = 1;
        direction = 0;
        if (any(isnan(y)))
            value = 0;
        end
        
    end
end


