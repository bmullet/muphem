function [A] = initA(A)
%% USED FOR FAILURE MODELS


%INITA Sets up struct A with constants and containers
% MODEL FLAGS (used to quickly change model)
%VISCOSITY_MODEL_FLAG = 'Hess and Dingwell';
VISCOSITY_THETA_G_FLAG = 'None';
%VISCOSITY_MODEL_FLAG = 'Hess and Dingwell'; % Rhyolite
VISCOSITY_MODEL_FLAG = 'Whittington et al.'; % Dacite
%VISCOSITY_THETA_G_FLAG = 'Bagdassarove-Dingwell';
VISCOSITY_THETA_C_FLAG = 'Costa';
%VISCOSITY_THETA_C_FLAG = 'None';
FRAGMENTATION_CONDITION = 'phi';
CRYSTAL_GROWTH = false;

if (~exist('A'))
    % set up new struct
    
    A.u0 = 9;
    
    A.useForchheimer = false;

    % Set some constants and containers
    A.lambda = @stressratio;
    
    A.Patm_ = 1.013e5;     % pascals
    A.vchamber_ = [];
    A.Ptop_ = [];
    A.r = 50;             % Conduit radius
    A.depth = 5000;          % Length
    A.chamber_fac = (2*A.lambda(-A.depth) + 1)/3;
    A.Bchm = 1e-10;      % Chamber compressibility (sphere)
    %A.Bchm = 1e-7;       % Chamber compressibility (sill)
    A.Vchm = 3e10;        % Chamber volume
    A.nb = 1e15;        % Bubble concentration

    
    A.gamma = 1.29;
    
    
    % Henry's law constants
    %A.hs = 4.1e-6;
    A.hs = 4.11e-6;
    A.hb = 0.5;
    A.hg = 0.046; %total volatile content
    A.Pcrit = (A.hg/A.hs)^(1/A.hb);   % Pcrit is pressure when volatiles first exsolve
    
    %A.lam = 1-A.hg; % melt mass fraction (1-total volatile mass fraction)
    
    % EOS Values
    %A.rhom0 = 2500;   % kg/m^3 (reference density of melt) (ANDESITE)
    %A.rhom0 = 2700;   % kg/m^3 (reference density of melt) (BASALT)
    %A.rhohc = 741;      % kg/m^3 (density of dissolved volatiles)
    A.rhom0 = 2500; % test vis a vis koyaguchi
    %A.rhom0 = 2600; % test vis a vis koyaguchi
    
    A.g = 9.81;         % N/kg (force/mass for gravity)
    %A.Pchamber = 134814680; % Chamber pressure (Pa)
    %A.Pchamber = 140e6;
    
    
    
    %A.T = 1200;         % K (temperature of magma) (BASALT)
    
    %A.T = 1100;          % K (temperature of magma) (ANDESITE)
    A.T = 1100;          % K (temperature of magma) (DACITE) (1000-1300) (was 1150)
    
    %A.T = 900;          % K (temperature of magma) (RHYOLITE)
    %A.T = 886;
    
    %A.mu = 100;         % liquid viscosity, Pa s (BASALT)    
    %A.mu = 1e4;         % liquid viscosity, Pa s (ANDESITE)
    A.mu = 1e5;         % liquid viscosity, Pa s (DACITE)
    %A.mu = 1e9;         % liquid viscosity, Pa s (RHYOLITE)
    %A.mu = 1000;         % test vis a vis eric
    
    %A.Rw = 455.59;       % J/kg-K Gas constant for water
    A.Rw = 461.4;
    %A.Rw = 322;
    %A.Plam0 = 10e6;     % reference pressure (for melt density), Pa
    %A.Blam = 1e-10;     % melt compressibility, Pa^-1
    
    %A.rhohc0 = 741;     % reference condensed density, kg/m^3
    %A.Phc0 = 0.1e6;     % refemce pressure (for condensed volatiles), Pa
    %A.Bhc = 2e-10;      % compressibility of dissolved volatiles, Pa^-1
      
    % Kirsch Equations: Conduit Wall Rock 
    A.k.p = .5;         % Poisson's ratio
    A.k.K = A.k.p/(1-A.k.p); % K constant
    A.k.rho = 2700;     % Rock Density
    A.k.P0 = []   ;     % Pore pressure
    
    % Mohr Coulomb Failure
    %A.mc.C = 9e6;
    %A.mc.C = 8e6;
    %A.mc.C = 5e6;

    A.mc.C = @cohesion;
    
    %A.mc.phi = deg2rad(15);
    A.mc.phi = deg2rad(35);
    
    % Fragmentation
    A.f0 = 0.0075;  % Darcy-Weisbach friction factor
    A.phi0 = .80; % critical gas volume fraction for fragmentation
    %A.mug = 1e-5; % gas viscosity
    A.mug = 1e-5; % gas viscosity
    A.Rash = 0.001; % ash radius
    A.dragC = 0.8; % drag coefficient
    
    % For strain rate fragmentation
    A.fragcond = FRAGMENTATION_CONDITION;
    A.ezz.k = 0.01;
    A.ezz.Ginf = 3e9;
        
    % Lateral gas loss
    A.kw0 = 1e-13; % permability constant (0 - 1e-12) ??? Bounds from kozono/kayaguchi but they have kw = 1e-15 in the paper
    A.kw0 = 0;
    A.Pstar = 20e6; % pressure constant (2.5 - 20 MPa)
    A.rhow = 1e3; % density of water
    A.kc = 2e-19; % Magma permeability
    
    % Forchheimer's Law
    A.ftb = 0.1; % Throat bubble ratio [0.05-0.5]
    A.m = 3.5; % Tortuosity factor
    A.Ff0 = 10;
    A.rb0 = (A.phi0/(4/3*pi*A.nb))^(1/3);
    
%     %%%%% GRIMSVOTN PARAMS %%%%%
%     A.H = 1700;                % Length for Grimsvotn
%     A.Vchm = 4/3*pi*(3.5e3)^3; % Chamber volume inferred at Grimsvotn
%     A.mu = 100;                % Basaltic viscosity
%     %A.Bchm = 1e-10;           % Chamber compressibility (sphere)
%     %A.Bchm = 1e-7;            % Chamber compressibility (sill)
%     A.Bchm = 1e-9;             % Somewhere in between
     
    A.mu0 = A.mu;
    
    
    xc0 = 0.30; % Crystal content (was 0.30)
    xcmax = 0.30;
    A.xcmax = xcmax;
    A.xc0 = xc0;
    A.xcexp = -0.5226;
    
end
    % now re-build functions
    A.Pchamber = (1.01e5+A.depth*A.g*A.k.rho)*A.chamber_fac; % Chamber pressure (Pa)
    
    
    % Henry's law
    w = @(p) min(A.hg, A.hs*p.^A.hb);

    if CRYSTAL_GROWTH
        A.xc = @(p) min(A.xcmax, A.xc0 + 0.55*(0.58815*(p/1e6).^(-0.5226)));
    else
        A.xc = @(p) A.xc0*ones(size(p));
    end
    
    switch VISCOSITY_THETA_G_FLAG
        case 'Bagdassarove-Dingwell'
            b = 22.4;
            theta_g = @(phi) 1./(1+b*phi);
        case 'Ducamp-Raj'
            b = -3;
            theta_g = @(phi) exp(b*phi/(1-phi));
        case 'Mackenzie'
            theta_g = @(phi) 1 - 5/3*phi;
        case 'Taylor'
            theta_g = @(phi) 1 + phi;
        otherwise
            % Default to no phi dependence
            theta_g = @(phi) ones(size(phi));
    end
    
    switch VISCOSITY_THETA_C_FLAG
        case 'Costa'
            c1 = 0.9995;
            c2 = 0.4;
            c3 = 1;
            B = 2.5;
            theta_c = @(xc) (1 - c1*erf(sqrt(pi)/2 * xc .* (1 + c2./(1-xc).^c3))).^-(B/c1);
            disp('Costa!')
        otherwise
            theta_c = @(xc) ones(size(xc));     
    end
    
   
    switch VISCOSITY_MODEL_FLAG
        case 'Hess and Dingwell'
            mufunc = @(w) 10.^(-3.545 + 0.833*log(w) + (9601 - 2368*log(w))./(A.T - (195.7 + 32.25*log(w))));
            A.mu = @(phi,p) mufunc(w(p)*100).*theta_g(phi).*theta_c(A.xc(p));
        case 'Whittington et al.'
            mufunc = @(w) 10.^(-4.43 + (7618.3 - 17.25*log10(w + 0.26))./(A.T - (406.1 - 292.6*log10(w + 0.26))));
            A.mu =  @(phi,p) mufunc(w(p)*100).*theta_g(phi).*theta_c(A.xc(p));
        otherwise
            A.mu = @(phi,p) A.mu0.*theta_g(phi).*theta_c(A.xc(p));
    end
    
    A.mu0l = A.mu(0,A.Pchamber);
    
    
   
end

function C = cohesion(zvec)
    FOS = 1;
    C_surface = 5e6;
    C_bottom = 20e6;
    C_transition_depth = -3000;

    m = (C_surface - C_bottom)/(0 - C_transition_depth);
    C = (zvec > C_transition_depth).*(C_surface + m*zvec);
    C = C + (zvec <= C_transition_depth)*(C_bottom);
end

function k = stressratio(zvec)
    k_surface = 0.6;
    k_bottom = 0.6;
    k_transition_depth = -2000;

    m = (k_surface - k_bottom)/(0 - k_transition_depth);
    k = (zvec > k_transition_depth).*(k_surface + m*zvec);
    k = k + (zvec <= k_transition_depth)*(k_bottom);

end
