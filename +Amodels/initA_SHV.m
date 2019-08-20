function [A] = initA(A)
%% USED FOR FAILURE MODELS

%INITA Sets up struct A with constants and containers
% MODEL FLAGS (used to quickly change model)
%VISCOSITY_MODEL_FLAG = 'Hess and Dingwell';
%VISCOSITY_THETA_G_FLAG = 'None';
VISCOSITY_MODEL_FLAG = 'Hess and Dingwell';
VISCOSITY_THETA_G_FLAG = 'None';
VISCOSITY_THETA_C_FLAG = 'Costa';
%VISCOSITY_THETA_C_FLAG = 'None';
CRYSTAL_GROWTH = true;

    A.useForchheimer = true;

    % Set some constants and containers
    A.Patm_ = 1.013e5;     % pascals
    A.vchamber_ = [];
    A.Ptop_ = [];
    A.r = 22.5;             % Conduit radius
    A.depth = 5000;          % Length
    A.Bchm = 1e-10;      % Chamber compressibility (sphere)
    %A.Bchm = 1e-7;       % Chamber compressibility (sill)
    A.Vchm = 3e10;        % Chamber volume
    A.nb = 1e9;        % Bubble concentration
    A.Pchamber = 120e6; % Chamber pressure (Pa)
    A.gamma = 1.29;

    
    % Henry's law constants
    A.hs = 4.1413e-6;
    A.hb = 0.5;
    A.hg = 0.046; %total volatile content
    A.Pcrit = (A.hg/A.hs)^(1/A.hb);   % Pcrit is pressure when volatiles first exsolve
    
    %A.lam = 1-A.hg; % melt mass fraction (1-total volatile mass fraction)
    
    % EOS Values
    %A.rhom0 = 2500;   % kg/m^3 (reference density of melt) (ANDESITE)
    %A.rhom0 = 2700;   % kg/m^3 (reference density of melt) (BASALT)
    %A.rhohc = 741;      % kg/m^3 (density of dissolved volatiles)
    A.rhom0 = 2450; % test vis a vis koyaguchi
    
    A.g = 9.81;         % N/kg (force/mass for gravity)
    
    %A.T = 1200;         % K (temperature of magma) (BASALT)
    A.T = 1123;          % K (temperature of magma) (ANDESITE)
    %A.T = 1000;          % K (temperature of magma) (DACITE)
    %A.T = 900;          % K (temperature of magma) (RHYOLITE)
    
    %A.mu = 100;         % liquid viscosity, Pa s (BASALT)    
    %A.mu = 1e6;         % liquid viscosity, Pa s (ANDESITE)
    A.mu = 1e5;         % liquid viscosity, Pa s (DACITE)
    %A.mu = 1e9;         % liquid viscosity, Pa s (RHYOLITE)
    %A.mu = 1000;         % test vis a vis eric
    
    %A.Rw = 455.59;       % J/kg-K Gas constant for water
    A.Rw = 461.4;
    %A.Plam0 = 10e6;     % reference pressure (for melt density), Pa
    %A.Blam = 1e-10;     % melt compressibility, Pa^-1
    
    %A.rhohc0 = 741;     % reference condensed density, kg/m^3
    %A.Phc0 = 0.1e6;     % refemce pressure (for condensed volatiles), Pa
    %A.Bhc = 2e-10;      % compressibility of dissolved volatiles, Pa^-1
      
    % Kirsch Equations: Conduit Wall Rock 
    A.k.p = .5;         % Poisson's ratio
    A.k.K = A.k.p/(1-A.k.p); % K constant
    A.k.rho = 2500;     % Rock Density
    A.k.P0 = []   ;     % Pore pressure
    
    % Mohr Coulomb Failure
    A.mc.C = 9e6;
    %A.mc.C = 8e6;
    A.mc.phi = deg2rad(15);
    %A.mc.phi = deg2rad(35);
    
    % Fragmentation
    A.f0 = 0.0075;  % Darcy-Weisbach friction factor
    A.phi0 = .80; % critical gas volume fraction for fragmentation
    A.phiforce = .85; % start of transition period (should be less than phi0)
    %A.mug = 1e-5; % gas viscosity
    A.mug = 1.5e-2; % gas viscosity
    A.Rash = 0.001; % ash radius
    A.dragC = 0.8; % drag coefficient
    
    % Lateral gas loss
    A.kw0 = 1e-13; % permability constant (0 - 1e-12) ??? Bounds from kozono/kayaguchi but they have kw = 1e-15 in the paper
    A.kw0 = 0;
    A.Pstar = 20e6; % pressure constant (2.5 - 20 MPa)
    A.rhow = 1e3; % density of water
    A.kc = 2e-19; % Magma permeability
    
    % Forchheimer's Law
    A.ftb = 0.3; % Throat bubble ratio [0.05-0.5]
    A.m = 2.2; % Tortuosity factor
    A.Ff0 = 10;
    
%     %%%%% GRIMSVOTN PARAMS %%%%%
%     A.H = 1700;                % Length for Grimsvotn
%     A.Vchm = 4/3*pi*(3.5e3)^3; % Chamber volume inferred at Grimsvotn
%     A.mu = 100;                % Basaltic viscosity
%     %A.Bchm = 1e-10;           % Chamber compressibility (sphere)
%     %A.Bchm = 1e-7;            % Chamber compressibility (sill)
%     A.Bchm = 1e-9;             % Somewhere in between
     
    A.mu0 = A.mu;
    
    % Henry's law
    w = @(p) min(A.hg, A.hs*p.^A.hb);
    
    xc0 = 0.45; % Crystal content
    xcmax = 0.6;
    
    if CRYSTAL_GROWTH
        xc = @(p) real(min(xcmax, xc0 + 0.55*(0.58815*(p/1e6).^(-0.5226))));
      
    else
        xc = @(p) xc0;
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
            theta_g = @(phi) 1;
    end
    
    switch VISCOSITY_THETA_C_FLAG
        case 'Costa'
            c1 = 0.9995;
            c2 = 0.4;
            c3 = 1;
            B = 2.5;
            
            
            theta_c = @(xc) (1 - c1*erf(sqrt(pi)/2 * xc .* (1 + c2./(1-xc).^c3))).^-(B/c1);
            
            
        otherwise
            theta_c = @(xc) 1;     
    end
    
   
    switch VISCOSITY_MODEL_FLAG
        case 'Hess and Dingwell'
            mufunc = @(w) 10.^(-3.545 + 0.833*log(w) + (9601 - 2368*log(w))./(A.T - (195.7 + 32.25*log(w))));
            A.mu = @(phi,p) mufunc(w(p)*100).*theta_g(phi).*theta_c(xc(p));
        case 'Whittington et al.'
            mufunc = @(w) 10.^(-4.43 + (7618.3 - 17.25*log10(w + 0.26))./(A.T - (406.1 - 292.6*log10(w + 0.26))));
            A.mu =  @(phi,p) mufunc(w(p)).*theta_g(phi).*theta_c(xc(p));
        otherwise
            A.mu = @(phi,p) A.mu0.*theta_g(phi).*theta_c(xc(p));
    end
    
    
    
   
end

