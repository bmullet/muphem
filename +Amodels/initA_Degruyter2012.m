function [A] = initA(A)
%% USED FOR FAILURE MODELS

%INITA Sets up struct A with constants and containers
% MODEL FLAGS (used to quickly change model)
VISCOSITY_MODEL_FLAG = 'Hess and Dingwell';
%VISCOSITY_THETA_G_FLAG = 'Bagdassarove-Dingwell';
%VISCOSITY_MODEL_FLAG = 'None';
VISCOSITY_THETA_G_FLAG = 'None';

    % Set some constants and containers
    A.Patm_ = 1.013e5;     % pascals
    A.vchamber_ = [];
    A.Ptop_ = [];
    A.r = 30;             % Conduit radius
    A.depth = 4000;          % Length
    A.Bchm = 1e-10;      % Chamber compressibility (sphere)
    %A.Bchm = 1e-7;       % Chamber compressibility (sill)
    A.Vchm = 3e10;        % Chamber volume
    A.nb = 1e15;        % Bubble concentration
    A.Pchamber = 140000000; % Chamber pressure (Pa)
    A.gamma = 4/3;
    
    % Henry's law constants
    A.hs = 4.109999e-6;
    A.hb = 0.5;
    A.hg = 0.035; %total volatile content
    A.Pcrit = (A.hg/A.hs)^(1/A.hb);   % Pcrit is pressure when volatiles first exsolve
    
    %A.lam = 1-A.hg; % melt mass fraction (1-total volatile mass fraction)
    
    % EOS Values
    %A.rhom0 = 2500;   % kg/m^3 (reference density of melt) (ANDESITE)
    %A.rhom0 = 2700;   % kg/m^3 (reference density of melt) (BASALT)
    %A.rhohc = 741;      % kg/m^3 (density of dissolved volatiles)
    A.rhom0 = 3000; % test vis a vis koyaguchi
    
    A.g = 9.81;         % N/kg (force/mass for gravity)
    
    %A.T = 1200;         % K (temperature of magma) (BASALT)
    A.T = 1159;          % K (temperature of magma) (ANDESITE)
    %A.T = 1000;          % K (temperature of magma) (DACITE)
    %A.T = 900;          % K (temperature of magma) (RHYOLITE)
    
    %A.mu = 100;         % liquid viscosity, Pa s (BASALT)    
    %A.mu = 1e6;         % liquid viscosity, Pa s (ANDESITE)
    A.mu = 1e7;         % liquid viscosity, Pa s (DACITE)
    %A.mu = 1e9;         % liquid viscosity, Pa s (RHYOLITE)
    %A.mu = 1000;         % test vis a vis eric
    
    A.Rw = 461.5;       % J/kg-K Gas constant for water
    %A.Plam0 = 10e6;     % reference pressure (for melt density), Pa
    %A.Blam = 1e-10;     % melt compressibility, Pa^-1
    
    %A.rhohc0 = 741;     % reference condensed density, kg/m^3
    %A.Phc0 = 0.1e6;     % refemce pressure (for condensed volatiles), Pa
    %A.Bhc = 2e-10;      % compressibility of dissolved volatiles, Pa^-1
      
    % Kirsch Equations: Conduit Wall Rock 
    A.k.p = .3;         % Poisson's ratio
    A.k.K = A.k.p/(1-A.k.p); % K constant
    A.k.rho = 3000;     % Rock Density
    A.k.P0 = []   ;     % Pore pressure
    
    % Mohr Coulomb Failure
    A.mc.C = 9e6;
    %A.mc.C = 8e6;
    A.mc.phi = deg2rad(15);
    %A.mc.phi = deg2rad(35);
    
    % Fragmentation
    A.f0 = 0.01;  % Darcy-Weisbach friction factor
    A.phi0 = .8; % critical gas volume fraction for fragmentation
    A.phiforce = .85; % end of transition period
    %A.mug = 1e-5; % gas viscosity
    A.mug = 1.5e-2; % gas viscosity
    A.Rash = 0.01; % ash radius
    A.dragC = 0.8; % drag coefficient
    
    % Lateral gas loss
    A.kw0 = 1e-13; % permability constant (0 - 1e-12) ??? Bounds from kozono/kayaguchi but they have kw = 1e-15 in the paper
    A.Pstar = 20e6; % pressure constant (2.5 - 20 MPa)
    A.rhow = 1e3; % density of water
    A.kc = 2e-19; % Magma permeability
    
%     %%%%% GRIMSVOTN PARAMS %%%%%
%     A.H = 1700;                % Length for Grimsvotn
%     A.Vchm = 4/3*pi*(3.5e3)^3; % Chamber volume inferred at Grimsvotn
%     A.mu = 100;                % Basaltic viscosity
%     %A.Bchm = 1e-10;           % Chamber compressibility (sphere)
%     %A.Bchm = 1e-7;            % Chamber compressibility (sill)
%     A.Bchm = 1e-9;             % Somewhere in between
     
    A.mu0 = A.mu;
    
    switch VISCOSITY_THETA_G_FLAG
        case 'Bagdassarove-Dingwell'
            b = 22.4;
            theta_g = @(phi) 1/(1+b*phi);
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
   
    switch VISCOSITY_MODEL_FLAG
        case 'Hess and Dingwell'
            mufunc = @(w) 10^(-3.545 + 0.833*log(w) + (9601 - 2368*log(w))/(A.T - (195.7 + 32.25*log(w))));
            A.mu = @(phi,p) mufunc(w(p))*theta_g(phi);
        otherwise
            A.mu = @(phi,p) A.mu0*theta_g(phi);
    end
    
end

