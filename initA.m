function [A] = initA(A)
%INITA Sets up struct A with constants and containers

    % Set some constants and containers
    A.Patm_ = 1.01e5;     % pascals
    A.vchamber_ = [];
    A.Ptop_ = [];
    A.r = 30;             % Conduit radius
    A.depth = 4000;          % Length
    A.Bchm = 1e-10;      % Chamber compressibility (sphere)
    %A.Bchm = 1e-7;       % Chamber compressibility (sill)
    A.Vchm = 3e10;        % Chamber volume
    A.nb = 1e15;        % Bubble concentration
    
    % Henry's law constants
    A.hs = 4e-6;
    A.hb = 0.5;
    A.hg = 0.03; %total volatile content
    A.Pcrit = (A.hg/A.hs)^(1/A.hb);   % Pcrit is pressure when volatiles first exsolve
    
    %A.lam = 1-A.hg; % melt mass fraction (1-total volatile mass fraction)
    
    % EOS Values
    %A.rhom0 = 2500;   % kg/m^3 (reference density of melt) (ANDESITE)
    %A.rhom0 = 2700;   % kg/m^3 (reference density of melt) (BASALT)
    %A.rhohc = 741;      % kg/m^3 (density of dissolved volatiles)
    A.rhom0 = 3000;      % test vis a vis koyaguchi
    
    A.g = 9.81;         % N/kg (force/mass for gravity)
    
    %A.T = 1400;         % K (temperature of magma) (BASALT)
    %A.T = 1200;          % K (temperature of magma) (ANDESITE)
    %A.T = 1000;          % K (temperature of magma) (DACITE)
    A.T = 900;          % K (temperature of magma) (RHYOLITE)
    
    %A.mu = 100;         % liquid viscosity, Pa s (BASALT)    
    A.mu = 1e5;         % liquid viscosity, Pa s (ANDESITE)
    %A.mu = 1e6;         % liquid viscosity, Pa s (DACITE)
    %A.mu = 1e8;         % liquid viscosity, Pa s (RHYOLITE)
    %A.mu = 1000;         % test vis a vis eric
    
    A.Rw = 461.5;       % J/kg-K Gas constant for water
    %A.Plam0 = 10e6;     % reference pressure (for melt density), Pa
    %A.Blam = 1e-10;     % melt compressibility, Pa^-1
    
    %A.rhohc0 = 741;     % reference condensed density, kg/m^3
    %A.Phc0 = 0.1e6;     % refemce pressure (for condensed volatiles), Pa
    %A.Bhc = 2e-10;      % compressibility of dissolved volatiles, Pa^-1
      
    % Kirsch Equations: Conduit Wall Rock 
    A.k.p = .5;         % Poisson's ratio
    A.k.K = A.k.p/(1-A.k.p); % K constant
    A.k.rho = 3000;     % Rock Density
    A.k.P0 = []   ;     % Pore pressure
    
    % Mohr Coulomb Failure
    A.mc.C = 9e6;
    %A.mc.C = 8e6;
    A.mc.phi = deg2rad(35);
    %A.mc.phi = deg2rad(35);
    
    % Fragmentation
    A.f0 = 0.01;  % Darcy-Weisbach friction factor
    A.phi0 = .85; % critical gas volume fraction for fragmentation
    A.phiforce = .6; % start of transition period
    A.mug = 1e-5; % gas viscosity
    A.Rash = 0.001; % ash radius
    A.dragC = 0.8; % drag coefficient
    
    % Lateral gas loss
    A.kw0 = 1e-13; % permability constant (0 - 1e-12) ??? Bounds from kozono/kayaguchi but they have kw = 1e-15 in the paper
    A.Pstar = 20e6; % pressure constant (2.5 - 20 MPa)
    A.rhow = 1e3; % density of water
    A.kc = 2e-19; % Magma permeability
    
    %%%%% GRIMSVOTN PARAMS %%%%%
%     A.depth = 6000;                % Length for Grimsvotn
%     A.Vchm = 4/3*pi*(3.5e3)^3; % Chamber volume inferred at Grimsvotn
%     A.mu = 100;                % Basaltic viscosity
%     %A.Bchm = 1e-10;           % Chamber compressibility (sphere)
%     %A.Bchm = 1e-7;            % Chamber compressibility (sill)
%     A.Bchm = 1e-9;             % Somewhere in between
     

end

