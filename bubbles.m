%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUBBLE GROW
%
% Calculates bubble growth timescales and limits based on formulas as
% presented in Gonnermann and Manga (2012). (See text for details)
%
% Assumes:  
%           - One volatile phase (water)
%           
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set input parameters

% Melt
eta = 1e6;              % (constant) melt viscosity (Pa-S)
pm = 50e6;              % melt pressure (Pa)
dpdt = 80e6/(10*60);    % decompression rate (Pa/s)
T = 1200;               % melt temperature (K)

% Gas
Cw = 0.02;              % total water content (wt %)
S = 3e-6;               % distance between bubbles (m)
R = 1e-6;               % bubble radius (m)
diffp = 10e6;           % bubble overpressure


%% Calculate parameters
Dwr = Cw * exp(-17.14 - 10661/T - 1.772*pm/T/1e6);      %  Rhyolite diffusivity (Zhang and Behrens 2000)
Dwb = Cw * exp(-8.56 - 19910/T);                        %  Basalt diffusivity (Zhang et al 2007)

D = Dwr; % Use basalt

%% Calculate timescales
taudec = pm/dpdt;               % Decompression timescale
taudif = (S-R)^2/D;             % Diffusion timescale
tauvis = eta/diffp;


%% Compare

Pedif = taudif/taudec;

fprintf('-----------------------\n')
fprintf('Term           Value  \n')
fprintf('-----------------------\n')
fprintf('taudec         %.2e   \n', taudec)
fprintf('taudif         %.2e   \n', taudif)
fprintf('tauvis         %.2e   \n', tauvis)
