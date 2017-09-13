function [eosfunc] = eosf(dF)
eosfunc.rhogofp = @rhogofp;
eosfunc.chidofp = @chidofp;
if (dF)
    eosfunc.calcphi = @calcphi;
else
    eosfunc.calcphi = @calcphif;
end
eosfunc.calcvars = @calcvars;

function [ rhog, chi_d, phi ] = calcvars(A,um,p)
rhog = rhogofp(A,p);
chi_d = chidofp(A,p);
if (A.delF)
    phi = calcphi(A,um,chi_d);
else
    phi = calcphif(A,um);
end

function [ rhog ] = rhogofp (A,p)
rhog = p/(A.Rw*A.T);

function [ chi_d ] = chidofp (A,p)
chi_d = A.hs*p.^A.hb;

function [ phi ] = calcphi (A,um,chid)
um0 = A.v_chamber_i;
phi = 1 - (1-A.hg)./(1-chid).*(um0)./um;

function [ phi ] = calcphif (A,um,~) %calculate phi after fragmentation
phi = 1 - (1-A.phi0)*A.umf./um;

% 
% function [rhoha, c, phi] = eos1(A,P,phi)
%     % Henry's Law
%     pcrit = (A.hg/A.hs)^(1/A.hb);
%     hc = min(A.hs*P.^A.hb,A.hg);     % dissolved volatile mass fraction (can't be larger than total volatile mass fraction)
%     
%     % Ideal gas law to find density of water
%     rhog = P./(A.Rw*A.T);
%     drhogdp = 1/(A.Rw*A.T);
%     
%     % Melt density
%     rhol = A.rholam0;
% 
%     % bulk density
%     rho = (A.lam./rhol + ha./rhoha + hc./rhohc).^(-1);
%     
%     % drhodp for bulk compressibility/sound speed calc
%     % note all terms in bulk density are functions of P except for lambda,
%     % hence this is kind of ugly...   
%     drhodp = -(dhadp./rhoha - ha*drhogdp./rhoha.^2 + dhcdp./rhohc - hc.*drhohcdp./rhohc.^2 ...
%         - A.lam.*drholdp./rhol.^2).*(rho).^2;
%     
%     % bulk compressibility
%     beta = 1./rho.*drhodp;
%     
%     % bulk soundspeed
%     c = sqrt(1./drhodp);
%     
%     % exsolve gas mass fraction
%     phi = ha.*(rho./rhoha);
% %     
% % function [rho,phi,c,beta,rhog,rhol,ns,n] = eos1phase(p,A)
% % 
% %   % solubility law
% %   
%   p0 = A.Pcrit; % p>p0 => all gas in solution
%   ns = min(A.s*p.^A.m,A.n0); % dissolved gas mass fraction
%   n = A.n0-ns; % exsolved gas mass fraction
%   dndp = -A.m*A.s.*p.^(A.m-1).*heaviside(p0-p);
%  
%   % liquid
%   
%   rhol = A.rhol0*(1+(p-p0)/A.Kl); % liquid density
%   drholdp = A.rhol0/A.Kl;
%   
%   % mixture
%   
%   rho = 1./(n./rhog+(1-n)./rhol); % mixture density
% 
%   beta = -rho.*(dndp.*(1./rhog-1./rhol)-...
% 		n.*drhogdp./rhog.^2-...
% 		(1-n).*drholdp./rhol.^2); % mixture compressibility
% 
%   % sound speed
%   
%   c = 1./sqrt(rho.*beta);
%   
%   % volume fraction of gas
%   
%   gamma = (1-n)./n.*rhog./rhol;
%   phi = 1./(1+gamma); % gas volume fraction
%   



