function [ varargout ] = muphem( varargin )
%(MU)lti(PH)ase (E)ruption (M)odel: Solves multiphase system of equations
%for explosive volcanic eruption.
%
%   Solves steady state 2-phase 1-dimensional conduit equations
%   Syntax to call from command window: muphem('multiflow2',{varargin})

  [varargout{1:nargout}] = feval(varargin{:});

end

function [ vargout ] = multiflow2(op,A,varargin)
plot = true;

if length(varargin) > 0
    params = varargin{1};
    A.depth = params.depth;
    A.xc0 = params.xc;
    A.radius = params.radius;
    A.T = params.T;
    A.Pchamber = A.depth*9.8*2700;
    A = Amodels.initA_Degruyter2012(A); % rebuild functions
end

% Perform shooting method via fzero
%c0 = findc0(A); 
% %A.c0 = c0;
% vbounds = [(sqrt(A.r) - 5.2)/1.6];            % Set upper boundary at 10% speed of sound at critical pressure       
%vbounds = [sqrt(A.r*1.2)/2-2.7]+0.2;
vbounds = A.u0;
options = optimset('Display','iter');
%options = optimset();
v_fzero = fzero(@(v) matchPatm(v,A),vbounds,options);
%v_fzero = vbounds;
A.v_chamber_i = v_fzero;

% Collect Solution
A.fricfac = 2;
zvec = nan;

while isnan(zvec)
    A.fricfac = A.fricfac * 0.5;
    [zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A);    
end

% Output Solution
if (plot)
    plotmuphem(A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec)
end

% Test for failure
[Srr, Szz, Stt, Srz] = kirsch(zvec,pvec,A,ugvec,umvec,rhogvec,phivec,pvec);
[Smax,Sfail,failure,Sprincipal,sigmavals] = mcfailure(A,Srr,Szz,Stt,Srz,zvec);

% Plot balance equations
plotbalances(A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec);

% Plot slip directions
slip = plotslipsurfaces(zvec,Sprincipal,A,plot);
porep = -zvec*1000*9.8;
[shear, no_shear, failure_shear, failure_no_shear] = plotfailureprofiles(A,Srr,Szz,Stt,Srz,zvec,pvec,plot,porep);

if plot
if (failure)
    disp('we have a failure!')
else
    disp('no failure')
end

%plot_failure(zvec,pvec,phivec,A,Srr,Szz,Stt,Srz, Smax,Sfail)
disp('Pressure at conduit exit:')
disp(min(pvec));
end

vargout = {A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,shear,no_shear,failure_shear, failure_no_shear};

end

function [ c0 ] = findc0(A)
% First check that chamber pressure is greater than exsolution pressure

if A.Pchamber < A.Pcrit
    disp('Chamber pressure is')
    disp(A.Pchamber)
    disp('Exsolution pressure is')
    disp(A.Pcrit)
    error('Error: Chamber pressure is lower than exsolution pressure! Check that?')
end
eos = eosf(1);
rhog = eos.rhogofp(A,A.Pcrit);
c0 = ((A.rhom0 - rhog)/(1-A.hg)*A.hs*A.hb*A.Pcrit^(A.hb-1))^(-1/2);
end
