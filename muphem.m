function [ varargout ] = muphem( varargin )
%(MU)lti(PH)ase (E)ruption (M)odel: Solves multiphase system of equations
%for explosive volcanic eruption.
%
%   Solves steady state 2-phase 1-dimensional conduit equations
%   Syntax to call from command window: muphem('multiflow2',{varargin})

  [varargout{1:nargout}] = feval(varargin{:});

end

function [ vargout ] = multiflow2(op,A)
plot = true;

% Perform shooting method via fzero
%c0 = findc0(A); 
% %A.c0 = c0;
vbounds = [3];            % Set upper boundary at 10% speed of sound at critical pressure       
%v_fzero = fzero(@(v) matchPatm(v,A),vbounds,optimset('Display','iter'));
v_fzero = vbounds;
A.v_chamber_i = v_fzero;

% Collect Solution
[zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A);

% Output Solution
if (plot)
    plotmuphem(A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec)
end

% Test for failure
[Srr, Szz, Stt, Srz] = kirsch(zvec,pvec,A,ugvec,umvec,rhogvec,phivec,pvec);
[Smax,Sfail,failure,Sprincipal,sigmavals] = mcfailure(A,Srr,Szz,Stt,Srz,zvec);

% Plot slip directions
slip = plotslipsurfaces(zvec,Sprincipal,A,plot);

if (failure)
    disp('we have a failure!')
else
    disp('no failure')
end

plot_failure(zvec,pvec,phivec,A,Srr,Szz,Stt,Srz, Smax,Sfail)
disp('Pressure at conduit exit:')
disp(min(pvec));
vargout = {A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,failure,Sprincipal,slip,sigmavals};

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
