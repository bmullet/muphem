function [ resid ] = matchPatm(v,A)
%MATCHPATM Create a residual to help fzero match atmospheric boundary
%condition.
%   - If the v0 leads to a shock in the conduit, the solver will return a
%   solution that is defined only for depths < 0. In this case, penalize with
%   residual proportional to depth of singularity (make model decrease the v0)
%   - If the v0 leads to pressure =! atmospheric, penalize with pressure
%   difference.

A.v_chamber_i = v;
if v<=0
    resid = -20;
    return
end


A.fricfac = 2;
counter = 1;

zvec = nan;

while isnan(zvec)
    A.fricfac = A.fricfac * 0.5;
    [zvec,pvec,~,~,~,~,~,~,~,Ar] = incoodes(A);
    
    if counter > 20
        break
    end
    
    counter = counter + 1;
end

if isnan(zvec)
    resid = -100;
elseif max(zvec) < 0
    % did not make it to surface, must be choked!
    disp('choked!')
    resid = abs(max(zvec))*10;
else
    % made it to surface but pressure is too high
    resid = (A.Patm_-pvec(end))/1e5*40;
end

if imag(resid)>0
    resid = real(resid)*1e6;
    disp('Imaginary residual!!')
end

end