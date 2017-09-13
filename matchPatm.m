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
    
    
     [zvec,pvec,~,~,~,~,~,~,~,A] = incoodes(A);
    if max(zvec) < A.depth
        % did not make it to surface
        resid = A.depth-max(zvec);
    else
        % made it to surface but pressure is too high
        resid = (A.Patm_-pvec(end))/1e5;      
    end
    
    if imag(resid)>0
        resid = real(resid)*1e6;
        disp('Imaginary residual!!')
    end
    
end