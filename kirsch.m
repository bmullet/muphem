function [Srr, Szz, Stt, Srz] = kirsch (zvec,p,A,ugvec,umvec,rhogvec,phivec,pvec)
    z = abs(zvec); % Change into depth
    %P0 = 1.01e5+10000*z;       % Pore Pressure gradient = 5 kPa/m
    
    Srr = p;    
    Szz = 1.01e5+A.k.rho*(z)*A.g;
    
    S = A.k.K*Szz;     % Poisson's ratio stress condition
    %S = Szz;            % Isotropic stress condition
    Stt = 2*(S)-p;
    
    Srz = nan(size(Stt));
    mu = A.mu(phivec,pvec);
    Srz(zvec<A.fragdepth) = 4*mu(zvec<A.fragdepth).*umvec(zvec<A.fragdepth)/A.r;
    Srz(zvec>=A.fragdepth) = A.f0*rhogvec(zvec>=A.fragdepth).*ugvec(zvec>=A.fragdepth).^2./2;
    Srz = Srz;
    
end
