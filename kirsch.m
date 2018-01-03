function [Srr, Szz, Stt] = kirsch (zvec,p,A)
    z = A.depth - zvec; % Change into depth
    %P0 = 1.01e5+10000*z;       % Pore Pressure gradient = 5 kPa/m
    
    Srr = p;    
    Szz = 1.01e5+A.k.rho*(z)*A.g;
    
    S = A.k.K*Szz;     % Poisson's ratio stress condition
    %S = Szz;            % Isotropic stress condition
    Stt = 2*(S)-p;
    
    
end
