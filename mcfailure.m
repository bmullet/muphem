function [Smax,Sfail,failure] = mcfailure(A,Srr,Szz,Stt,zvec)
    ppore = (A.depth-zvec).*1000*A.g; %pore pressure gradient
    Smin = min([Srr Szz Stt],[],2);
    Smax = max([Srr Szz Stt],[],2); 
    
    qu = 2*A.mc.C*tan(deg2rad(45)+A.mc.phi/2);
    Sfail = qu + (Smin-ppore)*tan(deg2rad(45)+A.mc.phi/2)^2;
    
    if any((Smax-ppore)>Sfail)
       	failure = 1;
    else
        failure = 0;
    end
    
    Smax = Smax - ppore;
       
end