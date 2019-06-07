function [Smax,Sfail,failure] = mcfailure(A,Srr,Szz,Stt,Srz,zvec)
    ppore = (A.depth-zvec).*1000*A.g; %pore pressure gradient
%     Smin = min([Srr Szz Stt],[],2);
%     Smax = max([Srr Szz Stt],[],2); 
    
    S = zeros(length(zvec),3,3);
    S(:,1,1) = Srr;
    S(:,2,2) = Stt;
    S(:,3,3) = Szz;
    S(:,1,3) = Srz;
    S(:,3,1) = Srz;
    
    principalstress = nan(length(zvec),3);
    % rotate into principal stress directions
    for i=1:length(zvec)
       principalstress(i,:) = eig(squeeze(S(i,:,:)))';        
    end
    
    Smin = min(principalstress,[],2);
    Smax = max(principalstress,[],2);
    
    qu = 2*A.mc.C*tan(deg2rad(45)+A.mc.phi/2);
    Sfail = qu + (Smin-ppore)*tan(deg2rad(45)+A.mc.phi/2)^2;
    
    if any((Smax-ppore)>Sfail)
       	failure = 1;
    else
        failure = 0;
    end
    
    Smax = Smax - ppore;
     
end