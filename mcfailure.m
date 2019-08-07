function [Smax,Sfail,failure,S,failstress] = mcfailure(A,Srr,Szz,Stt,Srz,zvec)
    depth = abs(zvec);
    
    ppore = (depth).*1000*A.g; %pore pressure gradient
%     Smin = min([Srr Szz Stt],[],2);
%     Smax = max([Srr Szz Stt],[],2); 
    
    S = zeros(length(zvec),3,3);
    S(:,1,1) = Srr;
    S(:,2,2) = Stt;
    S(:,3,3) = Szz;
    S(:,1,3) = Srz;
    S(:,3,1) = Srz;
    [~,ii] = min(abs(zvec - (A.fragdepth*1.001)));
    disp('Srz at plot point')
    disp(Srz(ii));
    
    principalstress = nan(length(zvec),3);
    % rotate into principal stress directions
    for i=1:length(zvec)
       %principalstress(i,:) = eig(squeeze(S(i,:,:)))';        
        [eV,D] = eig(squeeze(S(i,:,:)));
        [c, ind]=sort(diag(D),'descend'); % store the indices of which columns the sorted eigenvalues come from
        eV=eV(:,ind); % arrange the columns in this order 
        
        % store values
        S(i,:,:) = eV;
        principalstress(i,:) = c;
    end
    
    Smin = min(principalstress,[],2);
    Smax = max(principalstress,[],2);
    
    qu = 2*A.mc.C*tan(pi/4+A.mc.phi/2);
    Sfail = qu + (Smin-ppore)*tan(pi/4+A.mc.phi/2)^2;
    
    if any((Smax-ppore)>Sfail)
       	failure = 1;
    else
        failure = 0;
    end
    
    Smax = Smax - ppore;
    
    failstress = principalstress(ii,:);
end