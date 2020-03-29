%%

A = Amodels.initA_MSH;
A.phi0 = .75;
A.phiforce = 0.80;
A.chamber_fac = 0.70;
A.lambda = 0.50;

A = Amodels.initA_MSH(A);

% do fzero

options = optimset('TolX',0.0005,'Display','iter');

lambdas = 0.5:.01:1;
rvec = nan(size(lambdas));

for i = 1:length(lambdas)
   A.lambda = lambdas(i);
   
   A.chamber_fac = (1+A.lambda)/2;
   A = Amodels.initA_MSH(A);
   
   rvec(i) = fzero(@(r) failure(r,A), [20, 150], options);
    
end


function [stressdiff] = failure(r,A)

A.r = r;

out = muphem('multiflow2',0,A);

stressdiff = out{12};
failure = out{11};

if (stressdiff < 0) && failure
    % failure but not at fragmentation
    stressdiff = 1;
    
end

end