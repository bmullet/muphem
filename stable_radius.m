%%

A = Amodels.initA_MSH;
A.phi0 = .75;
A.phiforce = 0.80;
A.chamber_fac = 0.70;
A.lambda = 0.50;

A = Amodels.initA_MSH(A);

% do fzero

options = optimset('TolX',0.0005,'Display','iter');

lambdas = 0.6:.01:0.8;
rvec = nan(size(lambdas));

M = 4; % max num of workers

fail = zeros(size(rvec));



parfor (i = 1:length(lambdas),M)
   
   B = A;
   
   B.lambda = lambdas(i);
   
   B.chamber_fac = (1+B.lambda)/2;
   B = Amodels.initA_MSH(B);
   try
       rvec(i) = fzero(@(r) failure(r,B), [100, 160], options);
   catch ME
       disp('FAILED')
       fail(i) = 1;
   end
end

%csvwrite('no_shear.csv',[lambdas; rvec])


function [stressdiff] = failure(r,A)

shear = true;

A.r = r;

out = muphem('multiflow2',0,A);

if shear 
stressdiff = out{11};
f = out{13};

else
    
stressdiff = out{12}; % no shear is 13, shear is 12
f = out{14};

end

if (stressdiff < 0) && f
    % failure but not at fragmentation
    stressdiff = 1;
    
end

end
