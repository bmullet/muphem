% %%
% 
% A.lambda = 0.6;
% A = Amodels.initA_MSH(A);
% r = fzero(@(r) failure(r,A), [50, 180], options);

%%

A = Amodels.initA_MSH;
A.phi0 = .70;
A.phiforce = 0.75;
A.chamber_fac = 0.80;
A.lambda = 0.50;

A = Amodels.initA_MSH(A);

% do fzero

options = optimset('TolX',0.0005,'Display','iter');

%lambdas = 0.5:.01:1;
lambdas = [0.5:.01:1];
rvec = nan(size(lambdas));

M = 6; % max num of workers

fail = zeros(size(rvec));




%%


parfor (i = 1:length(lambdas),M)
%for i = 1:length(lambdas)
     
   B = A;
   
   B.lambda = lambdas(i);
   
   %B.chamber_fac = (1+B.lambda)/2;
   B = Amodels.initA_MSH(B);
   
   if lambdas(i) < 0.58
    x = 71.58;
   else
    x = 84.53 + (lambdas(i) - 0.61)*471.40
   end
   
   try
       rvec(i) = fzero(@(r) failure(r,B), x, options);
       disp('Found one!')
       disp(rvec(i))
   catch ME
       disp('FAILED')
       fail(i) = 1;
   end
end

%csvwrite('no_shear.csv',[lambdas; rvec])


function [stressdiff] = failure(r,A)

shear = false;

A.r = r;

out = muphem('multiflow2',0,A);

if shear 
stressdiff = out{11};
f = out{13};

else
    
stressdiff = out{12}; % no shear is 12, shear is 11
f = out{14};

end

if (stressdiff < 0) && f
    % failure but not at fragmentation
    stressdiff = 1;
    
end

end
