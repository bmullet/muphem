% %%
% 
% A.lambda = 0.6;
% A = Amodels.initA_MSH(A);
% r = fzero(@(r) failure(r,A), [50, 180], options);

%% test

%sd = failure(15, B);

%%

A = Amodels.initA_paper_dacite;
A.fragcond = 'phi';
%A.phi0 = .65;
%A.chamber_fac = 0.5;

% do fzero

options = optimset('TolX',0.0005,'Display','iter');

lambdas = [.6:.01:2];

rvec = nan(size(lambdas));

M = 6; % max num of workers

fail = zeros(size(rvec));

%%

parfor (i = 1:length(lambdas),M)
%for i = 1:length(lambdas)
     
   B = A;
   
   B.lambda = @(x) ones(size(x))*lambdas(i);
   B.chamber_fac = (2*lambdas(i) + 1)/3;
   B = Amodels.initA_paper_dacite(B);
   
   %x = lambdas(i)*30 + 10;
   %x = 15;
   x = -lambdas(i)*15 + 45
   try
       rvec(i) = fzero(@(r) failure(r,B), x, options);
       disp('Found one!')
       disp(rvec(i))
   catch ME
      disp('FAILED')
      disp(ME)
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
