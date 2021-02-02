% %%
% 
% A.lambda = 0.6;
% A = Amodels.initA_MSH(A);
% r = fzero(@(r) failure(r,A), [50, 180], options);

%% test

%sd = failure(15, B);

% Doing the maximum or minimum radius test
max_or_min = "min";

% 
vary = "k"; 


A = Amodels.initA_paper_dacite;
A.fragcond = 'phi';
A.phi0 = 0.75;
%A.phi0 = .65;
%A.chamber_fac = 0.5;

% do fzero

options = optimset('TolX',0.0005,'Display','iter');

lambdas = [0.5:.1:0.7];

rvec = nan(size(lambdas));
rz = nan(size(lambdas));
zr = nan(size(lambdas));
rt = nan(size(lambdas));
tr = nan(size(lambdas));
tz = nan(size(lambdas));
zt = nan(size(lambdas));
outs = cell(size(lambdas));


M = 6; % max num of workers

fail = zeros(size(rvec));

%%

parfor (i = 1:length(lambdas),M)
%for i = 1:length(lambdas)
     
   B = A;
   if vary == "k"
     B.lambda = @(x) ones(size(x))*lambdas(i);
     B.chamber_fac = (2*lambdas(i) + 1)/3;
     x = -lambdas(i)*(33) + 79.8; % with shear minimum
   else
     B.chamber_fac = lambdas(i);
   end
     B = Amodels.initA_paper_dacite(B);
   
   %x = lambdas(i)*50 + 10;
   %x = 15;
   %x = -lambdas(i)*15 + 45; % minimum
   %x = 30;
   %x = -lambdas(i)*(40) + 90; % no shear minimum
   
   %x = (lambdas(i) - 2)*(-96.42) + 15 % maximum
   try

   
       [r, mech] = fzero(@(r) failure(r,B,max_or_min), x, options);
       [stressdiff, mech, out] = failure(r,B,max_or_min);
       rvec(i) = r;
       outs{i} = out;
       switch mech
           case "r/z"
               rz(i) = r;
           case "z/r"
               zr(i) = r;
           case "z/t"
               zt(i) = r;
           case "t/z"
               tz(i) = r;
           case "r/t"
               rt(i) = r;
           case "t/r"
               tr(i) = r;
       end
       disp('Found one!')
       disp(r)
   catch ME
      disp('FAILED')
      disp(ME)
      fail(i) = 1;
   end
end

%csvwrite('no_shear.csv',[lambdas; rvec])


function [stressdiff, mech, out] = failure(r,A,max_or_min)

shear = true;

A.r = r;

out = muphem('multiflow2',0,A);

if shear 
    if max_or_min == "min"
        stressdiff = out{11};
        mech = out{19}{1};
    else % must be max
        stressdiff = out{13};
        mech = out{19}{2};
    end
else
    if max_or_min == "min"
        stressdiff = out{12}; % no shear is 12, shear is 11
        mech = out{19}{3};
    else % must be max
        stressdiff = out{14};
        mech = out{19}{4};
    end
end

end
