% test all vs

A = Amodels.initA_MSH;
A.phi0 = .70;
A.phiforce = 0.75;
A.chamber_fac = 0.80;
A.lambda = 0.50;

A = Amodels.initA_MSH(A);

A.r = 250;

% Search over vs

%vs = linspace(0.1,10,400);
vs = logspace(-1,1,400);
resids = nan(size(vs));

%%

for i = 1:length(vs)
   v = vs(i);
   A.v_chamber_i = v;
   
   [zvec,pvec,~,~,~,~,~,~,~,A] = incoodes(A);
   
   if max(zvec) < 0
       % did not make it to surface
       resid = abs(max(zvec))*10;
   else
       % made it to surface but pressure is too high
       resid = (A.Patm_-pvec(end))/1e5*40;
   end
   
   resids(i) = resid;
  
end

plot(vs,resids); hold on
plot(xlim,[0,0],'--r')