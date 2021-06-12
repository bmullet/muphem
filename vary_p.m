A = Amodels.initA_paper_MSH_A;
pmax = A.Pchamber + 10e6;
pmin = A.Pchamber - 10e6;
pvec = linspace(pmin, pmax, 10);


rvec = [20, 40, 60];

Q = zeros(length(pvec), length(rvec));

for j = 1:length(rvec)
parfor i = 1:length(pvec)
  B = A;
  B.r = rvec(j);
  B.Pchamber = pvec(i);
  out = muphem("multiflow2",0,B);
  Q(i,j) = mean(out{9} + out{10});
end
end

%%
R = Q./max(max(Q)); 
for j = 1:length(rvec)
   
   plot(pvec, R(:,j)); hold on;
   
end
legend("R=20", "R=40", "R=60")
title("Flux for different p_{chamber}, R")
ylabel("Normalized flux")
xlabel("p_{chamber}")
