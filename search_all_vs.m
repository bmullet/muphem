% test all vs
% 
% A = Amodels.initA_MSH;
% A.phi0 = .70;
% A.phiforce = 0.75;
% A.chamber_fac = 0.80;
% A.lambda = 0.50;
% 
% A = Amodels.initA_MSH(A);
% 
% % Search over vs
% 
% %vs = linspace(0.1,10,400);
% %vs = logspace(-1,log10(15),1000);

% should give three solutions

A = Amodels.initA_paper();
A.phi0 = 0.70;
A.phiforce = 0.701;
A.chamber_fac = 0.782;
A.r = 30;
A = Amodels.initA_paper(A);
phis = [0.782];

% try with Degruyter

% A = Amodels.initA_Degruyter2012;
% A.phi0 = 0.70;
% A.phiforce = 0.7001;
% A.chamber_fac = 0.782;
% A.r = 30;
% 
% phis = [1, 0.8, 0.7, 0.6];

vs = [0.01:.02:1]; 

resids = nan(length(phis),length(phis));
%%
for i = 1:length(phis)
    A.chamber_fac = phis(i);
    disp(A.phiforce)
    textprogressbar(sprintf('phif: %.2d\n',A.phiforce));
    A.Pchamber = (1.01e5+A.depth*A.g*A.k.rho)*A.chamber_fac;
    
for j = 1:length(vs) 
   textprogressbar(j/length(vs)*100);
   v = vs(j);
   A.v_chamber_i = v;
   
   [zvec,pvec,~,~,~,~,~,~,~,Ar] = incoodes(A);
   
   if max(zvec) < 0
       % did not make it to surface
       resid = abs(max(zvec))*10;
   else
       % made it to surface but pressure is too high
       resid = (A.Patm_-pvec(end))/1e5*40;
   end
   
   resids(i,j) = resid;
  
end
textprogressbar('done!')

plot(vs,resids(i,:),'DisplayName',num2str(phis(i))); hold on
plot(xlim,[0,0],'--r')
%ylim([-100, 100])
legend('show')
drawnow()

end

