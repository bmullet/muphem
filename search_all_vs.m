% test all vs

A = Amodels.initA_MSH;
A.phi0 = .70;
A.phiforce = 0.75;
A.chamber_fac = 0.80;
A.lambda = 0.50;

A = Amodels.initA_MSH(A);

% Search over vs

%vs = linspace(0.1,10,400);
%vs = logspace(-1,log10(15),1000);
vs = [0.01:.01:1,1.1:.05:15]; 

rs = [30, 70, 100, 150, 200, 250, 400];

resids = nan(length(rs),length(vs));
%%
for i = 1:2
    A.r = rs(i);
    disp(A.r)
    textprogressbar(sprintf('radius: %.2d\n',A.r));
for j = 1:length(vs) 
%    textprogressbar(j/length(vs)*100);
%    v = vs(j);
%    A.v_chamber_i = v;
%    
%    [zvec,pvec,~,~,~,~,~,~,~,A] = incoodes(A);
%    
%    if max(zvec) < 0
%        % did not make it to surface
%        resid = abs(max(zvec))*10;
%    else
%        % made it to surface but pressure is too high
%        resid = (A.Patm_-pvec(end))/1e5*40;
%    end
%    
%    resids(i,j) = resid;
  
end
textprogressbar('done!')

plot(vs,resids(i,:),'DisplayName',num2str(rs(i))); hold on
plot(xlim,[0,0],'--r')
ylim([-100, 100])
legend('show')
drawnow()

end

