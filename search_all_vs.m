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
%A.chamber_fac = 1;
A.r = 30;
A = Amodels.initA_paper(A);
A.fricfac = 1;
phis = [.7, .77, .85, 1];

% try with Degruyter

% A = Amodels.initA_Degruyter2012;
% A.phi0 = 0.70;
% A.phiforce = 0.7001;
% A.chamber_fac = 0.782;
% A.r = 30;
%
% phis = [1, 0.8, 0.7, 0.6];
%phis = [1];



%%

%phis = [100];
vs = [0.01:.01:2];

resids = nan(length(phis),length(vs));
phi0s = nan(length(vs), 3);
fricfacs = nan(size(resids));


for i = 1:length(phis)
    A.chamber_fac = phis(i);
    %     disp(A.phiforce)
    textprogressbar(sprintf('phif: %.2d\n',A.phiforce));
    A.Pchamber = (1.01e5+A.depth*A.g*A.k.rho)*A.chamber_fac;
    
    for j = 1:length(vs)
        textprogressbar(j/length(vs)*100);
        v = vs(j);
        A.v_chamber_i = v;
        
        A.fricfac = 2;
        counter = 1;
        
        zvec = nan;
        
        
        while isnan(zvec)
            A.fricfac = A.fricfac * 0.5;
            [zvec,pvec,~,~,~,~,~,~,~,Ar] = incoodes(A);
            
            if counter > 20
                break
            end
            
        end
        if isnan(zvec)
            resid = nan;
        elseif max(zvec) < 0
            % did not make it to surface
            resid = abs(max(zvec))*10;
        else
            % made it to surface but pressure is too high
            resid = (A.Patm_-pvec(end))/1e5*40;
        end
        
        resids(i,j) = resid;
        fricfacs(i,j) = A.fricfac;
        
        %phi0s(j,:) = Ar.du0;
        
    end
    
    textprogressbar('done!')
    
    plot(vs,resids(i,:),'DisplayName',num2str(phis(i))); hold on
    plot(xlim,[0,0],'--r','HandleVisibility','off')
    %ylim([-100, 100])
    legend('show')
    drawnow()
    
end

