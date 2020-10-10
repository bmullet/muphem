
A = Amodels.initA_paper();
A.phi0 = 0.70;
A.phiforce = 0.701;
%A.chamber_fac = 1;
A.r = 30;
A = Amodels.initA_paper(A);
A.fricfac = 1;
A.chamber_fac = 0.77;

A.Pchamber = (1.01e5+A.depth*A.g*A.k.rho)*A.chamber_fac;

u0s = [0.084, 0.25, 0.52];

A.u0 = u0s(3);

out = muphem('multiflow2',0,A);

%% Compare momentum balance
A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};

%%
MassLHS = (A.LHS(:,3,:));
plot(sum(A.LHS(:,3,:),3),zvec, A.RHS(:,3),zvec)

%%
figure
subplot(1,3,1)
plot(sum(A.LHS(:,3,:),3),zvec,A.RHS(:,3),zvec)
subplot(1,3,2)
plot(sum(A.LHS(:,2,:),3),zvec,A.RHS(:,2),zvec)
subplot(1,3,3)
plot(sum(A.LHS(:,1,:),3),zvec,A.RHS(:,1),zvec)