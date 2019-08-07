% RADIUSTEST tests whether a growing radius changes direction of failure

load('+Amodels/conduit_collapse_A')

rvec = [50; 75; 100; 125; 150; 175; 200; 225; 250; 275; 300];
out = cell(size(rvec));

% Collect solution
for i = 1:length(rvec)
    A.r = rvec(i);
    out{i} = muphem('multiflow2',0,A);
    
end

% Plotting
% unpack solution
%% Distribute solution
n1 = nan(3,length(rvec));
n2 = nan(3,length(rvec));
fail = zeros(1,length(rvec));
fragdepth = zeros(1,length(rvec));
Sfailure = nan(length(rvec),3);


for i = 1:length(rvec)
    A = out{i}{1}; zvec = out{i}{2}; pvec = out{i}{3}; ugvec = out{i}{4}; umvec = out{i}{5};
    phivec = out{i}{6}; rhogvec = out{i}{7}; chidvec = out{i}{8}; Qmvec = out{i}{9}; Qgvec = out{i}{10}; failure = out{i}{11}; Sprincipal = out{i}{12};
    slip = out{i}{13}; Sfail = out{i}{14};
    
    n1(:,i) = slip{1};
    n2(:,i) = slip{2};
    fail(i) = failure;
    fragdepth(i) = A.fragdepth;
    Sfailure(i,:) = Sfail;
end

disp('Output distributed')

%% Plot
figure()
subplot(141)
plot(rvec,fragdepth,'ok')
xlabel('r');
ylabel('frag loc (positive = more shallow)')

subplot(142)

scatter(rvec,n1(1,:).^2,[],fail); hold on
scatter(rvec,n2(1,:).^2,[],fail);
title('%r')
xlabel('r')

subplot(143)
scatter(rvec,n1(2,:).^2,[],fail); hold on
scatter(rvec,n2(2,:).^2,[],fail);
title('%\theta')
xlabel('r')

subplot(144)
scatter(rvec,n1(3,:).^2,[],fail); hold on
scatter(rvec,n2(3,:).^2,[],fail);
xlabel('r')
title('%z')

figure()
plot(rvec,Sfailure(:,1),'or',rvec,Sfailure(:,2),'ok',rvec,Sfailure(:,3),'oc')
xlabel('r (m)')
ylabel('Pa')
legend('\sigma_1','\sigma_2','\sigma_3')