function [] = plotBalanceEquations(out)
%% Compare momentum balance
A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};

% Mass balance
%plot(sum(A.LHS(:,3,:),3),zvec, sum(A.RHS.eq3(:,3)),zvec)
LHS = squeeze(A.LHS(:,1,:));
figure; subplot(131)
plot(LHS(:,1),zvec,'DisplayName','dpdz'); hold on;
plot(LHS(:,2),zvec,'DisplayName','dphidz');
plot(LHS(:,3),zvec,'DisplayName','ddeltaudz');
legend
title('Mass balance')

% First momentum balance eq
LHS = squeeze(A.LHS(:,2,:));
RHS = -A.RHS.eq2;

subplot(132)
plot(LHS(:,1),zvec,'DisplayName','dpdz'); hold on;
plot(LHS(:,2),zvec,'DisplayName','dphidz');
plot(LHS(:,3),zvec,'DisplayName','ddeltaudz');
plot(RHS(:,1),zvec,'DisplayName','gravity'); hold on;
plot(RHS(:,2),zvec,'DisplayName','wall friction (melt)');
plot(RHS(:,3),zvec,'DisplayName','wall friction (gas)');
legend
title('Momentum balance')

% Second momentum balance eq (with interphase friction)
LHS = squeeze(A.LHS(:,3,:));
RHS = -A.RHS.eq3;

subplot(133)
plot(LHS(:,1),zvec,'DisplayName','dpdz'); hold on;
plot(LHS(:,2),zvec,'DisplayName','dphidz');
plot(LHS(:,3),zvec,'DisplayName','ddeltaudz');
plot(RHS(:,1),zvec,'DisplayName','gravity'); hold on;
plot(RHS(:,2),zvec,'DisplayName','wall friction (melt)');
plot(RHS(:,3),zvec,'DisplayName','wall friction (gas)');
plot(RHS(:,4),zvec,'DisplayName','interphase force');
legend
title('Momentum balance')



end

