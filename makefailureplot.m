%% Make failure plot
%% Create pcolor plots from representative data
set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 20, ...
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1, ...
      'DefaultLineLineWidth',3) ;
  
close all;
load('failureandviscosity3.mat');
pvvecplt = pvvec/1e6; %MPa
figure
Qtot = vbot*3000;

itotal = length(muvvec)-2;
legplt = cell(itotal,1);
for i = 1:itotal
    xlim([-50 30]);
    ind = (failurevec(i,:)<1);
    Qtot(i,ind)
    plot(pvvecplt(ind),Qtot(i,ind));
    hold on
   
end
legend('1x10^3','3x10^3','8x10^3','2x10^4','6x10^4','2x10^5','5x10^4','1x10^6')
