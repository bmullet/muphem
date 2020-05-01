% demo for smoothing

eta = 0.2;

f1 = @(x) -1./(x-eta-eta/10);
f2 = @(x) (x+eta).^2;

% f1 = @(x) x^2 + 4;
% f2 = @(x) x;

t = @(x) (x - (-eta))/(2*eta);

xv = -1:.0001:1;
f = nan(size(xv));

for i = 1:length(xv)
  x = xv(i);
  if x < -eta
      f(i) = f1(x);
  elseif (x >= -eta) && (x < eta)
      f(i) = f1(x)^(1-t(x))*f2(x)^t(x);
  else
      f(i) = f2(x);
  end
end

plot(xv,f); hold on

%do tanh smoothing

cubic = @(x) -2*(t(x)^3 - 3/2*t(x)^2);

for i = 1:length(xv)
  x = xv(i);
  if x < -eta
      f(i) = f1(x);
  elseif (x >= -eta) && (x < eta)
      f(i) = f1(x)*(1-cubic(x)) + f2(x)*(cubic(x)); 
  else
      f(i) = f2(x);
  end
end

plot(xv,f); hold on