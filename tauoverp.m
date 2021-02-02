% Figure for max tau/p

phis = [0.60, 0.70, 0.80];
sol = cell(size(phis));

for j = 1:length(phis)

A = Amodels.initA_paper_dacite;
A.phi0 = phis(j);

rvec = 10:100;
taus = cell(size(rvec));
ps = cell(size(rvec));

parfor i = 1:length(rvec)
    B = A;
    B.r = rvec(i);
    disp(["Working on " num2str(i)])
    try
        out = muphem('multiflow2',0,B);
    catch ME
        disp(["Failed at " num2str(rvec(i))])
    end
    taus{i} = out{18};
    ps{i} = out{17};
end

soln.taus = taus;
soln.ps = ps;
sol{j} = soln;


end

%%
tp = nan(size(rvec));
for j = 1:length(phis)
    taus = sol{j}.taus;
    ps = sol{j}.ps;
for i = 1:length(rvec) 
    tau = taus{i};
    p = ps{i};
    tp(i) = max(tau./p);
end
semilogy(rvec, tp); hold on;
end
legend("0.6", "0.7", "0.8")

%%
% Figure for max tau/p
% vary phi

phis = 0.6:.01:0.8;

rvec = [20,50,100];
sol = cell(size(rvec));

for j = 1:length(rvec)

A = Amodels.initA_paper_dacite;
A.r = rvec(j);


taus = cell(size(phis));
ps = cell(size(phis));

parfor i = 1:length(phis)
    B = A;
    B.phi0 = phis(i);
    disp(["Working on " num2str(i)])
    try
        out = muphem('multiflow2',0,B);
    catch ME
        disp(["Failed at " num2str(phis(i))])
    end
    taus{i} = out{18};
    ps{i} = out{17};
    zs{i} = out{2};
end

soln.taus = taus;
soln.ps = ps;
soln.zs = zs;
sol{j} = soln;


end

%%
figure

tstar = @(c,q, pp)  ((1./((c + q + pp.*(1-q)).*(1-c+pp.*(q-1)))).^(1/2).*(c + pp + q.*(1-pp) - 1))/2;
phi = 38/180*pi;
q = tan(pi/4 + 1/2*phi).^2;
pp = 1/2.7;

tp = nan(size(phis));
ts = nan(size(phis));
colors = jet(3)*.9;
handles = [];

for j = 1:length(rvec)
    taus = sol{j}.taus;
    ps = sol{j}.ps;
    zs = sol{j}.zs;
for i = 1:length(phis) 
    tau = taus{i};
    p = ps{i};
    [tp(i), k] = max(tau./p);
    z = zs{i};
    zfrag = z(k);
    cohesion = A.mc.C(zfrag);
    Szz = -zfrag*2700*9.8;
    mu = tan(phi);
    c = 2*cohesion*((mu^2 + 1)^(1/2) + mu)./Szz;
    ts(i) = tstar(c,q,pp);
end
h = plot(phis, tp,"Color",colors(j,:)); hold on;
handles = [handles h];
plot(phis, ts,"--","Color",colors(j,:)); hold on;
end
legend(handles, "R=20", "R = 50", "R = 100","Location","northwest")
%ht = text(0.7, 0.5, {'{\color{[0 0 0]} o } Red', '{\color{blue} o } Blue', '{\color{black} o } Black'}, 'EdgeColor', 'k');
ylabel("\tau_{frag}/p_{frag}, \tau*")
xlabel("\phi_{frag}")
xlim([0.6,0.80])
grid on


