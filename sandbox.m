%% Test 1 - Test for uniqueness
A = initA();
z = 1:A.depth;
plith = z*A.k.rho*A.g;
A.Pchamber = max(plith)-60e6;
p = [];
%vvec = [.01:.01:3 3:.1:11 11:.002:12];
vvec = [.01:.01:.2 .2:.1:10];
%vvec = [.001:.002:.05 .2:.1:10];

z = [];
for i = vvec
    A.v_chamber_i = i;
    [zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A);
    p = [p min(pvec)];
    z = [z max(zvec)];
    disp(i)
end
plot(vvec,p)
hold on
xlabel('v')
ylabel('p top')
plot(xlim,[1.01e5 1.01e5],'--r')

figure
plot(vvec,z)
xlabel('v')
ylabel('top z')
%% Test 2 - Vary k
k = logspace(-10,-13,10);
vi = [];
pt = [];
for i = k
    disp(['k = '])
    disp(i)
    out = muphem('multiflow2',10e6,i);
    p = out{3};
    vg = out{4};
    vi = [vi min(vg)];
    pt = [pt min(p)];
end
figure
plot(k,pt)
xlabel('k')
ylabel('pressure')
set(gca,'XScale','log')
figure
plot(k,vi)
xlabel('k')
ylabel('chamber v')
set(gca,'XScale','log')

%% Test permeability law
kw0 = 1e-13;
z = 0:10000;
A = initA();
ppore = A.rhow*A.g*z+A.Patm_; % hyrdostatic pore pressure
kw1 = A.kw0*exp(-ppore/A.Pstar);

kw2 = 10.^(-14 - 3.2*log10(z/1000+.2));
plot(kw1,z,kw2,z)
legend('k1','k2')
set(gca,'XScale','log','Ydir','reverse')

%% Distribute solution
A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
disp('Output distributed')

%% Test 3 - Test for uniqueness
A = initA();
z = 1:A.depth;
plith = z*A.k.rho*A.g;
vvec = [.05:.2:15];
op = [-120e6 -100e6 -80e6 -40e6  -20e6 0];

p = zeros(length(op), length(vvec));
z = zeros(length(op), length(vvec));
for j = 1:length(op)
    A.Pchamber = max(plith)+op(j);
    for i = 1:length(vvec)
        A.v_chamber_i = vvec(i);
        try
            [zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A);
            p(j,i) = min(pvec);
            z(j,i) = max(zvec);
        catch
            p(j,i) = 10e6;
            z(j,i) = 500;
        end
        disp(vvec(i))
    end
    figure
    plot(vvec,p(j,:))
    hold on
    xlabel('v')
    ylabel('p top')
    plot(xlim,[1.01e5 1.01e5],'--r')
    title(['Overpressure = ' num2str(op(j)/1e6) ' MPa'])
    figure
    plot(vvec,z(j,:))
    xlabel('v')
    ylabel('top z')
    title(['Overpressure = ' num2str(op(j)/1e6) ' MPa'])
    drawnow()

end

%% Plot results of test three

figure
for i = 1:length(op)
    plot(vvec,p(i,:))
    hold on
end
xlabel('v')
ylabel('p top')
plot(xlim,[1.01e5 1.01e5],'--r')
legend('-120Mpa', '-100', '-80', '-40', '-20e6', '0')
title('Pressure at surface vs. initial velocity')


%% Test for interphase transition
ed = A.exdepth;
fd = A.fragdepth;
A.delF = 1;
bgi = (zvec>=ed&zvec<fd);
A.delF = 1;
delF = 1;
RHS1 = zeros(nnz(bgi),3);
phibg = phivec(bgi);
umbg = umvec(bgi);
ugbg = ugvec(bgi);
rhogbg = rhogvec(bgi);
pbg = pvec(bgi);
Qmbg = Qmvec(bgi); chidbg = chidvec(bgi);
zbg = zvec(bgi);
Fmgex = zeros(nnz(bgi),1);
F1 = Fmgex;
rbvec = Fmgex;
phitestex = Fmgex;
Is = Fmgex;
Fw = Fmgex;
Fmg1ex = Fmgex; Fmg2ex = Fmgex;
for i = 1:nnz(bgi)
    [Fmgex(i), phitestex(i), Fmg1ex(i), Fmg2ex(i)] = testfmg(zbg(i),[pbg(i), ugbg(i), umbg(i)],A);
end
fi = (zvec>=fd);
A.delF = 0;
delF = 0;
phif = phivec(fi);
umf = umvec(fi);
ugf = ugvec(fi);
rhogf = rhogvec(fi);
pf = pvec(fi);
Qmf = Qmvec(fi); chidf = chidvec(fi);
Fmgfr = zeros(nnz(fi),1);
zf = zvec(fi);
Fmgfr = zeros(nnz(fi),1);
F1 = Fmgfr;
Is = Fmgfr;
Fw = Fmgfr;
phitestfr = Fmgfr;
Fmg1fr = Fmgfr;
Fmg2fr = Fmgfr;
j = nnz(bgi);
for i=1:nnz(fi)
    [Fmgfr(i), phitestfr(i), Fmg1fr(i), Fmg2fr(i)] = testfmg(zf(i),[pf(i), ugf(i), umf(i)],A);
end
Fmg = [Fmgex; Fmgfr];
phitest = [phitestex; phitestfr];
Fmg1 = [Fmg1ex; Fmg1fr];
Fmg2 = [Fmg2ex; Fmg2fr];
testi = bgi | fi;

%% Bring down overpressure 
for p = 0e6:-1e6:-50e6
    disp(p)
    out = muphem('multiflow2',p);
end

%% Test for poissons ratio
A = initA();
for p = .1:.1:.5
    A.k.p = p;
    A.k.K = A.k.p/(1-A.k.p);
    out = muphem('multiflow2',0e6,A);
end


%% Plot results
%plotmuphem(A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec)

% Test for failure
[Srr, Szz, Stt] = kirsch(zvec,pvec,A);
[Smax,Sfail,failure] = mcfailure(A,Srr,Szz,Stt,zvec);
if (failure)
    disp('we have a failure!')
else
    disp('no failure')
end
plot_failure(zvec,pvec,phivec,A,Srr,Szz,Stt,Smax,Sfail,op)
disp('Pressure at conduit exit:')
disp(min(pvec));
vargout = {A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,failure};

%% Plot stability tests
load('findstablerhgfunc.mat')
plot(hgvec,stabler(2,:),'Color',[148/255, 20/255, 155/255]);
hold on;
plot(hgvec,stabler(3,:),'Color',[32/255, 89/255, 204/255]);
plot(hgvec,stabler(4,:),'Color',[32/255, 204/255, 49/255]);
title('Stable Radius v. Water Content')
xlabel('Initial Water Mass Fraction')
ylabel('Minimum Stable Radius (m)')
legend('10^5','10^6','10^7')
set(gca,'yscale','log')
xlim([0.03,0.07])

figure
load('GOLD2.mat')
plot(muvec,stabler(:,1));
hold on;
plot(muvec,stabler(:,2));
plot(muvec,stabler(:,3));
title('Stable Radius v. Water Content')
xlabel('Initial Water Mass Fraction')
ylabel('Minimum Stable Radius (m)')
legend('0.03','0.05','0.07')
set(gca,'yscale','log')
set(gca,'xscale','log')
