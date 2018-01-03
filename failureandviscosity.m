%% Make plots of various viscosities and failure points

A = initA();
muvvec = logspace(3,7,10);
r = 75;
pvvec = [-60e6:2e6:20e6];
A.r = r;

Nmu = length(muvvec);
Np = length(pvvec);

it = 0;
N = Nmu*Np;   

vbot = nan(Nmu,Np);
vgtop = nan(Nmu,Np);
vmtop = nan(Nmu,Np);
ptop = nan(Nmu,Np);
failurevec = nan(Nmu,Np);

for i = 1:length(muvvec)
        for k = 1:length(pvvec)
            
            A.mu = muvvec(i);
            p = pvvec(k);
            
            disp([num2str(it) ' out of ' num2str(N)]);
            disp(['Progress = ' num2str(it/N*100)]);
            disp(['mu = ' num2str(A.mu)]);
            disp(['p_ch = ' num2str(p)]);
            disp(['R = ' num2str(A.r)]);
            try
                out = muphem('multiflow2',p,A);
                A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
                phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
                disp('Output distributed')
                vbot(i,k) = min(umvec);
                vgtop(i,k) = max(ugvec);
                vmtop(i,k) = max(umvec);
                ptop(i,k) = min(pvec);
                failurevec(i,k) = failure;
            catch
                disp('Code failed!')
            end
        end
end
save('failureandviscosity3.mat')

