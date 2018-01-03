%% Bring down overpressure 

muvvec = [1e3, 1e4, 1e6, 1e8];
rvvec = [5, 25, 50, 100];
pvvec = [60e6:-30e6:-60e6];

Nmu = length(muvvec);
Nr = length(rvvec);
Np = length(pvvec);

it = 0;
N = Nmu*Nr*Np;   

vbot = nan(Nmu,Nr,Np);
vgtop = nan(Nmu,Nr,Np);
vmtop = nan(Nmu,Nr,Np);
ptop = nan(Nmu,Nr,Np);

for i = 1:length(muvvec)
    for j = 1:length(rvvec)
        for k = 1:length(pvvec)
            it = it+1;
            
            A.mu = muvvec(i);
            A.r = rvvec(j);
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
                vbot(i,j,k) = min(umvec);
                vgtop(i,j,k) = max(ugvec);
                vmtop(i,j,k) = max(umvec);
                ptop(i,j,k) = min(pvec);
            catch
                disp('Code failed!')
            end
        end
    end
end
save('explore.mat')

%%
% Plot up the results
%load('explore.mat')
for i = 1:length(muvvec)
    for j = 1:length(rvvec)
        plot(pvvec,reshape(vbot(i,j,:),1,5))
        disp('mu = ')
        disp(muvvec(i));
        disp('r = ')
        disp(rvvec(j));
        pause
    end
end