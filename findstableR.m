function [stabler] = findstableR()
A = initA();
N = 10;
muvec = logspace(4,7,N);
hgvec = linspace(.03,.07,N);

[rmaxvec, rminvec] = setrmaxmin(muvec);

tol = 1; % resolve to the nearest meter

stabler = zeros(N,N);

for i = 1:N
    % Loop through viscosities
    mu = muvec(i);
    rmax = rmaxvec(i);
    rmin = rminvec(i);
    
    for j = 1:N
        % Loop through water contents
        hg = hgvec(j);
        A.hg = hg;
        A.mu = mu;
        stabler(i,j) = findr(A,rmax,rmin,tol);
        
        
    end
end




  
end


function [rmaxvec, rminvec] = setrmaxmin(muvec)
    muset = [4 5 6 7];
    lb = [3 20 70 220];
    ub = [15 35 110 350];
    logmu = log10(muvec);
    rmaxvec = interp1(muset,ub,logmu);
    rminvec = interp1(muset,lb,logmu);

end


  function [r] = findr(A,rmax,rmin,tol)
        % First find bounds on stability
        A.r = rmax;
        
        try
            out = muphem('multiflow2',0e6,A);
            failure = out{11};
        catch
            failure = 1;
        end
        

        if (failure)
            r = 999;
            disp('rmax fails!')
            return
        end

        
        A.r = rmin;
        try
            out = muphem('multiflow2',0e6,A);
            failure = out{11};
        catch
            failure = 1;
        end
        
        if (~failure)
            disp('stable under 5m')
            r = 5;
            return
        end
        
        d = rmax-rmin;
        rtop = rmax;
        rbot = rmin;
        
        while d>tol
           rtest = rbot + d/2;
           A.r = rtest;
           try
               out = muphem('multiflow2',0e6,A);
               failure = out{11};
           catch
               failure = 1;
           end
           
           if (failure)
               %means we have a new rbot
               rbot = rtest;
           else
               rtop = rtest;
           end
           d = rtop-rbot;
           disp(d)
            
        end
        
        r = rtest;
        
    end
