    function Fmg = interphase(duv, pv, phiv, rhogv,A)
        Fmg = nan(size(duv));
        
        
            for n = 1:length(duv)
                du = duv(n);
                phi = phiv(n);
                rhog = rhogv(n);
                p = pv(n);
                
                rb = (3*phi/(4*pi*A.nb))^(1/3); %bubble radius
              
                
                rb = (3*phi/(4*pi*A.nb))^(1/3); %bubble radius
                rbnd = rb/A.C.rb0;
                
                if (A.useForchheimer)
                    k1 = (phi/A.phi0)^(A.m+2/3);
                    k2 = (phi/A.phi0)^((1+3*A.m)/2+1/3);
                    Fmg1 = -1/A.C.St*(1+A.C.Fo*k1/k2*abs(du)*rhog)*phi*(1-phi)/k1*du;
                    
                else
                    Fmg1 = -1/A.C.Reb * A.C.rc/A.C.rb0 * 1/rbnd^2 * (9/2*phi*(1-phi)*A.mu(phi,p*A.C.p0)/A.C.mu0*(du));
                end
                
                
                Fmg2 = -3/8*phi*(1-phi)*A.dragC*A.C.rc/A.Rash*rhog*abs(du)*(du)*A.C.delta*100;
                
                
                if A.delF
                    Fmg(n) = Fmg1;
                else
                    Fmg = Fmg2;
                end
                
            end