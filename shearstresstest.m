% Compare wall shear stress and normal stress

% Get shear stresses

shearstress = nan(size(umvec));
ss_bbl = 8*mu*umvec./A.r;
ss_frag = A.f0.*rhogvec.*ugvec.^2/A.r;

shearstress(zvec<A.fragdepth) = ss_bbl(zvec<A.fragdepth);
shearstress(zvec>=A.fragdepth) = ss_frag(zvec>=A.fragdepth);

plot(shearstress, zvec, pvec, zvec);
legend('Shear','Normal')