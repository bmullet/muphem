function [] = plotfailureprofiles(A,Srr,Szz,Stt,Srz,zvec,pvec)


phi = 38/180*pi;
cohesion = 5e6; mu = tan(phi); C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

fail_rz = @(p,tau,c,q)   .5*(p+1+((p-1)^2 + 4*(tau*p)^2)^0.5) - c - q*0.5*(p+1-((p-1)^2 + 4*(tau*p)^2)^0.5);
fail_rt = @(p,s,tau,c,q) -2*s + p + c + q*0.5*(p + 1 - ((p-1)^2 + 4*(tau*p)^2)^.5);
fail_tz = @(p,s,tau,c,q) 0.5*(p + 1 + ((p-1)^2 + 4*(tau*p)^2)^.5) - c - q*(2*s-p);

frz_low  = nan(length(zvec));
frz_high = nan(length(zvec));
frt = nan(length(zvec));
ftz = nan(length(zvec));

s = 1; % Horizontal stress gradient

for i = 1:length(zvec)
    
    
    
end

figure
subplot(121)
plot(Srz,zvec); hold on;
plot(Srr,zvec);
legend('\tau','p')
ylabel('z')
xlabel('Pa')

subplot(122)
plot(Srz./Srr,zvec);
ylabel('z')
xlabel('\tau/p')

figure
plot(Srr./Szz, zvec);
xlabel('p/S_{zz}')
ylabel('z');

end