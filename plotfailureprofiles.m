function [] = plotfailureprofiles(A,Srr,Szz,Stt,Srz,zvec,pvec)

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