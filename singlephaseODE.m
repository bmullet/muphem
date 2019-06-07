function [ dydz ] = singlephaseODE( z,y,u0,A )

p = y(1,:);
N = size(y,1); dydz = zeros(1,N);
rho = A.rhom0;
tau = 4*A.mu(0,p)*u0/A.r; %laminar flow in a pipe

rhs = -(2*tau./A.r+rho*A.g);
dydz(1,:) = rhs;
