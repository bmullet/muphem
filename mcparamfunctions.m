%% Helper functions for converting between MC parameters

% Coefficient of internal friction from angle of internal friction
mu = @(phi) tan(phi);

% Angle of internal friction from coefficient of internal friction
phi = @(mu) atan(mu);

% Angle between sigma_1 and the fault (assuming MC)
theta = @(phi) pi/4 + 1/2*phi;

% Compressive strength from cohesion and angle coefficient of internal
% friction
C0 = @(c, mu) 2*c*((mu+1)^.5 + mu);
