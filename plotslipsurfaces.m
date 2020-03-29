function [slip] = plotslipsurfaces(zvec,S,A,plot)

[~,i] = min(abs(zvec - (A.fragdepth*1.001)));

eV = squeeze(S(i,:,:,:));

% Rotate
pp = [1 0 0]; % point for plotting
omega = pi/4+deg2rad(30/2);

Rpos = @(omega) [[cos(omega) 0 sin(omega)];
     [0 1 0];
     [-sin(omega) 0 cos(omega)]];
Rneg = @(omega) [[cos(-omega) 0 sin(-omega)];
     [0 1 0];
     [-sin(-omega) 0 cos(-omega)]];

n1 = [1 0 0]'; % Same direction as sigma1
n1 = Rpos(omega)*n1; % Rotate
n1 = eV*n1; % Change basis to r, z, theta space

n2 = [1 0 0]'; % Same direction as sigma1
n2 = Rneg(omega)*n2; % Rotate
n2 = eV*n2; % Change basis to r, z, theta space

w = null(n1'); % Find two orthonormal vectors which are orthogonal to v
[P,Q] = meshgrid(-1:.1:1); % Provide a gridwork (you choose the size)
X1 = pp(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y1 = pp(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z1 = pp(3)+w(3,1)*P+w(3,2)*Q;

w = null(n2'); % Find two orthonormal vectors which are orthogonal to v
[P,Q] = meshgrid(-1:.1:1); % Provide a gridwork (you choose the size)
X2 = pp(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y2 = pp(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z2 = pp(3)+w(3,1)*P+w(3,2)*Q;

% Re-do rotations to find slip directions
omega = 3*pi/4+deg2rad(30/2);
n1 = [1 0 0]'; % Same direction as sigma1
n1 = Rpos(omega)*n1; % Rotate
n1 = eV*n1; % Change basis to r, z, theta space

n2 = [1 0 0]'; % Same direction as sigma1
n2 = Rneg(omega)*n2; % Rotate
n2 = eV*n2; % Change basis to r, z, theta space

om = linspace(0,2*pi,30);
x = cos(om);
y = sin(om);

if(plot)
fprintf('-----------SLIP DIRECTIONS---------\n')
fprintf('      %%r        %%z         %%theta\n')
fprintf('-----------------------------------\n')
fprintf('      %0.2f      %.02f       %.2f\n',n1(1)^2, n1(3)^2, n1(2)^2)
fprintf('      %0.2f      %.02f       %.2f\n',n2(1)^2, n2(3)^2, n2(2)^2)


    figure()
    
    ax1 = subplot(1,2,1);
    % Plot circle, slip directions, fault surfaces
    plot3(x,y,zeros(size(x))); hold on;
    quiver3([1 1],[0 0],[0 0],[n1(1) n2(1)],[n1(2) n2(2)],[n1(3) n2(3)])
    surf(X1,Y1,Z1)
    surf(X2,Y2,Z2)
    
    
    ax2 = subplot(1,2,2);
    % Plot a circle, slip directions, fault surfaces
    plot3(x,y,zeros(size(x))); hold on;
    quiver3([1 1],[0 0],[0 0],[n1(1) n2(1)],[n1(2) n2(2)],[n1(3) n2(3)])
    
    % Set views
    view(ax2,[0 0]);
    view(ax1,[0 90]);
    subplot(1,2,2); zlim([-1 1]);
    subplot(1,2,1); zlim([-1 1]);
end

% Return slip vectors
slip = {n1,n2};

end