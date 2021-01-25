% Failure with shear tractions

% taup is tau/p, tauz is tau/z. all other variables are normalized by
% sigmaz, except for q.

syms k p C q taup tauz pp

S1 = 1/2*(p + 1) + sqrt((1/2*(p-1))^2 + (taup*p)^2);
S3 = 1/2*(p + 1) - sqrt((1/2*(p-1))^2 + (taup*p)^2);

Stheta = 2*k - p;

e1 = Stheta - pp == C + q*(S3 - pp);
sol1 = solve(e1, p);

rightplus = -(2*C - 4*k + 2*pp + q + C*q - 2*k*q - pp*q + q*(4*C^2*taup^2 + C^2 - 16*C*k*taup^2 - 4*C*k - 8*C*pp*q*taup^2 - 2*C*pp*q + 8*C*pp*taup^2 + 2*C*pp + 4*C*q*taup^2 + 2*C*q + 2*C + 16*k.^2*taup^2 + 4*k.^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 16*k*pp*taup^2 - 4*k*pp - 8*k*q*taup^2 - 4*k*q - 4*k + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 4*pp*q*taup^2 + 2*pp + q^2 + 2*q + 1).^(1/2) - pp*q^2 + q^2)./(2*(- q^2*taup^2 + q + 1));
rightminus = -(2*C - 4*k + 2*pp + q + C*q - 2*k*q - pp*q - q*(4*C^2*taup^2 + C^2 - 16*C*k*taup^2 - 4*C*k - 8*C*pp*q*taup^2 - 2*C*pp*q + 8*C*pp*taup^2 + 2*C*pp + 4*C*q*taup^2 + 2*C*q + 2*C + 16*k.^2*taup^2 + 4*k.^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 16*k*pp*taup^2 - 4*k*pp - 8*k*q*taup^2 - 4*k*q - 4*k + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 4*pp*q*taup^2 + 2*pp + q^2 + 2*q + 1).^(1/2) - pp*q^2 + q^2)./(2*(- q^2*taup^2 + q + 1));
tausingular = (q + 1)^(1/2)/q; % when the above solutions go singular
% the most far right value of k we will need if using both plus and minus
% solutions
kright = (C + pp + q - pp*q - 2*taup*(q^2*taup^2 - q - 1)^(1/2) + 4*C*taup^2 + 4*pp*taup^2 + 2*q*taup^2 - 4*pp*q*taup^2 + 1)/(2*(4*taup^2 + 1));


e2 = S1 - pp == C + q*(S3 - pp);
sol2 = solve(e2,p);
top = (C*q - pp - C + 2*pp*q - pp*q^2 + q*(4*C^2*taup^2 + C^2 - 8*C*pp*q*taup^2 - 2*C*pp*q + 8*C*pp*taup^2 + 2*C*pp + 4*C*q*taup^2 + 2*C*q - 4*C*taup^2 - 2*C + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1).^(1/2) + q^2 + (4*C^2*taup^2 + C^2 - 8*C*pp*q*taup^2 - 2*C*pp*q + 8*C*pp*taup^2 + 2*C*pp + 4*C*q*taup^2 + 2*C*q - 4*C*taup^2 - 2*C + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1).^(1/2) + 1)/(2*(q^2*taup^2 + 2*q*taup^2 + q + taup^2));
bottom = -(C + pp - C*q - 2*pp*q + pp*q^2 + q*(4*C^2*taup^2 + C^2 - 8*C*pp*q*taup^2 - 2*C*pp*q + 8*C*pp*taup^2 + 2*C*pp + 4*C*q*taup^2 + 2*C*q - 4*C*taup^2 - 2*C + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1)^(1/2) - q^2 + (4*C^2*taup^2 + C^2 - 8*C*pp*q*taup^2 - 2*C*pp*q + 8*C*pp*taup^2 + 2*C*pp + 4*C*q*taup^2 + 2*C*q - 4*C*taup^2 - 2*C + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1)^(1/2) - 1)/(2*(q^2*taup^2 + 2*q*taup^2 + q + taup^2));
top2 = (C + pp + q - pp*q - ((C + pp + q - pp*q - 1)*(C + pp + q - pp*q + 4*C*taup^2 + 4*pp*taup^2 - 4*pp*q*taup^2 - 1))^(1/2) - 1)/(2*taup^2*(q - 1));
bottom2 = (C + pp + q - pp*q + ((C + pp + q - pp*q - 1)*(C + pp + q - pp*q + 4*C*taup^2 + 4*pp*taup^2 - 4*pp*q*taup^2 - 1))^(1/2) - 1)/(2*taup^2*(q - 1));

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
e3 = S1 - pp == C + q*(Stheta - pp);
sol3 = solve(e3,p);
leftplus = (C + pp - q + 2*C*q + 2*k*q + pp*q + (4*C^2*taup^2 + C^2 + 16*C*k*q*taup^2 + 4*C*k*q - 8*C*pp*q*taup^2 - 2*C*pp*q + 8*C*pp*taup^2 + 2*C*pp - 2*C*q - 4*C*taup^2 - 2*C + 16*k.^2*q^2*taup^2 + 4*k.^2*q^2 - 16*k*pp*q^2*taup^2 - 4*k*pp*q^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 4*k*q^2 - 8*k*q*taup^2 - 4*k*q + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 + 2*pp*q^2 + 4*pp*q*taup^2 - 4*pp*taup^2 - 2*pp + q^2 + 2*q + 1).^(1/2) + 4*k*q^2 - 2*pp*q^2 - 1)./(2*(q^2 + q - taup^2));
leftminus = (C + pp - q + 2*C*q + 2*k*q + pp*q - (4*C^2*taup^2 + C^2 + 16*C*k*q*taup^2 + 4*C*k*q - 8*C*pp*q*taup^2 - 2*C*pp*q + 8*C*pp*taup^2 + 2*C*pp - 2*C*q - 4*C*taup^2 - 2*C + 16*k.^2*q^2*taup^2 + 4*k.^2*q^2 - 16*k*pp*q^2*taup^2 - 4*k*pp*q^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 4*k*q^2 - 8*k*q*taup^2 - 4*k*q + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 + 2*pp*q^2 + 4*pp*q*taup^2 - 4*pp*taup^2 - 2*pp + q^2 + 2*q + 1).^(1/2) + 4*k*q^2 - 2*pp*q^2 - 1)./(2*(q^2 + q - taup^2));

e4 = top == bottom;
taustar = solve(e4, taup);

%% paraboloa minimums
kright = solve(rightplus == rightminus, k);
% first solution from above
kright = (C + pp + q - pp*q - 2*taup*(q^2*taup^2 - q - 1)^(1/2) + 4*C*taup^2 + 4*pp*taup^2 + 2*q*taup^2 - 4*pp*q*taup^2 + 1)/(2*(4*taup^2 + 1));

solve((2*(- q^2*taup^2 + q + 1)) == 0, taup)


%% vary tau
phi = 35/180*pi;
C = 0; % this is C0/sigz
q = tan(pi/4 + 1/2*phi)^2;
%pp = 1/2.7;
pp = 0;

%taus = .365:.0001:.367;
%taus = 0:.1:.36;
%taus =  0.6:.05:1;
taus = [0.6, 0.65];
k = 0:.01:5;
kk = [k k]; % used for solutions with more than one k value
%k = 0.5:.001:1.5;
colors = parula(length(taus));
handles = [];
figure

disp("Max tau:"); disp(eval(subs(taustar)));

for i = 1:length(taus)
    taup = taus(i);
    tsingular = eval(subs(tausingular));
    c = colors(i,:);
    
    topval = eval(subs(top));
    bottomval = eval(subs(bottom));
    
    right = [eval(subs(rightminus)) eval(subs(rightplus))];
    left = eval(subs(leftminus));
    
    topline = ones(size(k))*topval;
    bottomline = ones(size(k))*bottomval;

    if taup > tsingular
        % need both solutions
        ki = k < eval(subs(kright));
        righti = [ki ki];
        righti = righti & (right > bottomval) & (right < topval);

    else
        % only need the minus solution
        righti = [ones(size(k)) zeros(size(k))];
        righti = righti & (right > bottomval) & (right < topval);
      
    end
    
    
    lefti = (left > bottomline) & (left < topline);
    
    toprightk = max(k(righti));
    topleftk = max(k(lefti));
    bottomrightk = min(k(righti));
    bottomleftk = min(k(lefti));
    
    topi = (k < toprightk) & (k > topleftk);
    bottomi = (k < bottomrightk) & (k > bottomleftk);

    %plot(k, subs(rightplus)); hold on
    h = plot(kk(righti), right(righti), 'Color', c, 'DisplayName', string(taup)); hold on;
    handles = [handles h];
    plot(k(topi), topline(topi), 'Color', c); hold on;
    plot(k(bottomi), bottomline(bottomi),'Color',  c);

    %plot(k, subs(leftplus)); hold on
    plot(k(lefti), left(lefti),'Color',  c); hold on;
    %xlim([0,4])
    %ylim([0,4])
    grid on
    xlabel('k = S/\sigma_{zz}')
    ylabel('p/\sigma_{zz}')
    
    
end

legend(handles)
title('MC failure envelope with varying \tau/p')

%% vary pp
phi = 35/180*pi;
C = 0; % this is C0/sigz
q = tan(pi/4 + 1/2*phi)^2;
taup = 0.0;

pps = 0;


colors = parula(length(pps));
%colors = [0 0 0];
handles = [];
figure

for i = 1:length(pps)
    pp = pps(i);
    c = colors(i,:);
    k = 0:.01:5;

    right = eval(subs(rightminus));
    left = eval(subs(leftminus));
    topline = ones(size(k))*eval(subs(top));
    bottomline = ones(size(k))*eval(subs(bottom));

    righti = (right > bottomline) & (right < topline);
    lefti = (left > bottomline) & (left < topline);
    topi = zeros(size(topline));
    bottomi = zeros(size(bottomline));
    topi(find(lefti,1,'last'): find(righti,1,'last')) = true;
    bottomi(find(lefti,1,'first'): find(righti,1,'first')) = true;

    topi = topi > 0.5;
    bottomi = bottomi > 0.5;


    %plot(k, subs(rightplus)); hold on
    h = plot(k(righti), right(righti), 'Color', c, 'DisplayName', string(pp)); hold on;
    handles = [handles h];
    plot(k(topi), topline(topi), 'Color', c); hold on;
    plot(k(bottomi), bottomline(bottomi),'Color',  c);

    %plot(k, subs(leftplus)); hold on
    plot(k(lefti), left(lefti),'Color',  c); hold on;
    xlim([0,4])
    ylim([0,4])
    grid on
    xlabel('k = S/\sigma_{zz}')
    ylabel('p/\sigma_{zz}')
    
    
end

legend(handles)
title('MC failure envelope with varying p_p/\sigma_{zz}')

%% Taustar with pp
tausf = @(c, pp, q) 1/2 * (c + pp + q.*(1-pp) - 1)./((c+q+pp.*(1-q)).*(1-c+pp.*(q-1))).^(1/2);

pp = 0:.001:.4;
Cs = [0:.2:1];
colors = parula(length(Cs));

for i = 1:length(Cs)
    c = colors(i,:);
    semilogy(pp, tausf(Cs(i), pp, q), 'Color', c); hold on;
end
legend(string(Cs));
grid on;
xlabel('p_p/\sigma_{zz}')
ylabel('\tau*')
title('\tau* v. pore pressure for varied C_0/\sigma_{zz}')


%% Get right boundary right

phi = 35/180*pi;
C = 0; % this is C0/sigz
q = tan(pi/4 + 1/2*phi)^2;
%pp = 1/2.7;
pp = 0;

%taus = .365:.0001:.367;
%taus = 0:.1:.36;
%taus =  0.6:.05:1;
taus = [0.57, 0.59];
k = 0:.02:5;
kk = [k k]; % used for solutions with more than one k value
%k = 0.5:.001:1.5;
%colors = parula(length(taus));
handles = [];
figure

for i = 1:length(taus)
    taup = taus(i);
    c = colors(i,:);
    
    right = [eval(subs(rightminus)) eval(subs(rightplus))];
    righti = imag(right) == 0;
    plot(kk(righti), right(righti), '*'); hold on;
end

legend(handles)
title('MC failure envelope with varying \tau/p')

%% Can top p/sigz be greater the right inflection point?

d = (2*C + 2*pp + q + C*q - (4*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) - pp*q - pp*q^2 + q^2 + q*(2*C + 2*pp + 2*q + 2*C*pp + 2*C*q - (4*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) + (4*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1)^2)/(8*(eps + (q + 1)^(1/2)/q)^2 + 2)^2 + 4*C^2*(eps + (q + 1)^(1/2)/q)^2 - 2*pp*q^2 - 2*pp^2*q + C^2 + 4*pp^2*(eps + (q + 1)^(1/2)/q)^2 + pp^2 + q^2 + pp^2*q^2 + (16*(eps + (q + 1)^(1/2)/q)^2*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1)^2)/(8*(eps + (q + 1)^(1/2)/q)^2 + 2)^2 - 4*pp*q^2*(eps + (q + 1)^(1/2)/q)^2 - 8*pp^2*q*(eps + (q + 1)^(1/2)/q)^2 - (4*C*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) + 4*pp^2*q^2*(eps + (q + 1)^(1/2)/q)^2 - 2*C*pp*q - (4*pp*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) - (4*q*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) + 8*C*pp*(eps + (q + 1)^(1/2)/q)^2 + 4*C*q*(eps + (q + 1)^(1/2)/q)^2 + 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + (4*pp*q*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) - 8*C*pp*q*(eps + (q + 1)^(1/2)/q)^2 - (16*pp*(eps + (q + 1)^(1/2)/q)^2*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) - (8*q*(eps + (q + 1)^(1/2)/q)^2*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) - (16*C*(eps + (q + 1)^(1/2)/q)^2*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) + (16*pp*q*(eps + (q + 1)^(1/2)/q)^2*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2) + 1)^(1/2) - (2*q*(C + pp + q + 4*C*(eps + (q + 1)^(1/2)/q)^2 - pp*q + 4*pp*(eps + (q + 1)^(1/2)/q)^2 + 2*q*(eps + (q + 1)^(1/2)/q)^2 - 2*(eps + (q + 1)^(1/2)/q)*(q^2*(eps + (q + 1)^(1/2)/q)^2 - q - 1)^(1/2) - 4*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1))/(8*(eps + (q + 1)^(1/2)/q)^2 + 2))/(2*q - 2*q^2*(eps + (q + 1)^(1/2)/q)^2 + 2) + (q*(2*C*pp - 2*pp - 2*q - 2*C + 2*C*q - 4*C*(eps + (q + 1)^(1/2)/q)^2 + 4*pp*q - 4*pp*(eps + (q + 1)^(1/2)/q)^2 - 4*q*(eps + (q + 1)^(1/2)/q)^2 + 4*C^2*(eps + (q + 1)^(1/2)/q)^2 - 2*pp*q^2 - 2*pp^2*q + C^2 + 4*pp^2*(eps + (q + 1)^(1/2)/q)^2 + pp^2 + q^2 + pp^2*q^2 - 4*pp*q^2*(eps + (q + 1)^(1/2)/q)^2 - 8*pp^2*q*(eps + (q + 1)^(1/2)/q)^2 + 4*pp^2*q^2*(eps + (q + 1)^(1/2)/q)^2 - 2*C*pp*q + 8*C*pp*(eps + (q + 1)^(1/2)/q)^2 + 4*C*q*(eps + (q + 1)^(1/2)/q)^2 + 8*pp*q*(eps + (q + 1)^(1/2)/q)^2 - 8*C*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1)^(1/2) - pp - C + C*q + (2*C*pp - 2*pp - 2*q - 2*C + 2*C*q - 4*C*(eps + (q + 1)^(1/2)/q)^2 + 4*pp*q - 4*pp*(eps + (q + 1)^(1/2)/q)^2 - 4*q*(eps + (q + 1)^(1/2)/q)^2 + 4*C^2*(eps + (q + 1)^(1/2)/q)^2 - 2*pp*q^2 - 2*pp^2*q + C^2 + 4*pp^2*(eps + (q + 1)^(1/2)/q)^2 + pp^2 + q^2 + pp^2*q^2 - 4*pp*q^2*(eps + (q + 1)^(1/2)/q)^2 - 8*pp^2*q*(eps + (q + 1)^(1/2)/q)^2 + 4*pp^2*q^2*(eps + (q + 1)^(1/2)/q)^2 - 2*C*pp*q + 8*C*pp*(eps + (q + 1)^(1/2)/q)^2 + 4*C*q*(eps + (q + 1)^(1/2)/q)^2 + 8*pp*q*(eps + (q + 1)^(1/2)/q)^2 - 8*C*pp*q*(eps + (q + 1)^(1/2)/q)^2 + 1)^(1/2) + 2*pp*q - pp*q^2 + q^2 + 1)/(2*q + 4*q*(eps + (q + 1)^(1/2)/q)^2 + 2*(eps + (q + 1)^(1/2)/q)^2 + 2*q^2*(eps + (q + 1)^(1/2)/q)^2);
eps = 0.00001;
phi = 35/180*pi;
C = 1; % this is C0/sigz
q = tan(pi/4 + 1/2*phi)^2;
pp = 1/2.7;
eval(subs(d))


