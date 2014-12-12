function Srr = ElasticTensionShooting(S_0)

E = 1e7;
tau0 = 5e-3;
sigma = 0.5;
R0 = 0.7;

S_initial = [S_0;0;0];
R_initial = linspace(0, 2, 2001);

options = odeset('Events',@ElasticTensionShootingEvent,'RelTol',1e-4);

[R,S_Stack] = ode23(@ElasticTensionUpdate,R_initial,S_initial,options,sigma);

% create dimensional pressure distribution

dS = diff(S_Stack(:,1))./diff(R);
p = (1./sqrt(1+dcorneal(R(2:end)).^2)).*( -S_Stack(2:end,1).*dcorneal(R(2:end))./R(2:end) ...
    - dS.*dcorneal(R(2:end)) ...
    -S_Stack(2:end,1).*d2corneal(R(2:end)) ...
    +S_Stack(2:end,1).*(dcorneal(R(2:end)).^2).*d2corneal(R(2:end))./(1+(dcorneal(R(2:end))).^2) );

p = (E*tau0/R0)*tau(S_Stack(:,3)).*[ - 2*S_0*d2corneal(0);p];

% create the dimensional radial strain

dR = [(R(2)-R(1))/(S_Stack(2,3)-S_Stack(1,3)); diff(R)./diff(S_Stack(:,3))];
RadialStrain = (abs(dR).*sqrt(1+ (dcorneal(R)).^2) - sqrt(1 + (drest(S_Stack(:,3))).^2))./(sqrt(1 + (drest(S_Stack(:,3))).^2));

% create the dimensional angular tension

AngularTension = E*tau0*tau(S_Stack(:,3)).*(sigma*S_Stack(:,1) + (R - S_Stack(:,3))./S_Stack(:,3));

save ThicknessVariationSmoothCornea_v1.mat R S_Stack S_0 p RadialStrain AngularTension R0 E tau0

Srr = S_Stack(end,1);  % S_Stack(end,1) = S_rr = 0

 