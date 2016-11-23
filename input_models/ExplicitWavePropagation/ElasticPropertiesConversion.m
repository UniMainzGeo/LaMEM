% Input:    density & seismic wave velocuty
% Output:   Elastic properties
clear


% % Salt properties
% rho = 2200          % density
% vp  = 4500          % P wave velocity [m/s]
% vs  = vp/sqrt(3)    % S wave velocity


% % Basement (granite) properties
% rho = 2900          % density
% vp  = 6000          % P wave velocity [m/s]
% vs  = vp/sqrt(3)    % S wave velocity


% Sediment1 (sandstone) properties
rho = 2900          % density
vp  = 5000          % P wave velocity [m/s]
vs  = vp/sqrt(3)    % S wave velocity

% % Sediment2 (sandstone) properties
% rho = 2900          % density
% vp  = 4000          % P wave velocity [m/s]
% vs  = vp/sqrt(3)    % S wave velocity


G       = vs^2*rho                      % shear, G
lamda   = vp^2*rho-2*G;
K       = rho*(vp^2 - (4/3)*vs^2)     % bulk


vp = sqrt((K+(4/3)*G)/rho);  % P wave velocity
vs = sqrt(G/rho);               % S wave velocity