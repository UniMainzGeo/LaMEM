% Elastic properties

vp  = 5000          % P wave velocity
vs  = vp/sqrt(3)    % S wave velocity
rho = 2200          % density

mu    = vs^2*rho;                       % shear, G
lamda = vp^2*rho-2*G;
K     = density*(vp^2 - (4/3)*vs^2)     % bulk



vp = sqrt((K+(4/3)*G)/density)  % P wave velocity
vs = sqrt(mu/rho)               % S wave velocity