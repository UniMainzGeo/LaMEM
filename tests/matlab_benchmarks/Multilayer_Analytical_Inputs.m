function [misfit,lam_vec,q_vec_SaltSediments, l_c] = Multilayer_Analytical_Inputs(mu, Ampl, lam_vec,LowerBC)
% Slight modification of DetachmentFolding_7layersOverburden_InputForOptimization.m of Boris Kaus
% Computational setup with a salt layer and an overburden that is composed of 7 layers
% with given thickness but of different viscosity that call the
% Biot_nlayer_powerlaw_time.m
% Thickness set up is equal to the numerical setup of DetachmentFolding_Zagros_Cos_Amplitude1.dat

mu_1 = mu(1);
mu_2 = mu(2);
mu_3 = mu(3);
mu_4 = mu(4);
mu_5 = mu(5);
mu_6 = mu(6);
mu_7 = mu(7);

Ampl_1 = Ampl(1);
Ampl_2 = Ampl(2);
Ampl_3 = Ampl(3);
Ampl_4 = Ampl(4);
Ampl_5 = Ampl(5);
Ampl_6 = Ampl(6);
Ampl_7 = Ampl(7);
Ampl_8 = Ampl(8);
Ampl_9 = Ampl(9);

n_lay   =   1;
n_mat   =   1;

for ilam=1:length(lam_vec)
    ilam;
    
    %======================================================================
    % Initial parameters in dimensional units
    
    %    Layer                7             6           5           4           3           2           1               (detachment layer)   
    mu              =   [     10^mu_7       10^mu_6     10^mu_5     10^mu_4     10^mu_3     10^mu_2     10^mu_1         1e19     	].';        %   define viscosity structure
    
    rho             =   [     2700          2700     	2700        2700        2700        2700        2700            2200        ].';        %   define density structure
    n               =   [     n_lay         n_lay       n_mat       n_lay      	n_lay       n_lay       n_lay           n_mat       ].';        %   powerlaw exponent
    g               =   10;                                                                                                                                                      %   gravity
    gam             =   zeros(size(mu));
    
    Ampl            =   [   Ampl_9      Ampl_8      Ampl_7      Ampl_6      Ampl_5      Ampl_4       Ampl_3      Ampl_2       Ampl_1	].';        %   define initial amplitudes of interface

    H_int           =   [  6.5e3        6.0e3       5.5e3       5e3         4.5e3       4e3         3.5e3           3e3       	0   ].';        %   mean height of interface
   
    ninterface      =   length(H_int);
    Ind_interface   =   8;                  % The main perturbed interface

    Bound           =   lower({'Free Slip','No Slip','fast erosion','Free surface'});
    lower_bound     =   Bound{LowerBC};
    upper_bound     =   Bound{4};

    str_bg          =   -1e-15;
    Ampl_max        =   1.001;

    k_e             =   0;          % erosional diffusivity
    
    %======================================================================
    
    %======================================================================
    % Perform non-dimensionalization
    
    % Select characteristic scales
    mu_c                            =   1e20;                               % characteristic viscosity
%     l_c                             =   10e3;                               % characteristic length
    l_c = max(H_int);
    t_c                             =   1/1e-15;                            % makes 1e-15 characteristic strainrate
    T_c                             =   873;                                % characteristic temperature

    % derived scales
    vel_char                        =   (l_c/t_c);                              %   Velocity        [m/s   ]
    Stress_char                     =   mu_c/t_c;                               %   Stress/Pressure [Pa=N/m2]
    force_char                      =   Stress_char*l_c^2;                      %   Force,          [N   = kg * m/s2]
    watt_char                       =   force_char*l_c/t_c;                     %   Watt,           [J/s = m2???kg???s-3]
    joule_char                      =   force_char*l_c;                         %   Joule,          [Nm = m2???kg???s-2]
    kg_char                         =   Stress_char*l_c*t_c^2;                  %   kg
    rho_char                        =   kg_char/l_c^3;                          %   kg/m3
    
    % Non-dimensionalization
    mu                              =   mu/mu_c;
    rho                             =   rho/rho_char;
    n                               =   n;
    g                               =   g/(l_c/t_c^2);
    Ampl                            =   Ampl/l_c;
    H_int                           =   H_int/l_c;
    str_bg                          =   str_bg/(1/t_c);
    Ampl_max                        =   Ampl_max/l_c;
    SecYear                         =   3600*24*365.25;
    %
    %======================================================================
    
    lambda          =   lam_vec(ilam);
    
   [q_dom,q, Vel,Amplitude,Time, Vel_Erosion_Time] = Biot_nlayer_powerlaw_time(lambda,ninterface,mu ,n, gam, rho, g, Ampl,H_int,lower_bound,upper_bound, str_bg, Ampl_max, k_e, Ind_interface);
     
    q_vec_SaltSediments(ilam) = -q_dom;
    q_vec_Surface(ilam)       = q(1);
    
end

[q_max,id]  = max(q_vec_SaltSediments);

lam_dom     = lam_vec(id);

wavelength_km = lam_dom*l_c/1e3;

misfit         = abs(wavelength_km-131.625)/5 + abs(q_max-76.3)/20; % error with observed wavelength
