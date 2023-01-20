% 
clear

SecYear =   3600*24*365.25
Tau     =   0;
G       =   5e10;
eta     =   1e22;
T       =   900;

% Flow law parameters
AD      =   2.5e-17; %1/Pa^n/s, 
n       =   3.5;
Ea      =   532000; %J/mol 
R       =   8.314; %J/mol/K

% Define correction coefficient F2 
% for a strain rate based viscosity formulation
F2=1/2^((n-1)/n)/3^((n+1)/2/n)

% Strain Rate
eps=1e-15; %1/s


% Compute and check activation energy exponent
eterm=Ea/n/R/(T+273.15);
if (eterm>100)
    eterm=100;
end
eterm=exp(eterm);
% Compute viscosity
% eta=F2/AD^(1/n)/eps^((n-1)/n)*eterm;

dt = 0.0001*1e6*SecYear;

Tau_vec     = 0;
Time_vec    = 0;
for itime=2:1000
    
    Tau_new =   Tau;
    dTau    =   realmax;
    it = 1;
    while dTau>1e-10
        e_vis       =   Tau_new/2/eta;               % viscous strainrate
        Tau_new1    =   Tau + 2*G*dt*(eps-e_vis);    % update stress
         eta         =   F2/AD^(1/n)/e_vis^((n-1)/n)*eterm;
%         eta         =   F2/AD^(1/n)/eps^((n-1)/n)*eterm;            % no local iterations
        

%       Tau_new1    =   Tau + 2*G*dt*(str-e_vis);               # update stress
%       eta         =   F2/AD**(1/n)/e_vis**((n-1)/n)*eterm;    # local iterations; employ viscous strain rate here (the wrong approach is to use the full strainrate)
%       e_vis       =   Tau_new1/2/eta;                          # viscous strainrate
%       
        if eta>1e28
            eta=1e28;
        end
        
        dTau        = Tau_new1-Tau_new;
        Tau_new     = Tau_new1;
        it=it+1;
        
    end
    it
    Tau             = Tau_new;
    
    Tau_vec(itime)  =   Tau_new;
    Time_vec(itime) =   Time_vec(itime-1) + dt;
end


figure(1), hold on
plot(Time_vec/SecYear/1e6,Tau_vec/1e6,'r-');
xlabel('Time [myrs]')
ylabel('Tau [MPa]')

