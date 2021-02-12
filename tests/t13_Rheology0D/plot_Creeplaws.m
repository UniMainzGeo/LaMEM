% This takes a LaMEM input script with complex (multiple) rheological creeplaws & 
% runs it for different strainrates. It reads out the stress & average
% viscosities
%
% Change the creep law within Rheology_CombinedCreepLaws_0D.dat, to test
% different ones
clear


Strainrate_vec = -22:.5:-12;

T_vec = 400:50:900;
T_vec = 650;
system('mkdir ./Out_DisCr')
for j=1:length(T_vec)
    T = T_vec(j);
    for i=1:length(Strainrate_vec)
        
        Ebg         = Strainrate_vec(i);     % Get strain rate
        
        % Command to be run
        command = ['../../bin/opt/LaMEM -ParamFile Rheology_CombinedCreepLaws_0D.dat -exx_strain_rates  -',num2str(10^Ebg)];
        
        % Add temperature (C):
        command = [command, ' -cstTemp[0] ',num2str(T),' -temp_top ',num2str(T),' -temp_bot ',num2str(T)];
        
        % Run LaMEM:
        system(command)
        
        % run python command, that computes average T2nd, viscosity etc. and
        % writes them into a file
        command =   '/Users/kausb/anaconda3/bin/python TransferData_MATLAB.py'; system(command);
        out     =   csvread('output.txt');        % read python output
        
        T2nd(i,j)     =   out(1); % stress (MPa)
        E2nd(i,j)     =   out(2); % Strainrate
        Temp(i,j)     =   out(3); % T [C]
        P(i,j)        =   out(4); % Pressure (total) [MPa]
        log10Eta(i,j) =   out(5); % Creep Viscosity
        
    end
    
end


if length(T_vec)==1
    
    figure(1), clf
    subplot(121)
    loglog(E2nd,T2nd,'o-')
    
    xlabel('Strainrate [1/s]')
    ylabel('Stress [MPa]')
    title(['T = ',num2str(T),' C'])
    
    
    subplot(122)
    loglog(E2nd,10.^log10Eta,'o-')
    
    xlabel('Strainrate [1/s]')
    ylabel('Eta [Pas]')
    
    
else
    figure(1), clf
    subplot(121)
    contourf(Temp, log10(E2nd),log10(T2nd),25);
    colorbar
    title('log10(T2nd)')
    
    subplot(122)
    contourf(Temp, log10(E2nd),log10Eta,25);
    colorbar
    title('log10(Eta)')
    
    
    
    
end


if 1==1
    % compute this in an independent manner
    eII         = E2nd;
    flow_choice = 1;                % 1- Diff; 2-Disl; 9-Diffusion + Dislocation Creep Assumed
    
    PPa         = 300*1e6*0;          % Pressure in general in Pa
%     pphase      = ['an_wet'];
    pphase      = ['an_dry'];
    
    gsiz        = 100;
    TK          = Temp+273.15;
    
    
    
    [mu_d]        = fun_visc(1,eII,flow_choice,gsiz,TK,PPa,pphase);
    Tau_anal    = 2.*mu_d.*eII/1e6;   % in MPa

    subplot(121)
    hold on
    loglog(eII,Tau_anal,'+','MarkerSize',10)
  
    subplot(122)
    hold on
    loglog(eII,mu_d,'+','MarkerSize',10)
    
    flow_choice = 2;                % 1- Diff; 2-Disl; 9-Diffusion 
    [mu_n]        = fun_visc(1,eII,flow_choice,gsiz,TK,PPa,pphase);
    Tau_anal    = 2.*mu_n.*eII/1e6;   % in MPa
    
    
    subplot(121)
    hold on
    h=loglog(eII,Tau_anal,'ks','MarkerSize',10)
  
    subplot(122)
    hold on
    loglog(eII,mu_n,'ks','MarkerSize',10)
  
    
    legend('LaMEM','Matlab routine diffusion creep', 'Matlab routine dislocation creep')
%     if  flow_choice == 1  
%         title('diffusion creep')
%     elseif  flow_choice == 2
%         title('dislocation creep')
%     elseif  flow_choice == 9
%         title('harmonic average of dffusion & dislocation creep')
%     end
        
        
        
    
    
end




