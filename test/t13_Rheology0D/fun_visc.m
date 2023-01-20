function [mu]  = fun_visc(input_choice,input_var,flow_choice,gsiz,TK,PPa,string_phs)
%
% function to calculate effective viscosity from lab experiments
% EXPERIMENTS = UNIAXIAL SHORTENING
% --------------------------------------------------------------
% phases included:
%===================================
% qtz_wet (Rutter & Brodie 2004)
% an_wet  (Rybacki et al 2006)
% an_dry  (Rybacki et al 2006)
% ol_wet  (Karato & Jung 2003)
% ol_dry  (Karato & Jung 2003)
%===================================
%input_choice= 1; %For Strain Rate (1/s)
%input_choice= 2; %For Stress      (MPa)
%input_var        %Stress in MPa or strain rate in 1/s depending on choice
%flow_choice = 1; %For Diffusion Creep
%flow_choice = 2; %For Dislocation Creep
%flow_choice = 9; %For Diffusion + Dislocation Creep - Harmonic Average
%
% %Example Input
% input_choice= 1;           %Viscosity as a function of strain rate
% input_var   = 1e-12;       %Strain Rate
% flow_choice = 9;           %Diffusion + Diflocation Creep Assumed
% gsiz        = 100;         %grain size in microns %
% TK          = 1000;        %Temperature K
% PPa         = 100;         %Pressure Pa
% string_phs  = ['qtz_wet']; %for wet quartz rheology
%[mu]  = fun_visc(input_choice,input_var,flow_choice,gsiz,TK,PPa)
if strcmp(string_phs,'qtz_wet')==1
    % Data table from Burgman & Dresen 2008 - Experiment after Rutter and Brodie 2004 ---
    %index   dif   disl
    logA   = [-0.4  -4.9]; %Logarithm of pre-exponential factor
    npow   = [   1     3]; %Power law exponent Rutter Brodie give 2.97
    Qact   = [ 220   242]; %Activation energy (KJ)
    m_gr0  = [   2     0]; %Grain size Exponent (will convert to negative)
    r_fug  = [   0     1]; %Exponent of Fugacity
    Vact   = [   0     0]; %Activation Volume cm-3
    fugH   = [ 300   300]; %Fugacity of water MPa (from Experiment RB04)
elseif strcmp(string_phs,'an_wet')==1
    % Data table from Burgman & Dresen 2008 - Experiment after Rybacki et al 2006 ---
    %index   dif   disl
    logA   = [-0.7   0.2]; %Logarithm of pre-exponential factor
    npow   = [   1     3]; %Power law exponent
    Qact   = [ 159   345]; %Activation energy (KJ)
    m_gr0  = [   3     0]; %Grain size Exponent (will convert to negative)
    r_fug  = [   1     1]; %Exponent of Fugacity
    Vact   = [  38    38]; %Activation Volume cm-3
    fugH   = [ 10^2.2 10^2.2]; %Fugacity of water MPa (Actually they do not provide - Average Value - their figure 3)
elseif strcmp(string_phs,'an_dry')==1
    % Data table from Burgman & Dresen 2008 - Experiment after Rybacki et al 2006 ---
    %index   dif   disl  
    logA   = [12.1  12.7]; %Logarithm of pre-exponential factor
    npow   = [   1     3]; %Power law exponent
    Qact   = [ 460   641]; %Activation energy (KJ)
    m_gr0  = [   3     0]; %Grain size Exponent (will convert to negative)
    r_fug  = [   0     0]; %Exponent of Fugacity
    Vact   = [  24    24]; %Activation Volume cm-3
    fugH   = [  1      1]; %Fugacity of water MPa 
    %(Dry samples are assumed to have fugacity exponent r_fug = 0, cf Rybacki et al 2006)
elseif strcmp(string_phs,'ol_dry')==1
    % Data table from Burgman & Dresen 2008 - Exper after Karato & Jung 2003 (disl)
    %index   dif   disl                    and Hirth & Kohlstedt (diff)
    logA   = [9.2    6.1]; %Logarithm of pre-exponential factor
    npow   = [   1     3]; %Power law exponent
    Qact   = [ 375   510]; %Activation energy (KJ)
    m_gr0  = [   3     0]; %Grain size Exponent (will convert to negative)
    r_fug  = [   0     0]; %Exponent of Fugacity
    Vact   = [  10    14]; %Activation Volume cm-3
    fugH   = [   1     1]; %Fugacity of water MPa
elseif strcmp(string_phs,'ol_wet')==1
    % Data table from Burgman & Dresen 2008 - Exper after Karato & Jung 2003 (disl)
    %index   dif   disl                       and Hirth & Kohlstedt (diff)
    logA   = [7.4    2.9]; %Logarithm of pre-exponential factor
    npow   = [   1     3]; %Power law exponent
    Qact   = [ 375   470]; %Activation energy (KJ)
    m_gr0  = [   3     0]; %Grain size Exponent (will convert to negative)
    r_fug  = [   1   1.2]; %Exponent of Fugacity
    Vact   = [  20    24]; %Activation Volume cm-3
    fugH   = [  1      1]; %Fugacity of water MPa (An average value from Karato book- rigorously not correct - avoid)
else
    error('Phase not found, check spelling')
end
% Conversion Factors and constants ---------------------------------------------------
R      = 8.314; %Gas Constant
MPa2Pa = 1e6;   %MPa  -> Pa
cm32m3 = 1e-6;  %cm3  -> m3
J2kJ   = 1e-3;  %Joul -> kJoule
%Choice of Flow law
if flow_choice == 1
    i_flow = 1;    %Diffusion or dislocation creep
    i_run  = 1;    %Run the whole thing once
elseif flow_choice ==2
    i_flow = 2;    %Choose the second flowlaw
    i_run  = 1;    %Run the whole thing once
elseif flow_choice ==9
    i_run  = 2;    %Run The whole thing and then harmonic average
end
for ii = 1:i_run
    if flow_choice ==9 & ii==1
        i_flow = 1;
    elseif flow_choice ==9 & ii==2
        i_flow = 2;
    end
    % Preprocess
    A0     = 10.^(logA);
    m_gr   = -m_gr0;
    PMPa   =  PPa/MPa2Pa;
    %Geometrical Factor for PURE SHEAR (0.5 if want to check vs experiments i.e. use sigma_d -
    %   i.e. no gemometric factor correction depends on the setup of the experiment)
    FG_e   = 1./(2.^((npow(i_flow)-1)./npow(i_flow)).*3.^((npow(i_flow)+1)./(2.*npow(i_flow))))
    FG_s   = 1./(3.^((npow(i_flow)+1)./2));    
    if input_choice ==1
        eII    = input_var;
        mu1    =         FG_e.*eII.^(1/npow(i_flow)-1).*A0(i_flow).^(-1           /npow(i_flow))...
            .*gsiz.^(-m_gr(i_flow)./npow(i_flow)).*fugH(i_flow).^(-r_fug(i_flow)/npow(i_flow))...
            .*exp((Qact(i_flow)+PMPa.*MPa2Pa.*Vact(i_flow).*cm32m3.*J2kJ)./(R.*J2kJ.*TK.*npow(i_flow)));
    end    
    if input_choice ==2
        sII    = input_var;
        mu1    =        FG_s.*sII.^(1-npow(i_flow)).*A0(i_flow).^(-1).*fugH(i_flow).^(-r_fug(i_flow))...
            .*gsiz.^(-m_gr( i_flow))...
            .*exp((Qact(i_flow)+PMPa.*MPa2Pa.*Vact(i_flow).*cm32m3.*J2kJ)./(R.*J2kJ.*TK));
    end    
    if flow_choice ==9
        if ii == 1
            mu_temp = mu1;
        else
            mu = 1./(1./mu_temp + 1./mu1).*MPa2Pa; %In Pa.s
            %mu = min(mu_temp,mu1).*MPa2Pa; %Uncomment this if you do not
            %                               %want harmonic averaging
        end
    else %if Flow Choice 1 or 2
            mu = mu1.*MPa2Pa;  %In Pa.s
    end
end
return;







% %Examples for Deformation Maps ----------------------------------------
% clear
% clc
% close all
% 
% 
% eII         = 1e-12;
% flow_choice = 9;         %Diffusion + Dislocation Creep Assumed
% PPa         = 300*1e6;   %Pressure in general in Pa
% 
% TK1 = linspace(773,1273);
% GS1 = logspace(0,4);
% 
% [gsiz TK] = ndgrid(GS1,TK1);
% 
% [mu]  = fun_visc(1,eII,flow_choice,gsiz,TK,PPa);
% 
% stress = 2*mu*eII*1e-6; %In MPa
% 
% figure(1)
% pcolor(log10(gsiz),log10(stress),TK-273),shading flat,colorbar
% xlabel('log10 of grain size in microns')
% ylabel('log10 of stress in MPa')
% title('Temperature C')
% 
% figure(2)
% pcolor(log10(gsiz),log10(stress),log10(mu)),shading flat,colorbar
% xlabel('log10 of grain size in microns')
% ylabel('log10 of stress in MPa')
% title('Viscosity log10')
% 
% figure(3)
% plot(log10(mu(:)),TK(:),'.')
% xlabel('log10(visc)')
% ylabel('TK')
% 
% %Check sensitivity with Strain Rate
% eII1 = logspace(-15,-9);
% gsiz = 10;
% [eII TK] = ndgrid(eII1,TK1);
% [mu]  = fun_visc(1,eII,flow_choice,gsiz,TK,PPa);
% 
% figure(4)
% sMPa  = 2*mu.*eII*1e-6;
% 
% pcolor(log10(mu),TK,log10(sMPa)),shading flat,colorbar
% xlabel('log10 of viscosity')
% title('log10 of stress in MPa')
% ylabel('Temperature (K)')
% 
% figure(5)
% pcolor(TK,log10(sMPa),log10(mu)),shading flat,colorbar
% title('log10 of viscosity')
% ylabel('log10 of stress in MPa')
% xlabel('Temperature (K)')
% 
% figure(6)
% pcolor(log10(eII),log10(sMPa),log10(mu)),shading flat,colorbar
% xlabel('log10 of strain rate')
% ylabel('log10 of stress in MPa')
% title('log10 of viscosity')
