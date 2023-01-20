clear,clc,close all

eII         = 1e-13;
flow_choice = 9;                % Diffusion + Dislocation Creep Assumed

PPa         = 300*1e6;          % Pressure in general in Pa
pphase      = ['an_dry'];

TK1 = linspace(773,1273);       % in K
GS1 = logspace(0,4);            % in microns

[gsiz TK]   = ndgrid(GS1,TK1);


[mu]        = fun_visc(1,eII,flow_choice,gsiz,TK,PPa,pphase);


stress = 2*mu*eII*1e-6; %In MPa



figure(1)
pcolor(log10(gsiz),log10(stress),TK-273),shading flat,colorbar
xlabel('log10 of grain size in microns')
ylabel('log10 of stress in MPa')
title('Temperature C')

figure(2)
pcolor(log10(gsiz),log10(stress),log10(mu)),shading flat,colorbar
xlabel('log10 of grain size in microns')
ylabel('log10 of stress in MPa')
title('Viscosity log10')

figure(3)
plot(log10(mu(:)),TK(:),'.')
xlabel('log10(visc)')
ylabel('TK')

%Check sensitivity with Strain Rate
eII1 = logspace(-15,-9);
gsiz = 10;
[eII TK] = ndgrid(eII1,TK1);
[mu]  = fun_visc(1,eII,flow_choice,gsiz,TK,PPa,pphase);

figure(4)
sMPa  = 2*mu.*eII*1e-6;

pcolor(log10(mu),TK,log10(sMPa)),shading flat,colorbar
xlabel('log10 of viscosity')
title('log10 of stress in MPa')
ylabel('Temperature (K)')

figure(5)
pcolor(TK,log10(sMPa),log10(mu)),shading flat,colorbar
title('log10 of viscosity')
ylabel('log10 of stress in MPa')
xlabel('Temperature (K)')

figure(6)
pcolor(log10(eII),log10(sMPa),log10(mu)),shading flat,colorbar
xlabel('log10 of strain rate')
ylabel('log10 of stress in MPa')
title('log10 of viscosity')
