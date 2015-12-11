% Plots convergenve results

clear


% fname{1}   = 'FB_convergence_BlockMG_VC1e0.dat';
% fname{2}   = 'FB_convergence_BlockMG_VC1e1.dat';
% fname{3}   = 'FB_convergence_BlockMG_VC1e2.dat';
% fname{4}   = 'FB_convergence_BlockMG_VC1e3.dat';
% fname{5}   = 'FB_convergence_BlockMG_VC1e4.dat';
% fname{6}   = 'FB_convergence_BlockMG_VC1e5.dat';

% fname{1}   = 'FB_convergence_CoupledMG_VC1e0.dat';
% fname{2}   = 'FB_convergence_CoupledMG_VC1e1.dat';
% fname{3}   = 'FB_convergence_CoupledMG_VC1e2.dat';
% fname{4}   = 'FB_convergence_CoupledMG_VC1e3.dat';
% fname{5}   = 'FB_convergence_CoupledMG_VC1e4.dat';
% fname{6}   = 'FB_convergence_CoupledMG_VC1e5.dat';



% % Test the GMG method we use to solve Au=f as part of the block factorization method
% % (with GAMG coarse grid solver, and chebyshev/jacobi(3,3) smoothers),
% fname{1}   = 'A_solve_BlockMG_64by64by64_3MGLevels.dat';
% fname{2}   = 'A_solve_BlockMG_128by128by128_4MGLevels.dat';
% fname{3}   = 'A_solve_BlockMG_256by256by256_5MGLevels.dat';        % using 5 GMG levels

% Test the GMG method we use to solve Au=f as part of the fully coupled method
% (with GAMG coarse grid solver, and richardson/jacobi(3,3)/0.5 smoothers),
fname{1}   = 'A_solve_CoupledMG_64by64by64_3MGLevels.dat';
fname{2}   = 'A_solve_CoupledMG_128by128by128_4MGLevels.dat';
fname{3}   = 'A_solve_CoupledMG_256by256by256_5MGLevels.dat';        % using 5 GMG levels




for num=1:length(fname)
    clear a
    % Import data
    a       =   importdata(fname{num},' ',1);
    for i=2:size(a.textdata,1)
        Residual.Iter(i-1)                  =   str2num(a.textdata{i,1});
        Residual.UnpreconditionedNorm(i-1)	=   str2num(a.textdata{i,6});
        Residual.TrueNorm(i-1)             	=   str2num(a.textdata{i,10});
        Residual.Relative(i-1)              =   a.data(i-1);
    end
    
    ResidualResults{num} = Residual;
end


lnwidth=2;


figure(1), clf
syms = {'k','r','b','g','m','y'}

for i=1:length(fname)
    
    
%     h=semilogy( ResidualResults{i}.Iter,ResidualResults{i}.Relative,syms{i} )
     h=semilogy( ResidualResults{i}.Iter,ResidualResults{i}.UnpreconditionedNorm/ResidualResults{i}.UnpreconditionedNorm(1),syms{i} )
     
    set(h,'linewidth',lnwidth)
    hold on
    
end

xlabel('# iterations','fontsize',20)
ylabel('Rel. Residual','fontsize',20)
% title('Block MG (uncoupled) - Falling block test','fontsize',20)

% title('Coupled MG - Falling block test','fontsize',20)

axis([0 100 1e-6 1e1])

% h=legend('Viscosity contrast 1','Viscosity contrast 10','Viscosity contrast 10^2','Viscosity contrast 10^3','Viscosity contrast 10^4','Viscosity contrast 10^5','Location','Southeast')
% set(h,'fontsize',12)


% Title and legend for Ku=f test
title('Scalability of solving Ku=f (velocity block)','fontsize',20)

h=legend('velocity\_DOF\_0.8M','velocity\_DOF\_6.3M','velocity\_DOF\_50.5M','Location','Northeast')
set(h,'fontsize',12)
set(gca,'fontsize',12)
text(0,4e-7,'Solver details: GMRES, with PC=GMG with Chebyshev/Jacobi(3,3) smoothers and GAMG as coarse grid solver','fontsize',8)
text(40,2e-6,'Multiple falling spheres setup with viscosity matrix=1 and sphere=1e3','fontsize',8)
 




% Title and legend for Ku=f in coupled MG problem
title('Scalability of solving |K VP; PV 0| |u; p|=|f; g| (coupled MG)','fontsize',20)

h=legend('velocity\_DOF\_0.8M','velocity\_DOF\_6.3M','velocity\_DOF\_50.5M','Location','Northeast')
set(h,'fontsize',12)
set(gca,'fontsize',12)
text(0,4e-7,'Solver details: GMRES, with PC=GMG with Chebyshev/Jacobi(3,3) smoothers and GMRES/GAMG as coarse grid solver','fontsize',8)
text(40,2e-6,'Multiple falling spheres setup with viscosity matrix=1 and sphere=1e3','fontsize',8)
 