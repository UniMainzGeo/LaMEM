%WeakScalingTests_MatVec_JUQUEEN
%
% This summarizes scaling tests performed on JUQUEEN to test the Mat-Vec
% multiplication performance of LaMEM. This test is useful as the code
% cannot scale better than the performance of Mat-Vec products.

clear

fnt_size = 20;
fnd_name = 'arial';

% TEST 1:
% 
% Falling Block test, FDSTAG, 5000 MatVec products

Cores               = [1        8       64      512      1024      2048    4096   8192    16384];
Time_FDSTAG         = [78.4     93.4    95.52   95.56    95.74     95.822  95.829 95.72   94.87]
Time_Q1P0           = [312.9    336.5   349.6   353.6    348.8     353.0   NaN    NaN     NaN]
Time_Q2Pm1          = [1242.9   1465    1574    1569.18  1569.12   1575.1  1579   NaN     NaN]


relTime_FDSTAG      =   Time_FDSTAG(2)./Time_FDSTAG*100;
relTime_Q1P0        =   Time_Q1P0(2)./Time_Q1P0*100;
relTime_Q2Pm1       =   Time_Q2Pm1(2)./Time_Q2Pm1*100;


% Create plot
figure(1), clf
semilogx(Cores,Time_FDSTAG,'o-',...
         Cores,Time_Q1P0,'o-',...
         Cores,Time_Q2Pm1,'o-')
set(gca,'XTick', [1 8 64 512 1024 4096 16832])
xlabel('Number of cores','Fontsize',fnt_size,'FontName',fnd_name)
ylabel('Time for 5000 MatVec products [s]','Fontsize',fnt_size,'FontName',fnd_name)

legend('FDSTAG','Q1P0','Q2Pm1')
title('WeakScalingTest MatVec products, Falling Block setup with 32x32x32 nodes/core','Fontsize',fnt_size,'FontName',fnd_name)
set(gca,'Fontsize',fnt_size,'FontName',fnd_name)


figure(2), clf
semilogx(Cores(2:end),relTime_FDSTAG(2:end),'o-',...
         Cores(2:end),relTime_Q1P0  (2:end),'o-',...
         Cores(2:end),relTime_Q2Pm1 (2:end),'o-')

set(gca,'XTick', [8 64 512 1024 4096 16832])
xlabel('Number of cores','Fontsize',fnt_size,'FontName',fnd_name)
ylabel('Relative efficiency [%]','Fontsize',fnt_size,'FontName',fnd_name)
legend('FDSTAG','Q1P0','Q2Pm1')


title('WeakScalingTest MatVec products, FDSTAG with 32x32x32/core','Fontsize',fnt_size,'FontName',fnd_name)
set(gca,'Fontsize',fnt_size,'FontName',fnd_name)
