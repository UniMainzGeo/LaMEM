%LaMEM_MATLAB_BENCHMARKS
%
% Perform all available MATLAB Benchmarks

clc, clear;


if ismac
    CreatePlot = logical(1);
else
    CreatePlot = logical(0);
end


% TEST 1: RAYLEIGH-TAYLOR BENCHMARK
TestNumber              = 1;
Test_name{TestNumber}   = 'TWO LAYER Rayleigh-Taylor ';
[TotalError_RT, AllowedError_RT] = LaMEM_Benchmark_2D_RayleighTaylor(CreatePlot);
if TotalError_RT<AllowedError_RT
    Test(TestNumber) = 1;    % passed
else
    Test(TestNumber) = 0;    % failed
end

% TEST 2: Single Layer folding benchmark with Newtonian Rheologies
TestNumber              = 2;
Test_name{TestNumber}   = '2D NEWTONIAN SINGLE LAYER FOLDING ';
[TotalError_SF1, AllowedError_SF1 ] = LaMEM_Benchmark_2D_SingleLayer_Folding(CreatePlot, 1, 1);
if TotalError_SF1<AllowedError_SF1
    Test(TestNumber) = 1;    % passed
else
    Test(TestNumber) = 0;    % failed
end


% TEST 3: Single Layer folding benchmark with Newtonian Rheologies
TestNumber              = 3;
Test_name{TestNumber}   = '2D POWERLAW SINGLE LAYER FOLDING ';
[TotalError_SF2, AllowedError_SF2 ] = LaMEM_Benchmark_2D_SingleLayer_Folding(CreatePlot, 3, 4);
if TotalError_SF2<AllowedError_SF2
    Test(TestNumber) = 1;    % passed
else
    Test(TestNumber) = 0;    % failed
end


% TEST 4: Multilayer folding benchmark with Newtonian Rheologies &
% perturbation at the salt-sediment interface
TestNumber              = 4;
Test_name{TestNumber}   = '2D MULTILAYER LAYER FOLDING WITH PERTURBATION AT SALT-SEDIMENT INTERFACE';
[TotalError_MF1, AllowedError_MF1 ] = LaMEM_Benchmark_MultilayerFolding(CreatePlot, 0);
if TotalError_MF1<AllowedError_MF1
    Test(TestNumber) = 1;    % passed
else
    Test(TestNumber) = 0;    % failed
end

% TEST 5: Multilayer folding benchmark with Newtonian Rheologies &
% perturbation at all  interface
TestNumber              = 5;
Test_name{TestNumber}   = '2D MULTILAYER LAYER FOLDING WITH PERTURBATION AT ALL INTERFACES';
[TotalError_MF2, AllowedError_MF2 ] = LaMEM_Benchmark_MultilayerFolding(CreatePlot, 1);
if TotalError_MF2<AllowedError_MF2
    Test(TestNumber) = 1;    % passed
else
    Test(TestNumber) = 0;    % failed
end

% DISPLAY A SUMMARY TO THE SCREEN:
clc
disp(['==============================================='])
for i=1:length(Test)
    str = ['TEST ',num2str(i),' OF ',num2str(length(Test)),'  : ',Test_name{i}];
    if Test(i)==0
        str = [str,' FAILED'];
    else
        str = [str,' PASSED'];
    end
    
    
    disp(str);
    
end
disp(['==============================================='])



%CLEANUP THE DIRECTORY
if isunix
   !rm *.pvd
   !rm -rf ./Timestep_00000*
   !rm -rf ./InitialMesh
   !rm -rf ./InitialParticles   
end



