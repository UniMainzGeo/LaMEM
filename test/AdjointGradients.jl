# This tests various adjoint gradient cases
dir = "t08_AdjointGradients";
if test_superlu
   # t8_Adjoint_rho_SensitivityKernel
   keywords   = (  "|Div|_inf",
                  "|Div|_2",
                  "|mRes|_2")

   acc        = (  (rtol=1e-7, atol=1e-6), 
                  (rtol=1e-5, atol=1e-5), 
                  (rtol=1e-5, atol=1e-5), 
               );

   # Perform tests

   ParamFile = "t8_AdjointGradients.dat";
   @test perform_lamem_test(dir,ParamFile,"t8_AdjointGradients_Sphere_ND_all.expected",
                           args="",
                           keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec)
end
if test_superlu
   # t8_AdjointGradients_Sphere_ND_all
   keywords   = (  "|Div|_inf",
                  "|Div|_2",
                  "|mRes|_2",
                  "|           delta(rho)[  1]",
                  "|                  eta[  0]",
                  "|   Velocity check            :",
                  "|  adjoint     2:          eta[ 0]")

   acc        = (  (rtol=1e-7, atol=1e-6), 
                  (rtol=1e-8, atol=1e-5), 
                  (rtol=1e-8, atol=1e-5), 
                  (rtol=1e-6, atol=1e-5), 
                  (rtol=1e-6, atol=1e-5), 
                  (rtol=1e-6, atol=1e-5), 
                  (rtol=1e-6, atol=1e-5), 
               );

   ParamFile = "t8_AdjointGradients.dat";
   @test perform_lamem_test(dir,ParamFile,"t8_AdjointGradients_Sphere_ND_all.expected",
                           args="",
                           keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec)
end

if test_superlu
   # t8_AdjointGradients_CompareGradients_1
   keywords   = (  "|Div|_inf",
                  "|Div|_2",
                  "|mRes|_2",
                  "|       FD     1:          eta[ 1]",
                  "|  adjoint     2:          eta[ 1]",
                  "|       FD     3:          eta[ 0]",
                  "|  adjoint     4:          eta[ 0]")

   acc        = (  (rtol=1e-7, atol=1e-6), 
                  (rtol=1e-8, atol=1e-5), 
                  (rtol=1e-8, atol=1e-5), 
                  (rtol=1e-8, atol=1e-5), 
                  (rtol=1e-8, atol=1e-5), 
                  (rtol=1e-8, atol=1e-5), 
                  (rtol=1e-8, atol=1e-5), 
               );


   ParamFile = "t8_AdjointGradients_CompareGradients.dat";
   @test perform_lamem_test(dir,ParamFile,"t8_AdjointGradients_CompareGradients_1.expected",
                           args="",
                           keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec)
end
# t8_AdjointGradients_CompareGradients_geo
keywords   = (  "|       FD     1:          eta[ 1]",
                "|  adjoint     2:          eta[ 1]",
                "|       FD     3:          eta[ 0]",
                "|  adjoint     4:          eta[ 0]",
                "|           delta(rho)[  1]")

acc        = (  (atol=1e-30, ), 
                (atol=1e-30, ), 
                (atol=1e-28, ), 
                (atol=1e-28, ), 
                (atol=1e-3,  ),  
             );

ParamFile = "t8_AdjointGradients_CompareGradients_geo.dat";
@test perform_lamem_test(dir,ParamFile,"t8_AdjointGradients_CompareGradients_geo.expected",
                        args="",
                        keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec)



#t8_AdjointGradients_CompareGradients_2
keywords   = (  "|Div|_inf",
                "|Div|_2",
                "|mRes|_2",
                "|       FD     1:            n[ 0]",
                "|  adjoint     2:            n[ 0]",
                "|   Prefactor A               :")

acc        = (  (rtol=1e-7, atol=1e-6), 
                (rtol=1e-8, atol=1e-5), 
                (rtol=1e-8, atol=1e-5), 
                (rtol=1e-8, atol=1e-5), 
                (rtol=2e-6, atol=1e-5), 
                (rtol=1e-8, atol=1e-5), 
             );

ParamFile = "t8_AdjointGradients_CompareGradients_2.dat";
@test perform_lamem_test(dir,ParamFile,"t8_AdjointGradients_CompareGradients_2.expected",
                        args="",
                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)

# t8_Adjoint_Subduction2D_FreeSlip
keywords   = (  "|Div|_inf",
                "|Div|_2",
                "|mRes|_2",
                "|                  eta[  0]",
                "|           delta(rho)[  1]",
                "|           delta(rho)[  2]",
                "|                  eta[  2]",
                "|                  eta[  1]",
                "|      log10       eta[  0]",
                "|       FD     7:           fr[ 2]"
                )

acc        = (  (rtol=1e-7, atol=1e-6), 
                (rtol=1e-8, atol=1e-5), 
                (rtol=1e-8, atol=1e-5), 
                (rtol=1e-8, atol=1e-5), 
                (rtol=1e-8, atol=1e-5), 
                (rtol=1e-8, atol=1e-5),
                (rtol=1e-8, atol=1e-5),
                (rtol=1e-8, atol=1e-5),
                (rtol=1e-3, atol=1e-5),
                (rtol=1e-8, atol=1e-5),
             );

ParamFile = "t8_Subduction2D_FreeSlip_DirectSolver.dat";
@test perform_lamem_test(dir,ParamFile,"t8_Subduction2D_FreeSlip_DirectSolver_p1.expected",
                        args="-nel_y 2",
                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)


# t8_Adjoint_PSD
keywords   = (  "|  adjoint     1:          rho[ 2]",
                "|       FD     2:          rho[ 2]",
                "| Current Cost function = "
                )

acc        = (  (rtol=1e-6, atol=1e-6), 
                (rtol=1e-6, atol=1e-5), 
                (rtol=1e-6, atol=1e-5), 
             );

ParamFile = "t8_FB_PSDTest.dat";
@test perform_lamem_test(dir,ParamFile,"t8_FB_PSDTest_p1.expected",
                        args="-nel_x 8 -nel_y 8 -nel_z 8 ",
                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)

# t8_Adjoint_rho_SensitivityKernel_PSD
if test_superlu
   keywords   = ( "|Div|_inf",
                  "|Div|_2",
                  "|mRes|_2",
                  "| Current Cost function = "
                  )

   acc        = (  (rtol=1e-7, atol=1e-6), 
                  (rtol=1e-5, atol=1e-5), 
                  (rtol=1e-4, atol=1e-5), 
                  (rtol=1e-6, atol=1e-5), 
               );   

   ParamFile = "t8_AdjointGradients_SensitivityKernel_PSD.dat";
   @test perform_lamem_test(dir,ParamFile,"t8_Adjoint_rho_SensitivityKernel_PSD_p2.expected",
                           args="",
                           keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec)
end
# t8_Adjoint_n_SensitivityKernel_PSD
keywords   = ( "|   Norm of field gradient vector :",
                )

acc        = (  (rtol=5e-1, atol=1e-6), 
             );
split_sign       = (":",)       

ParamFile = "t8_PSDKernelPaper.dat";
@test perform_lamem_test(dir,ParamFile,"t8_Adjoint_n_SensitivityKernelPaper_PSD.expected",
                        args="-nel_x 8  -nel_y 8 -nel_z 8 ",
                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)

# t8_Adjoint_eta0_SensitivityKernel_PSD
keywords   = ( "|   Norm of field gradient vector :",
                )

acc        = (  (rtol=5e-1, atol=1e-6), 
             );
split_sign       = (":",)       

ParamFile = "t8_PSDKernelPaper.dat";
@test perform_lamem_test(dir,ParamFile,"t8_Adjoint_eta0_SensitivityKernelPaper_PSD.expected",
                        args="-nel_x 8  -nel_y 8 -nel_z 8 -Type[0] eta0",
                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec, 
                        split_sign=split_sign)