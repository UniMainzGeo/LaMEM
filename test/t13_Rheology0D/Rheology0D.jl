using Statistics

"""
    t, τII =  StressTime_0D(FileName::String, DirName::String="")

This extracts τII-time curves from LaMEM results, by computing the average stress in the domain
"""
function StressTime_0D(FileName::String, DirName::String="")

    Timestep, _, _= Read_LaMEM_simulation(FileName, DirName);
    
    N = length(Timestep)
    τII_vec = Vector{Float64}(undef, N)
    t_vec = Vector{Float64}(undef, N)
    for (i,it) in enumerate(Timestep)
        data,t = read_LaMEM_timestep(FileName, it, DirName, fields=("j2_dev_stress [MPa]",));
        τII =    Float64(mean(data.fields.j2_dev_stress))
        τII_vec[i] = τII
        t_vec[i] = t[1]
    end

    return t_vec, τII_vec
end

"""
    τII = Viscoelastoplastic0D(G, η, ε, t, τ_yield=1e100)

Analytical stress buildup in 0D, for a given elastic shear module `G`, viscosity `η`, strainrate `ε` and timesteps `t`.
The returned stress is in [Pa]
"""
function Viscoelastoplastic0D(G, η, ε, t, τ_yield=1e100)
    
    SecYear   = 3600*24*365.25;   # sec/year
    τ_maxwell = η/G;              # Maxwell time
    time      = t*SecYear*1e6;    # t in seconds

    τII = 2*η*(1.0 .- exp.(-time/τ_maxwell))*ε
    τII[τII .> τ_yield] .= τ_yield # plasticity

    return τII # in Pa
end



function Plot_StressStrain(t_anal,τII_anal, t_num, τII_num, dir, filename="Analytics_vs_LaMEM.png"; τII_no_iter=nothing)

    # Open figure 
    f = Figure(resolution = (1500, 800))
    ax = Axis(f[1, 1],  xlabel = "time [Myrs]", ylabel = "τII [MPa]")
    lines!(ax, t_anal, τII_anal,  label = "Analytical") 
    if !isnothing(τII_no_iter)
        lines!(ax, t_anal, τII_no_iter,  label = "Analytical, no iterations") 
    end
    scatter!(ax, t_num, τII_num,  label = "LaMEM", color=:red) 
    

    axislegend()

    save(joinpath(dir,filename), f) 

    return f 
end 


function Viscoelastoplastic0D_dislocationcreep(T, ε, tmax,  τ_yield=1e100, G = 5e10, AD = 2.5e-17, n = 3.5, Ea = 532e3)
    
    SecYear   = 3600*24*365.25;   # sec/year
    t         = range(0,tmax, 200);
    t_s       = t*SecYear*1e6

    # Define material properties for each of the phases 
    
    # these are the parameters for Dry Olivine, to be consistent with Gerya's book; can be changed ofcourse
    R     = 8.3144621;     # J/mol/K

    # Define correction coefficient F2 
    # for a strain rate based viscosity formulation
    F2          = 1/2^((n-1)/n)/3^((n+1)/2/n);
    eterm       = exp(Ea/n/R/(T+273.15));

    τII = zeros(Float64,length(t))
    for i = 2:length(t)
        Δt = t_s[i] - t_s[i-1]

        # to do this correct, we need to perform local iterations
        τ_new = LocalIterations_DC(τII[i-1], Δt, G, ε, n, eterm, F2, AD, τ_yield)

        τII[i] = τ_new
    end

    η     = F2/AD^(1/n)/ε^((n-1)/n)*eterm;     # effective viscosity
    τ_M   = η/G;

    τII_noLocal  = 2*η*(1.0 .- exp.(-t_s/τ_M))*ε; 
    τII_noLocal[τII_noLocal .> τ_yield] .= τ_yield # plasticity

    return t, τII, τII_noLocal  # in Pa
end

# perform local no linear iterations
function LocalIterations_DC(τ, Δt, G, ε, n, eterm, F2, AD, YieldStress)
    
    ε_vis = ε/2
    τ_new = τ
    for it = 1:50
        τ_new1  =   τ + 2*G*Δt*(ε-ε_vis); 
        η       =   F2/AD^(1/n)/ε_vis^((n-1)/n)*eterm;
        ε_vis   =   τ_new1/2/η;     
        if η>1e28
            η = 1e28
        end
        if τ_new1>YieldStress
            τ_new1 = YieldStress
        end
        dTau    = τ_new1-τ_new;
        τ_new   = τ_new1;
    end

    return τ_new
end


"""
    τ = StressStrainrate0D_LaMEM(FileName, DirName, OutFile, ε_vec)

Computes stress-strainrate curves by running a LaMEM for different strain rates and reading the stress `τ` data back.
Input:
- `FileName`: LaMEM input *.dat file 
- `DirName`: Directory where the input file is located
- `OutFile`: Name of the output file (without the `*.pvd`)
- `ε_vec`: vector with strain rate values

Output:
- `τ`: vector with deviatoric stress values

"""
function StressStrainrate0D_LaMEM(FileName::String, DirName::String="", OutFile="Rheolog0D_linearViscous", ε_vec=[-1e-15 -1e-14])
    cur_dir = pwd();
    τ = zeros(length(ε_vec))
    cd(DirName)
    for (i,ε) in enumerate(ε_vec)
        out = run_lamem_local_test(FileName, 1, "-exx_strain_rates $ε"; opt=true, bin_dir="../../bin")  # run LaMEM
        data, t = read_LaMEM_timestep(OutFile, 2, pwd(), fields=("j2_dev_stress [MPa]",))   # read stress
        τ[i] = mean(data.fields.j2_dev_stress)  # store
    end
    cd(cur_dir)
    return τ
end


"""
"""
function Plot_StressStrainrate(ε, τ, τ_anal, dir, filename="t13_Stress_Strainrate.png")

    # Open figure 
    f = Figure(resolution = (1500, 800))
    ax = Axis(f[1, 1],  xlabel = "strainrate [1/s]", ylabel = "τII [MPa]"; yscale=log10, xscale=log10)
    lines!(ax, -ε[:], τ_anal[:],  label = "Analytical") 

    scatter!(ax, -ε[:], τ[:],  label = "LaMEM", color=:red) 
    

    axislegend(position=:rb)

    save(joinpath(dir,filename), f) 

    return f 
end 


function AnalyticalSolution_DislocationCreep(type="DryOlivine", T=1000, ε=[1e-15, 1e-16])

    if type=="DryOlivine"
        AD    = 2.5e-17;   # 1/Pa^n/s, 
        n     = 3.5;       # dimensionless
        Ea    = 532000;    # J/mol 
        R     = 8.3144621;     # J/mol/K

        # Define correction coefficient F2 
        # for a strain rate based viscosity formulation
        F2    = 1/2^((n-1)/n)/3^((n+1)/2/n);
    else
        error("not implemented")
    end

    eterm   =  exp(Ea/n/R/(T+273.15));
    η       =  F2/AD^(1/n)./abs.(ε).^((n-1)/n)*eterm
    τ       =  2.0.*η.*abs.(ε);

    return τ
end