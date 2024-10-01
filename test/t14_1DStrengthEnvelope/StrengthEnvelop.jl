using Dates

function Plot_StrengthEnvelop(filename::String, dir::String, z::Vector, DataSets::NTuple{N,Vector}, Names::NTuple{N,String}) where {N}
    # Open figure 
    f = Figure(size = (800, 1500))
    datetime = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    ax = Axis(f[1, 1],  xlabel = "τII [MPa]", ylabel = "Depth [km]", title=datetime)
    
    for i=1:length(DataSets)
        τII = DataSets[i]
        lines!(ax, τII, z,  label = Names[i]) 
    end

    axislegend(position=:rb)

    save(joinpath(dir,filename), f) 
    return f 
end 


function Analytical_StrengthEnvelop(phase, T, P, τy)
    # list of parameters for analytical solution
    # Phases: 0: Air, 1: WetQuartzite, 2: Granite, 3: DryOlivine
    Bn      = [5e-19, 1.55371e-17, 1.67675e-25, 1.48058e-16]
    En      = [0, 154e3, 186.5e3, 532e3]
    n       = [1, 2.3, 3.3, 3.5]
    Vn      = [0,0,0,1.7e-5]
    R       = 8.314463          # gas constant
    ε       = 1e-15             # strainrate

    I       = round.(Int, phase) .+ 1       # phase of rock
    T_K     = T .+ 273.15;                  # in K
    τ       = 1e-6*Bn[I].^(-1.0 ./ n[I]) .* ε.^(1.0./n[I]) .* exp.( (En[I] + Vn[I].*P)./(n[I].*R.*T_K ) )    # in MPa
    ind     = findall( τ .> τy)
    τ[ind]  = τy[ind]

    return τ
end