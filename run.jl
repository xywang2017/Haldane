using PyPlot
using DataFrames
using CSV
include("haldane.jl")
fpath = pwd()

function compute_spectrum()
    df = DataFrame() 

    ϕs = collect(1:15) .// 61
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        q = denominator(ϕ)
        p = numerator(ϕ)
        hd = Haldane()
        fname = joinpath(fpath,"data/p$(p)q$(q).txt")
        constructHaldane(hd;ϕ=ϕ,lk=1,fname=fname)
        df[!,"$(p)_$(q)"] = hd.spectrum[:]
    end

    out_name = joinpath(fpath,"data/q61_chern0.csv")
    CSV.write(out_name,df,delim=",")
end
compute_spectrum()
##
function plot_ll(flag::Bool=false)
    df = DataFrame(CSV.File(joinpath(fpath,"data/q97.csv")))
    ϕs = [parse(Int,split(ϕnames,"_")[1])//parse(Int,split(ϕnames,"_")[2]) for ϕnames in names(df)]

    fig = figure(figsize=(8,6))
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        q = denominator(ϕ)
        p = numerator(ϕ)
        ϵ = df[:,"$(p)_$(q)"]
        plot(ones(length(ϵ))*p/q,ϵ,"b.",ms=1)
    end
    xlabel(L"ϕ/ϕ_0")
    ylabel(L"ϵ")
    tight_layout()
    display(fig)
    if flag==true 
        savefig("haldane_spectrum.png",transparent=false)
    end
    close(fig)
    return nothing
end

plot_ll(false)

# plot orbital magnetization for fixed particle number
function plot_orbital_magnetization_from_firsttwopoints_fixedN(flag=false)
    fig, ax = subplots(figsize=(4,3))
    qs = [31;61;173]
    colors = ["r","b","g"]
    cnt = 1
    for q in qs
        df = DataFrame(CSV.File(joinpath(fpath,"data/q$(q)_chern0.csv")))
        ϕs = [1//q; 2//q]
        νs = collect(0:size(df,1)) ./(size(df,1)/2)
        moment = zeros(Float64,length(νs))
        for iν in eachindex(νs)
            Eν = zeros(Float64,2)
            if (iν>1)
                for iϕ in eachindex(ϕs)
                    p = numerator(ϕs[iϕ])
                    q = denominator(ϕs[iϕ])
                    ϵ = sort(df[:,"$(p)_$(q)"])
                    Eν[iϕ] = sum(ϵ[1:(iν-1)]) /(length(ϵ)÷2)
                end
                moment[iν] = (Eν[2] - Eν[1])/(ϕs[1]-ϕs[2]) *0.218629 # in units of Bohr magneton
                # moment[iν] = (Eν[2] - Eν[1])/(ϕs[1]-ϕs[2]) *sqrt(3)/(4π) # in units of e/hbar
            end
        end
        ax.plot(νs,moment,".",c=colors[cnt],label="LL q$(q)",ms=2)
        cnt += 1
    end

    ax.plot(hd.νvec,hd.M*1.5862,"k:",label="Niu")
    # ax.plot(hd.νvec,hd.M*2,"r-",label="Niu x 2")
    ax.legend()
    ax.set_xlabel(L"ν")
    ax.set_ylabel(L"M_z (μ_B)")
    # ax.set_title(L"ϕ/ϕ_0=1.5/%$(q)")
    tight_layout()
    display(fig)
    if (flag ==true)
        fname = "M_vs_filling_chern_0_phase_LL_comparison_Niu.pdf"
        savefig(fname,dpi=500)
    end
    close(fig)
end

plot_orbital_magnetization_from_firsttwopoints_fixedN(true)