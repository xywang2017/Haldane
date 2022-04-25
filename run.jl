using PyPlot
using DataFrame
using CSV
include("haldane.jl")
fpath = pwd()

function compute_spectrum()
    ϕs = collect(1:96) .// 97
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        q = denominator(ϕ)
        p = numerator(ϕ)
        hd = Haldane()
        fname = joinpath(fpath,"data/p$(p)q$(q).txt")
        constructHaldane(hd;ϕ=ϕ,lk=4,fname=fname)
    end

    # data processing 
    df = DataFrame() 
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        q = denominator(ϕ)
        p = numerator(ϕ)
        df[!,"$p_$q"] = readdlm(joinpath(fpath,"data/p$(p)q$(q).txt"))
    end
    out_name = joinpath(fpath,"data/q$q.csv")
    CSV.write(out_name,df,delim=",")
end

##
function plot_ll(flag::Bool=false)
    fig = figure(figsize=(8,6))
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        q = denominator(ϕ)
        p = numerator(ϕ)
        fname = joinpath(fpath,"data/p$(p)q$(q).txt")
        ϵ = readdlm(fname)
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