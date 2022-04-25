using PyPlot
using DataFrames
using CSV
include("haldane.jl")
fpath = pwd()

function compute_spectrum()
    df = DataFrame() 

    ϕs = collect(1:96) .// 97
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        q = denominator(ϕ)
        p = numerator(ϕ)
        hd = Haldane()
        fname = joinpath(fpath,"data/p$(p)q$(q).txt")
        constructHaldane(hd;ϕ=ϕ,lk=4,fname=fname)
        df[!,"$(p)_$(q)"] = hd.spectrum[:]
    end

    out_name = joinpath(fpath,"data/q97.csv")
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