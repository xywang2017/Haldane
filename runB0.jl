using PyPlot
using DataFrames
using CSV
include("haldaneB0.jl")
fpath = pwd()

hd = HaldaneB0()
constructHaldaneB0(hd;lk=128)

#
function plot_Mμ(hd)
    fig = figure(figsize=(3,3))
    plot(hd.μvec,hd.M,"k--",label="M")
    plot(hd.μvec,hd.Msr,"m-",label="Msr")
    plot(hd.μvec,hd.Mc,"g-",label="Mc")
    legend()
    xlabel("μ")
    ylabel("M")
    tight_layout()
    display(fig)
    savefig("M_chern_-1_phase.pdf",transparent=true)
    close(fig)
end

plot_Mμ(hd)

#
function plot_Mν(hd)
    fig = figure(figsize=(3,3))
    plot(hd.νvec,hd.M,"k--",label="M")
    plot(hd.νvec,hd.Msr,"m-",label="Msr")
    plot(hd.νvec,hd.Mc,"g-",label="Mc")
    legend()
    xlabel("ν")
    ylabel("M")
    tight_layout()
    display(fig)
    savefig("M_filling_chern_-1_phase.pdf",transparent=true)
    close(fig)
end

plot_Mν(hd)