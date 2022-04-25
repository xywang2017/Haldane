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
    plot(hd.μvec,hd.M*1.5862,"k--",label="M")
    plot(hd.μvec,hd.Msr*1.5862,"m-",label="Msr")
    plot(hd.μvec,hd.Mc*1.5862,"g-",label="Mc")
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
    # from e/hbar * a^2*eV to α * μB, α = 1.5862
    fig = figure(figsize=(4,3))
    plot(hd.νvec,hd.M*1.5862,"b-",label="M")
    # plot(hd.νvec,hd.Msr*1.5862,"m-",label="Msr")
    # plot(hd.νvec,hd.Mc*1.5862,"g-",label="Mc")
    legend()
    xlabel("ν")
    ylabel(L"M (μ_B)")
    tight_layout()
    display(fig)
    savefig("M_filling_chern_1_phase.pdf",transparent=true)
    close(fig)
end

plot_Mν(hd)