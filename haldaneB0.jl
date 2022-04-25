using Parameters
using LinearAlgebra
using DelimitedFiles

@with_kw mutable struct Params 
    t1::Float64 = 1.0
    t2::Float64 = 0.2
    θ::Float64 = π/2

    σ0::Matrix{ComplexF64} = ComplexF64[1 0;0 1]
    σx::Matrix{ComplexF64} = ComplexF64[0 1;1 0]
    σy::Matrix{ComplexF64} = ComplexF64[0 -1im;1im 0]
    σz::Matrix{ComplexF64} = ComplexF64[1 0;0 -1]
end

mutable struct HaldaneB0 
    params::Params
    lk::Int
    k1::Vector{Float64}
    k2::Vector{Float64}
    
    Uk::Array{ComplexF64,4}  # 2 x 2 x lk x lk
    ϵn::Array{Float64,3} # 2 x lk x lk

    μvec::Vector{Float64} # vector of chemical potentials to calculate magnetic moment 
    νvec::Vector{Float64} # vector of filling fractions corresponding to the chemical potential 
    Msr::Vector{Float64} #
    Mc::Vector{Float64} #
    M::Vector{Float64} #

    HaldaneB0() = new()
end

function constructHaldaneB0(hd::HaldaneB0;lk::Int=64)
    hd.lk = lk 
    hd.k1 = collect(0:(hd.lk-1)) ./ hd.lk 
    hd.k2 = collect(0:(hd.lk-1)) ./ hd.lk 
    hd.params = Params(θ=π/2)

    computeSpectrum(hd)
    ϵ_max = maximum(hd.ϵn)
    ϵ_min = minimum(hd.ϵn)
    hd.μvec = ϵ_min .+ (ϵ_max -ϵ_min)/100 * collect(0.5:1:99.5) 
    hd.νvec = [sum((sign.(hd.μvec[iμ] .- hd.ϵn[:]) .+ 1)./2)/hd.lk^2 for iμ in eachindex(hd.μvec)]
    computeOrbitalMagneticMoment(hd)

    return nothing
end

function computeSpectrum(hd::HaldaneB0)
    hd.Uk = zeros(ComplexF64,2,2,hd.lk,hd.lk)
    hd.ϵn = zeros(Float64,2,hd.lk,hd.lk)

    for i2 in eachindex(hd.k2), i1 in eachindex(hd.k1)
        F = eigen(Hermitian(H(hd.k1[i1],hd.k2[i2],hd.params)))
        hd.Uk[:,:,i1,i2] = F.vectors 
        hd.ϵn[:,i1,i2] = F.values
    end
    return nothing
end

function computeOrbitalMagneticMoment(hd::HaldaneB0)
    hd.Msr = zeros(Float64,length(hd.μvec))
    hd.Mc = zeros(Float64,length(hd.μvec))
    Msr, Mc = 0.0+0.0im, 0.0 + 0.0im
    for iμ in eachindex(hd.μvec)
        for i2 in eachindex(hd.k2), i1 in eachindex(hd.k1)
            k1,k2 = hd.k1[i1], hd.k2[i2]
            for idx_band in 1:2
                if (sign(hd.μvec[iμ] - hd.ϵn[idx_band,i1,i2])>0)
                    idx_band_pair = 3 - idx_band
                    Msr = - ( hd.Uk[:,idx_band,i1,i2]' * ∂xH(k1,k2,hd.params) * hd.Uk[:,idx_band_pair,i1,i2]) * 
                          ( hd.Uk[:,idx_band_pair,i1,i2]' * ∂yH(k1,k2,hd.params) * hd.Uk[:,idx_band,i1,i2]) /
                          (hd.ϵn[idx_band,i1,i2] - hd.ϵn[idx_band_pair,i1,i2])
                    Mc = -2 * ( hd.Uk[:,idx_band,i1,i2]' * ∂xH(k1,k2,hd.params) * hd.Uk[:,idx_band_pair,i1,i2]) * 
                          ( hd.Uk[:,idx_band_pair,i1,i2]' * ∂yH(k1,k2,hd.params) * hd.Uk[:,idx_band,i1,i2]) /
                          (hd.ϵn[idx_band,i1,i2] - hd.ϵn[idx_band_pair,i1,i2])^2  * (hd.μvec[iμ]-hd.ϵn[idx_band,i1,i2])
                    
                    hd.Msr[iμ] += imag(Msr)
                    hd.Mc[iμ] += imag(Mc)
                end
            end
        end
    end
    hd.Msr ./= hd.lk^2
    hd.Mc ./= hd.lk^2
    hd.M = hd.Msr + hd.Mc
    return nothing
end

function H(k1::Float64,k2::Float64,params::Params)
    ϵ = params.t1 * (1 + cos(2π*k1) + cos(2π*k2)) * params.σx + 
        params.t1 * (sin(2π*k1) + sin(2π*k2)) * params.σy + 
        2params.t2 * cos(params.θ) * (cos(2π*k1) + cos(2π*k2) + cos(2π*(k2-k1))) * params.σ0 + 
        2params.t2 * sin(params.θ) * (sin(2π*k1) - sin(2π*k2) + sin(2π*(k2-k1))) * params.σz
    return ϵ
end

function ∂1H(k1::Float64,k2::Float64,params::Params)
    d1H = 2π * params.t1 * (-sin(2π*k1)*params.σx + cos(2π*k1)*params.σy) +
          4π * params.t2 * cos(params.θ) * (-sin(2π*k1) + sin(2π*(k2-k1))) * params.σ0 + 
          4π * params.t2 * sin(params.θ) * (cos(2π*k1) - cos(2π*(k2-k1))) * params.σz
    return d1H
end

function ∂2H(k1::Float64,k2::Float64,params::Params)
    d2H = 2π * params.t1 * (-sin(2π*k2)*params.σx + cos(2π*k2)*params.σy) +
          4π * params.t2 * cos(params.θ) * (-sin(2π*k2) - sin(2π*(k2-k1))) * params.σ0 + 
          4π * params.t2 * sin(params.θ) * (-cos(2π*k2) + cos(2π*(k2-k1))) * params.σz
    return d2H
end

function ∂xH(k1::Float64,k2::Float64,params::Params)
    return sqrt(3)/(4π) * ∂1H(k1,k2,params)
end

function ∂yH(k1::Float64,k2::Float64,params::Params)
    return 1/(2π) * (0.5*∂1H(k1,k2,params)+∂2H(k1,k2,params))
end