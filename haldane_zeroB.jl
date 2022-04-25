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
    
    Msr::Array{Float64,3}
    Mc::Array{Float64,3}

    HaldaneB0() = new()
end

function H(k1::Float64,k2::Float64,params::Params)
    return 0
end

function ∂1H(k1::Float64,k2::Float64,params::Params)
    d1H = 2π * params.t1 * (-sin(2π*k1)*params.σx)
    return 0
end

function ∂2H(k1::Float64,k2::Float64,params::Params)
    return 0
end

function ∂xH(k1::Floata64,k2::Float64,params::Params)
    return sqrt(3)/(4π) * ∂1H(k1,k2,params)
end

function ∂yH(k1::Floata64,k2::Float64,params::Params)
    return 1/(2π) * (0.5*∂1H(k1,k2,params)+∂2H(k1,k2,params))
end