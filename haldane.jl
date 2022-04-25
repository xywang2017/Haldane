using Parameters
using LinearAlgebra
using DelimitedFiles

@with_kw mutable struct Params 
    t1::Float64 = 1.0
    t2::Float64 = 0.2
    θ::Float64 = π/2
end

mutable struct Haldane 
    params::Params
    p::Int 
    q::Int
    lk::Int
    k1::Vector{Float64}
    k2::Vector{Float64}
    H::Array{ComplexF64,6}  #qx2xqx2xl1xl2, 2 is sublattice label
    spectrum::Array{Float64,3}

    fname::String

    Haldane() = new()
end

function constructHaldane(hd::Haldane;ϕ::Rational=1//10,lk::Int=1,fname="data.txt")
    hd.params = Params() 
    hd.fname = fname
    hd.lk = lk 
    hd.q = denominator(ϕ)
    hd.p = numerator(ϕ)
    hd.k1 = collect(0:(hd.lk-1)) ./ (hd.lk * hd.q)
    hd.k2 = collect(0:(hd.lk-1)) ./ (hd.lk * hd.q)


    constructH(hd)
    computeSpectrum(hd)
    # writeout(hd)

    return nothing
end

function constructH(hd::Haldane)
    ϕ = hd.p/hd.q
    θ = hd.params.θ
    hd.H = zeros(ComplexF64,hd.q,2,hd.q,2,hd.lk,hd.lk)
    for i2 in 1:hd.lk, i1 in 1:hd.lk 
        k1 = hd.k1[i1]
        
        for r1 in 1:hd.q 
            k2 = hd.k2[i2] + (r1-1)/hd.q
            # BA
            hd.H[r1,2,r1,1,i1,i2] = hd.params.t1 * (exp(-1im*(π/6)*ϕ) + exp(1im*(π/6)*ϕ+1im*2π*k2)) * exp(-1im *2π*k2/3)
            hd.H[r1,2,mod(r1+hd.p-1,hd.q)+1,1,i1,i2] = hd.params.t1 * (exp(1im*2π*k1-1im*π*ϕ))* exp(-1im *2π*k2/3)
            
            # AA
            hd.H[r1,1,r1,1,i1,i2] = hd.params.t2 * (exp(1im*θ-1im*(0)*ϕ+1im*2π*k2)) 
            hd.H[r1,1,mod(r1+hd.p-1,hd.q)+1,1,i1,i2] = hd.params.t2 * (exp(-1im*θ - 1im*(π/2)*ϕ + 1im*2π*k1) + 
                                            exp(1im*θ - 1im*(3π/2)*ϕ + 1im*2π*(k1-k2)) )
            # BB
            hd.H[r1,2,r1,2,i1,i2] = hd.params.t2 * (exp(-1im*θ+1im*(2π/3)*ϕ+1im*2π*k2)) 
            hd.H[r1,2,mod(r1+hd.p-1,hd.q)+1,2,i1,i2] = hd.params.t2 * (exp(1im*θ - 1im*(π/6)*ϕ + 1im*2π*k1) + 
                                            exp(-1im*θ - 1im*(11π/6)*ϕ + 1im*2π*(k1-k2)) )
        end
    end
    H = reshape(hd.H,2hd.q,2hd.q,hd.lk,hd.lk)
    for i2 in 1:hd.lk, i1 in 1:hd.lk 
        H[:,:,i1,i2] = H[:,:,i1,i2] + H[:,:,i1,i2]'
    end
    return nothing
end

function computeSpectrum(hd::Haldane)
    hd.spectrum = zeros(Float64,2hd.q,hd.lk,hd.lk)
    H = reshape(hd.H,2hd.q,2hd.q,hd.lk,hd.lk)
    for i2 in 1:hd.lk, i1 in 1:hd.lk 
        hd.spectrum[:,i1,i2] = eigvals(Hermitian(H[:,:,i1,i2]))
    end
    return nothing
end

function writeout(hd::Haldane)
    writedlm(hd.fname,hd.spectrum[:])
    return nothing
end