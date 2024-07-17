using OpenArcDIM
using CairoMakie
using Revise
using LinearAlgebra
using BenchmarkTools
using KrylovKit


α = 0.0
β = 1.0
φ(t) = exp(cos(3*t)) 


#Sampling rate

p = 3  # This type of concentration is better because it puts points evenly at the end of the segment
v(s,p) = (1/p - 1/2)*(1-2*s)^3 + (1/p)*(2*s -1) + 1/2
w(s,p) = 2*( v(s,p)^p/(v(s,p)^p + v(1. -s,p)^p) )



#Geometry



γ(t) = ( -1im*( cos(t) + 1im*sin(t))  + 1im ) #- ( -8.190169853792138e-18 + 0.01542811005090392im)
γ′(t) = -1im*( sin(t) + 1im*cos(t))
#Approximated Geometry


γ_0 = γ(0.) 
# γ′_0 = sum(γ′.(T))/length(T)
γ′_0 = γ′(0.)#sum( (γ.(T) .- γ_0).*T )/sum(T.*T)
γ′_0 = γ′_0/norm(γ′_0)
γ̃(t) = real(conj(γ(t) - γ_0)*γ′_0)*γ′_0 

numpoles = 500
σ = 2π*sqrt(2)
T = sqrt.(1:numpoles) .- sqrt(numpoles)
Dist = 0.4*exp.(σ*T) 
Dist = Dist[Dist .> 10^(-15)]
numpoles = length(Dist)
M =(Int ∘ ceil)(2*sqrt(numpoles))
M = M +  M%2


Ycluster = Poles(γ_0, -γ′_0 * Dist  , 1im*γ′_0 * ones(numpoles)) 
Yfar = Poles(γ_0, 0.5*exp.(2π*1im*(1/2:M-1+1/2)/M), 1im*exp.(2π*1im*(0:M-1)/M) )
Yfar_conjugated = Poles(γ_0, conj.(Yfar.r), conj.(Yfar.r) )


#Samples

number_sample_points =  4*numpoles |> ceil |> Int
T_samples = (w.( collect(0:1:number_sample_points-1)/(number_sample_points-1), p))*0.1
r= γ.( T_samples)
N = γ′.(T_samples )
N = N./abs.(N)*1im

Samples = Poles(-γ_0, r, N )


fig = Figure(size = (600, 600))
ax = Axis(fig[1, 1], aspect = 1)
scatter!(ax, real.(Yfar.r), imag.(Yfar.r), color = :green, markersize = 10, marker =:x)
# scatter!(ax, real.(Yfar_conjugated.r), imag.(Yfar_conjugated.r), color = :purple)
scatter!(ax, real.(Ycluster.r), imag.(Ycluster.r), color = :blue, markersize = 10, marker =:x)
scatter!(ax, real.(Samples.Singularity .+ Samples.r), imag.(Samples.Singularity .+Samples.r), color = :black)
# scatter!(ax, real.(Samples_aprox.Singularity .+ Samples_aprox.r), imag.(Samples_aprox.Singularity .+Samples_aprox.r), color = :orange)
fig