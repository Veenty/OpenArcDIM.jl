using OpenArcDIM
using CairoMakie
using Revise
using LinearAlgebra
using BenchmarkTools


function ArcParametrization(θ, r, center)
    return center + r*exp(1im*θ)
end

sing_point1 = 1. + 0im
θ₁ = π/2*3 + π*0.0 
d̂₁ =  cos(θ₁) + 1im*sin(θ₁)
n̂₁ = 1im*d̂₁


numpoles = 500
σ = 2π*sqrt(2)
T = sqrt.(1:numpoles) .- sqrt(numpoles)
Dist = 0.5*exp.(σ*T) 
Dist = Dist[Dist .> 10^(-15)]
numpoles = length(Dist)
M = (Int ∘ ceil)(2*sqrt(numpoles))

Yext = Poles(sing_point1, d̂₁* Dist  , n̂₁ * ones(numpoles)) 
Yint = Poles(sing_point1,   (-conj.(Yext.r) - abs.(Yext.r).^2 ) ./abs.(1. .+ Yext.r).^2  , n̂₁ * ones(numpoles))



Sr_x = -2*sin.(-Dist/2).^2
Sr_y = sin.(-Dist) 
Yboundary  = Poles(sing_point1 , Sr_x + 1im*Sr_y, -1im*(sing_point1 .+ (Sr_x + 1im*Sr_y) )  ) 



#FarAway Poles

r = 0.5*exp.(2π*1im*(0:M-1)/M)
Yfar = Poles(sing_point1, r, 1im*exp.(2π*1im*(0:M-1)/M) )



Θ = 2π*(-100:100)/300
circle_points = ArcParametrization.(Θ, 1.0, 0.0)





number_sample_points = 1.5*numpoles .|> ceil .|> Int
cutoff = 15
step = cutoff /number_sample_points

p = 8  # This type of concentration is better because it puts points evenly at the end of the segment
v(s,p) = (1/p - 1/2)*(1-2*s)^3 + (1/p)*(2*s -1) + 1/2
w(s,p) = 2*( v(s,p)^p/(v(s,p)^p + v(1. -s,p)^p) )

DistSamples = w.(collect(1:number_sample_points)*0.5/number_sample_points, p)
DistSamples = 2π*(0.06)*DistSamples[DistSamples .> 10^(-14)]
number_sample_points = length(DistSamples)

r_x = -2*sin.(DistSamples/2).^2
r_y = sin.(DistSamples) 
Samples = Poles(sing_point1 , r_x + 1im*r_y, -1im*(r_x + 1im*r_y)  ) 


f1(θ) = sqrt(θ)*cos(θ)

f2(θ) = cos(θ)


Vint = VandermontG(Samples, Y, ones(length(Samples)), ones(length(Y)))
Vext = VandermontG(Samples, Ȳ, ones(length(Samples)), ones(length(Ȳ)))
Vfar = VandermontG(Samples, Ỹ, ones(length(Samples)), ones(length(Ỹ)))

Vboundary = Vandermont∇Gn(Samples, Yboundary , Yboundary.Singularity .+ Yboundary.r , ones(length(Samples)), abs.(Yboundary.r)) - Vandermontn∇G(Samples, Yboundary , Samples.Singularity .+ Samples.r , ones(length(Samples)), abs.(Yboundary.r)) 

U, S, V = svd(Vboundary)

heatmap(Vboundary)



V = [Vint+Vext Vfar]

b = f1.(DistSamples) 

plot(DistSamples, b)

μ = V\b
norm(V*μ - b)







fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], aspect = 1)
# fig[1, 1] = ax

scatter!(ax,Yint.Singularity_x .+ Yint.r_x, Yint.Singularity_y .+ Yint.r_y )
scatter!(ax,Yext.Singularity_x .+ Yext.r_x, Yext.Singularity_y .+ Yext.r_y )
scatter!(ax,Yfar.Singularity_x .+ Yfar.r_x, Yfar.Singularity_y .+ Yfar.r_y )
scatter!(ax,Samples.Singularity_x .+ Samples.r_x, Samples.Singularity_y .+ Samples.r_y )
# scatter!(ax,Yboundary.Singularity_x .+ Yboundary.r_x, Yboundary.Singularity_y .+ Yboundary.r_y , marker =:diamond)
lines!(ax, real.(circle_points ), imag.(circle_points ), color = :black)
fig

save("sources_circle.png", fig)



(1 .+ Yext.r[end])/abs(1 .+ Yext.r[end])
(1 .+ Yint.r[end])/abs(1 .+ Yint.r[end])