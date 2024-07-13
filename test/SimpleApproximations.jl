using OpenArcDIM
using CairoMakie
using Revise
using LinearAlgebra
using BenchmarkTools

#Here we approximate the functions

#f= sqrt(x+a)
#f = 1/sqrt(x+a)

a = 1.

f1(x) = sqrt.(real.(x+a))
f2(x) = 1/sqrt.(real.(x+a) )


#Cluster Poles
sing_point1 = -a + 0im
θ₁ = π 
d̂₁ =  cos(θ₁) + 1im*sin(θ₁)
n̂₁ = 1im*d̂₁

numpoles = 4000
σ = 2π*sqrt(2)
T = sqrt.(1:numpoles) .- sqrt(numpoles)
Dist = 4*exp.(σ*T) 
Dist = Dist[Dist .> 10^(-14)]
numpoles = length(Dist)
M = (Int ∘ ceil)(3*sqrt(numpoles))

Y = Poles(sing_point1, d̂₁* Dist  , n̂₁ * ones(numpoles)) 


#FarAway Poles

r = 1.5*exp.(2π*1im*(0:M-1)/M)
Ỹ = Poles(sing_point1, r, 1im*exp.(2π*1im*(0:M-1)/M) )



#Samples 

number_sample_points = 1.5*numpoles .|> ceil .|> Int
cutoff = 15
step = cutoff /number_sample_points

p = 8  # This type of concentration is better because it puts points evenly at the end of the segment
v(s,p) = (1/p - 1/2)*(1-2*s)+ (1/p)*(2*s -1) + 1/2
w(s,p) = 2*( v(s,p)^p/(v(s,p)^p + v(1. -s,p)^p) )


# DistSamples = 1.0*10 .^( cutoff*( sqrt.( (0:number_sample_points-1)/number_sample_points) .-1 ) )
DistSamples = w.(collect(1:number_sample_points)*0.5/number_sample_points, p)
DistSamples = DistSamples[DistSamples .> 10^(-14)]
number_sample_points = length(DistSamples)
Samples = Poles(sing_point1 , DistSamples  .+ 0im,  -1im*ones(number_sample_points ) )

# Vandermont Matrices Acumulation

VS = VandermontG(Samples, Y, ones(length(Samples)), ones(length(Y)))
VD = Vandermont∇G(Samples, Y, Y.n*1im , ones(length(Samples)),abs.(Y.r))
VN = Vandermont∇∇G(Samples, Y, Y.n, Samples.n, ones(length(Samples)),abs.(Y.r).^2)

#MFS Smooth Matrix

VMFS = VandermontG(Samples, Ỹ, ones(length(Samples)), ones(length(Y)))

V1 = [VS VMFS]
V2 = [VD VMFS]
V3 = [VN VMFS]

#RHS
b = f1(Samples)

μ1 = V1\b
μ2 = V2\b
μ3 = V3\b

norm(V1*μ1 - b )
norm(V2*μ2 - b )
norm(V3*μ3 - b )





#Test Interpolation error
TestPointsNumber = 500
Test_points = Poles( sing_point1, w.(collect(1:TestPointsNumber)*0.5/TestPointsNumber, p-1), -1im*ones(TestPointsNumber) ) 
Test_eval_1 = G_Evaluation(Test_points, Y, μ1[1:numpoles], ones( length(Y)))  + G_Evaluation(Test_points, Ỹ, μ1[numpoles+1:end], ones( length(Ỹ))) 
Test_eval_2 = ∇G_Evaluation(Test_points, Y, 1im*Y.n, μ2[1:numpoles], abs.( Y.r))  + G_Evaluation(Test_points, Ỹ, μ2[numpoles+1:end], ones( length(Ỹ)))
Test_eval_3 = ∇∇G_Evaluation(Test_points, Y, Y.n, Test_points.n, μ3[1:numpoles], abs.( Y.r).^2)  + G_Evaluation(Test_points, Ỹ, μ3[numpoles+1:end], ones( length(Ỹ)))

Test_eval_exact = f1(Test_points)


norm(Test_eval_1 - Test_eval_exact, Inf)
norm(Test_eval_2 - Test_eval_exact, Inf)
norm(Test_eval_3 - Test_eval_exact, Inf)


fig = Figure()
ax = Axis(fig[1, 1],yscale = log10, xscale = log10)
lines!(ax, real.(Test_points.r)  ,((Test_eval_1 - Test_eval_exact) .+ 10^(-16)) .|> abs)
lines!(ax, real.(Test_points.r)  ,((Test_eval_2 - Test_eval_exact) .+ 10^(-16)) .|> abs)
lines!(ax, real.(Test_points.r)  ,((Test_eval_3 - Test_eval_exact) .+ 10^(-16)) .|> abs)
lines!(ax, DistSamples[1]*ones(13), 10.0 .^( -16:0.5:-10) , linestyle = :dash, color = :black)
fig

#Test Extrapolation error, towards 0

number_samples_extra = 32
Test_extrapolation = Poles( sing_point1, 10.0 .^(-number_samples_extra :0.5:-14) , -1im*ones(2*(number_samples_extra - 14) + 1) )


Test_eval_1 = G_Evaluation(Test_extrapolation, Y, μ1[1:numpoles], ones( length(Y)))  + G_Evaluation(Test_extrapolation, Ỹ, μ1[numpoles+1:end], ones( length(Ỹ))) 
Test_eval_2 = ∇G_Evaluation(Test_extrapolation, Y, 1im*Y.n, μ2[1:numpoles], abs.( Y.r))  + G_Evaluation(Test_extrapolation, Ỹ, μ2[numpoles+1:end], ones( length(Ỹ)))
Test_eval_3 = ∇∇G_Evaluation(Test_extrapolation, Y, Y.n, Test_extrapolation.n, μ3[1:numpoles], abs.( Y.r).^2)  + G_Evaluation(Test_extrapolation, Ỹ, μ3[numpoles+1:end], ones( length(Ỹ)))
Test_eval_exact_extra = f1(Test_extrapolation)

norm(Test_eval_1 - Test_eval_exact_extra, Inf)
norm(Test_eval_2 - Test_eval_exact_extra, Inf)
norm(Test_eval_3 - Test_eval_exact_extra, Inf)

fig = Figure()
ax = Axis(fig[1, 1],yscale = log10, xscale = log10)
lines!(ax, real.(Test_extrapolation.r)  ,((Test_eval_1 - Test_eval_exact_extra) .+ 10^(-16)) .|> abs)
lines!(ax, real.(Test_extrapolation.r)  ,((Test_eval_2 - Test_eval_exact_extra) .+ 10^(-16)) .|> abs)
lines!(ax, real.(Test_extrapolation.r)  ,((Test_eval_3 - Test_eval_exact_extra) .+ 10^(-16)) .|> abs)
fig