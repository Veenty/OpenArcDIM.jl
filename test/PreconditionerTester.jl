using OpenArcDIM
using CairoMakie
using Revise
using LinearAlgebra
using BenchmarkTools
using KrylovKit


α = 1.0
β = 1.0
φ(t) = exp(cos(3*t)) 

#Geometry

γ(t) = ( -1im*( cos(t) + 1im*sin(t))  + 1im ) #- ( -8.190169853792138e-18 + 0.01542811005090392im)
γ′(t) = -1im*( sin(t) + 1im*cos(t))
#Approximated Geometry
T = collect(-0.3:0.001:0.3)#2π*(-100:100)/3000

γ_0 = sum(γ.(T))/length(T)
# γ′_0 = sum(γ′.(T))/length(T)
γ′_0 = sum( (γ.(T) .- γ_0).*T )/sum(T.*T)
γ′_0 = γ′_0/norm(γ′_0)
γ̃(t) = real(conj(γ(t) - γ_0)*γ′_0)*γ′_0 


#FarAway Poles
M_list = collect(2:2:40)
error_direct = zeros(length(M_list))
error_special = zeros(length(M_list))
error_special_aprox_L = zeros(length(M_list))
error_special_aprox_R = zeros(length(M_list))

for i in eachindex(M_list)
# i = 8
M = M_list[i]
r = exp(2π*1im*0.5/M)*0.5*exp.(2π*1im*(0:M-1)/M)
Yfar = Poles(0.0, r, r )

Yfar_conjugated = Poles(0.0, conj.(Yfar.r), conj.(Yfar.r) )



#Samples 

number_sample_points =  4*M |> ceil |> Int

p = 3  # This type of concentration is better because it puts points evenly at the end of the segment
v(s,p) = (1/p - 1/2)*(1-2*s)^3 + (1/p)*(2*s -1) + 1/2
w(s,p) = 2*( v(s,p)^p/(v(s,p)^p + v(1. -s,p)^p) )
T_samples = (w.( collect(0:1:number_sample_points-1)/(number_sample_points-1), p).-1)*0.3 
r= γ.( T_samples)
N = γ′.(T_samples )
N = N./abs.(N)*1im

Samples = Poles(-γ_0, r, N )

#Aproximated Samples points

r̃ = γ̃.( T_samples)
Ñ = γ′_0*1im
Samples_aprox = Poles(0.0, r̃, Ñ*ones(length(r̃)) )


#Right Hand Side

b̃ = φ.(T_samples)
b = [α*b̃ ; β*b̃ ]


#Direct Aproach

S_1 = VandermontG(Samples, Yfar, ones(length(Samples)), ones(length(Yfar)))
S_direct = [S_1 S_1]
∂S_1 = Vandermontn∇G(Samples, Yfar, Samples.n ,ones(length(Samples)), ones(length(Yfar)))
∂S_direct = [∂S_1 ∂S_1]
V_direct = [S_direct ; ∂S_direct ]


coef = V_direct\b

error_direct[i] = norm(V_direct*coef - b)
# norm(S_direct*coef - α*b̃)
# norm(∂S_direct*coef - β*b̃)

#Using Special Green function

S_D = VandermontG(Samples, Yfar_conjugated, ones(length(Samples)), ones(length(Yfar))) - VandermontG(Samples, Yfar, ones(length(Samples)), ones(length(Yfar)))
S_N = VandermontG(Samples, Yfar_conjugated, ones(length(Samples)), ones(length(Yfar))) + VandermontG(Samples, Yfar, ones(length(Samples)), ones(length(Yfar)))

∂S_D = Vandermontn∇G(Samples, Yfar_conjugated, Samples.n ,ones(length(Samples)), ones(length(Yfar))) - Vandermontn∇G(Samples, Yfar, Samples.n ,ones(length(Samples)), ones(length(Yfar)))
∂S_N = Vandermontn∇G(Samples, Yfar_conjugated, Samples.n ,ones(length(Samples)), ones(length(Yfar))) + Vandermontn∇G(Samples, Yfar, Samples.n ,ones(length(Samples)), ones(length(Yfar)))

V_special = [S_N S_D ; ∂S_N ∂S_D]

coef_special = V_special\b

error_special[i] = norm(V_special*coef_special - b)
# norm( [S_N S_D]coef_special - α*b̃)
# norm( [∂S_N ∂S_D]coef_special - β*b̃)

# S = [S_D S_N]
# norm(S*(S\(α*b̃)) - α*b̃)

#Using Special Green function with Aproximated Samples

#S̃_D = VandermontG(Samples_aprox, Yfar_conjugated, ones(length(Samples_aprox)), ones(length(Yfar))) - VandermontG(Samples_aprox, Yfar, ones(length(Samples_aprox)), ones(length(Yfar)))
S̃_N = VandermontG(Samples_aprox, Yfar_conjugated, ones(length(Samples_aprox)), ones(length(Yfar))) + VandermontG(Samples_aprox, Yfar, ones(length(Samples_aprox)), ones(length(Yfar)))

∂S̃_D = Vandermontn∇G(Samples_aprox, Yfar_conjugated, Samples_aprox.n ,ones(length(Samples_aprox)), ones(length(Yfar))) - Vandermontn∇G(Samples_aprox, Yfar, Samples_aprox.n ,ones(length(Samples_aprox)), ones(length(Yfar)))
#∂S̃_N = Vandermontn∇G(Samples_aprox, Yfar_conjugated, Samples_aprox.n ,ones(length(Samples_aprox)), ones(length(Yfar))) + Vandermontn∇G(Samples_aprox, Yfar, Samples_aprox.n ,ones(length(Samples_aprox)), ones(length(Yfar)))

# tol = 1e-13

# U_N, Σ_N, Vt_N = svd(S̃_N)
# Σ⁻¹_N = zeros(length(Σ_N))
# Σ⁻¹_N[Σ_N .> tol] =  1 ./ Σ_N[Σ_N .> tol]

# U_D, Σ_D, Vt_D = svd(∂S̃_D)
# Σ⁻¹_D = zeros(length(Σ_D))
# Σ⁻¹_D[Σ_D .> tol] =  1 ./ Σ_D[Σ_D .> tol]



Δ_N = S_N- S̃_N
∂Δ_D = ∂S_D - ∂S̃_D

# plot(T_samples, ∂Δ_D[:, 3])

# [S_N S_D ; ∂S_N ∂S_D]
C = [Δ_N  S_D ; ∂S_N ∂Δ_D]

# S̃_N

function precon_right(x) 

    N = Int(length(x) /2)
    λ₁ = S̃_N\(x[1:N]) 
    λ₂ = ∂S̃_D\(x[N+1:end])
    # λ₁ = Vt_N*(Σ⁻¹_N.*(U_N'*x[1:N]))
    # λ₂ = Vt_D*(Σ⁻¹_D.*(U_D'*x[N+1:end]))

    return x + C*[λ₁ ; λ₂]

end

function precon_left(x)

    Cx = C*x
    N = Int(length(Cx) /2)
    λ₁ = S̃_N\Cx[1:N]
    λ₂ = ∂S̃_D\(Cx[N+1:end])

    return x + [λ₁ ; λ₂]



end


b_precon_left = [ S̃_N\(α*b̃) ; ∂S̃_D\(β*b̃) ]
A⁻¹V = [I+S̃_N\(Δ_N)  S̃_N\S_D ; ∂S̃_D\∂S_N I+∂S̃_D\∂Δ_D]


x_left, info_left = linsolve(precon_left, b_precon_left, tol = 1e-9)
error_special_aprox_L[i] = norm(V_special*x_left - b)
# norm( precon_left( b_precon_left) - A⁻¹V*b_precon_left, Inf)

x_right, info_right = linsolve(precon_right, b, tol = 1e-13)

# plot( T_samples, (precon_right(x_right) - b .+ 10^(-16))[1:length(x_right)÷2] .|> abs .|> log10)
# plot( T_samples, (precon_right(x_right) - b .+ 10^(-16))[length(x_right)÷2+1:end] .|> abs .|> log10)

# plot( T_samples, α*x_right[1:length(x_right)÷2])


println("M = ", M)
println(info_right)
error_special_aprox_R[i] = norm( precon_right(x_right) - b  )

end

# vals_right, vecs, info_eig_right = eigsolve( precon_right, b, tol = 1e-15, 250; krylovdim = 250);
# vals_left, vecs, info_eig_left = eigsolve( precon_left, b_precon_left, tol = 1e-12, 100; krylovdim = 100);


fig = Figure(size = (600, 600))
ax = Axis(fig[1, 1], aspect = 1, yscale = log10)
plot!(ax, M_list, error_direct, color = :red, label = "Direct")
plot!(ax, M_list, error_special, color = :blue, label = "Special")
plot!(ax, M_list, error_special_aprox_L, color = :green, label = "Special Aprox Left")
plot!(ax, M_list, error_special_aprox_R, color = :purple, label = "Special Aprox Right")
Legend(fig[1,2], ax)
fig



#

# fig = Figure(size = (600, 600))
# ax = Axis(fig[1, 1], aspect = 1)
# scatter!(ax, real.(vals_left), imag.(vals_left), color = :red, markersize = 20)
# scatter!(ax, real.(vals_right), imag.(vals_right), color = :blue)
# fig




_, Σ, _ = svd(  S̃_N  )
sum(Σ .> 10^(-10))
plot(Σ)


fig = Figure(size = (600, 600))
ax = Axis(fig[1, 1], aspect = 1)
scatter!(ax, real.(Yfar.r), imag.(Yfar.r), color = :green, markersize = 40, marker =:x)
scatter!(ax, real.(Yfar_conjugated.r), imag.(Yfar_conjugated.r), color = :purple)
scatter!(ax, real.(Samples.Singularity .+ Samples.r), imag.(Samples.Singularity .+Samples.r), color = :black)
scatter!(ax, real.(Samples_aprox.Singularity .+ Samples_aprox.r), imag.(Samples_aprox.Singularity .+Samples_aprox.r), color = :orange)
fig


