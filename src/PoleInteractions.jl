export Poles, G, ∇G, ∇G, ∇G , ∇∇G
export G_Evaluation, ∇G_Evaluation, ∇∇G_Evaluation
export VandermontG, Vandermont∇Gn, Vandermontn∇G, Vandermont∇∇Gnn, Vandermontn∇∇Gn, Vandermontnn∇∇G

struct Poles{T}

    Singularity_x::T
    Singularity_y::T
    r_x::Array{T,1}  #this is the position of the pole minus the position of the singularity
    r_y::Array{T,1} 
    n_x::Array{T,1} 
    n_y::Array{T,1} 
    # Distance::Array{Float64,1}

end


function Poles(Singularity, r, n)
    Singularity_x = real(Singularity)
    Singularity_y = imag(Singularity)
    r_x = real(r)
    r_y = imag(r)
    n_x = real(n)
    n_y = imag(n)
    return Poles{typeof(r_x[1])}(Singularity_x, Singularity_y, r_x, r_y, n_x, n_y)
end


export Poles

function Base.getproperty(p::Poles, s::Symbol)
    if s == :Singularity
        return p.Singularity_x + 1im*p.Singularity_y
    elseif s == :r
        return p.r_x + 1im*p.r_y
    elseif s == :n
        return p.n_x + 1im*p.n_y
    else
        return getfield(p, s)
    end
end

function Base.length(p::Poles)
    return length(p.r_x)
end

import Base: +, -

function +(p::Poles, x)
    return  (p.Singularity+x)*( !( p.Singularity + x ≈ 0 ) ) .+ p.r
end

function -(p::Poles, x)
    return  (p.Singularity-x)*( !( p.Singularity - x ≈ 0 ) ) .+ p.r
end


function G(Δr)
    #Δr = x-y
    return log(abs(Δr))/(2π)

end

function ∇G(Δr, n₁)

    real( (1/Δr) * n₁ )/(2π)

end

function ∇∇G(Δr, n₁, n₂ )

    ∂xG = real( (1/Δr)^2 * n₁  )/(2π)
    ∂yG = real( 1im*(1/Δr)^2 * n₁  )/(2π)

    return ∂xG*real(n₂) + ∂yG*imag(n₂) 

end



function G_Evaluation(X::Poles{T}, Y::Poles{T}, μ, ωᵣ) where T

    Gμ⁺ = zeros(length(X))
    Gμ⁻ = zeros(length(X))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 

    for i = 1:length(X)
        for j = 1:length(Y)
            Δr =X.r[i] - Y.r[j] + offset  
            Gμ_j = G(Δr)*μ[j]*ωᵣ[j]
            Gμ⁻[i] += Gμ_j*( Gμ_j <0)
            Gμ⁺[i] += Gμ_j*( Gμ_j >0) 
   
        end
    end

    return Gμ⁺+Gμ⁻

end


function ∇G_Evaluation(X::Poles{T}, Y::Poles{T}, n,μ, ωᵣ) where T

    ∇Gμ⁺ = zeros(length(X.r))
    ∇Gμ⁻ = zeros(length(X.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 

    for i = 1:length(X)
        for j = 1:length(Y)
            Δr = (X.r[i] - Y.r[j]) + offset
            ∇Gμ_j = ∇G(Δr, n[j])*μ[j]*ωᵣ[j]
            ∇Gμ⁻[i] += ∇Gμ_j*( ∇Gμ_j <0)
            ∇Gμ⁺[i] += ∇Gμ_j*( ∇Gμ_j >0) 
   
        end
    end

    return ∇Gμ⁺+∇Gμ⁻

end

function ∇∇G_Evaluation(X::Poles{T}, Y::Poles{T}, n₁, n₂ ,μ,ωᵣ) where T

    ∇∇Gμ⁺ = zeros(length(X.r))
    ∇∇Gμ⁻ = zeros(length(X.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 
    for i = 1:length(X)
        for j = 1:length(Y)
            Δr = (X.r[i] - Y.r[j]) + offset
            ∇∇Gμ_j = ∇∇G(Δr, n₁[j], n₂[i])*μ[j]*ωᵣ[j]
            ∇∇Gμ⁻[i] += ∇∇Gμ_j*( ∇∇Gμ_j <0)
            ∇∇Gμ⁺[i] += ∇∇Gμ_j*( ∇∇Gμ_j >0) 
   
        end
    end

    return ∇∇Gμ⁺+∇∇Gμ⁻

end


function VandermontG(X, Y, ωₗ, ωᵣ)

    VG = zeros(length(X.r), length(Y.r))

    for j = 1:length(Y)
        for i = 1:length(X)

            Δr = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  )   .+ (X.r[i] - Y.r[j])
            VG[i,j] = ωₗ[i]*G(Δr)*ωᵣ[j]

        end
    end

    return VG

end


function Vandermont∇G(X, Y, n, ωₗ, ωᵣ)

    V∇G = zeros(length(X.r), length(Y.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 
    for j = 1:length(Y)
        for i = 1:length(X)

            Δr = (X.r[i] - Y.r[j]) + offset
            V∇G[i,j] = ωₗ[i]*∇G(Δr, n[j])*ωᵣ[j]

        end
    end

    return V∇G

end

function Vandermont∇Gn(X, Y, n, ωₗ, ωᵣ)

    V∇G = zeros(length(X.r), length(Y.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 
    for j = 1:length(Y)
        for i = 1:length(X)

            Δr = (X.r[i] - Y.r[j]) + offset
            V∇G[i,j] = ωₗ[i]*∇G(Δr, n[j])*ωᵣ[j]

        end
    end

    return V∇G

end
function Vandermontn∇G(X, Y, n, ωₗ, ωᵣ)

    V∇G = zeros(length(X.r), length(Y.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 
    for j = 1:length(Y)
        for i = 1:length(X)

            Δr = (X.r[i] - Y.r[j]) + offset
            V∇G[i,j] = ωₗ[i]*∇G(Δr, n[i])*ωᵣ[j]

        end
    end

    return V∇G

end





function Vandermont∇∇G(X, Y, n₁, n₂, ωₗ, ωᵣ)

    V∇∇G = zeros(length(X.r), length(Y.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 

    for j = 1:length(Y)
        for i = 1:length(X)

            Δr = (X.r[i] - Y.r[j]) + offset
            V∇∇G[i,j] = ωₗ[i]*∇∇G(Δr, n₁[j], n₂[i])*ωᵣ[j]

        end
    end

    return V∇∇G

end

function Vandermont∇∇Gnn(X, Y, n₁, n₂, ωₗ, ωᵣ)

    V∇∇G = zeros(length(X.r), length(Y.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 

    for j = 1:length(Y)
        for i = 1:length(X)

            Δr = (X.r[i] - Y.r[j]) + offset
            V∇∇G[i,j] = ωₗ[i]*∇∇G(Δr, n₁[j], n₂[j])*ωᵣ[j]

        end
    end

    return V∇∇G

end


function Vandermontn∇∇Gn(X, Y, n₁, n₂, ωₗ, ωᵣ)

    V∇∇G = zeros(length(X.r), length(Y.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 

    for j = 1:length(Y)
        for i = 1:length(X)

            Δr = (X.r[i] - Y.r[j]) + offset
            V∇∇G[i,j] = ωₗ[i]*∇∇G(Δr, n₁[i], n₂[j])*ωᵣ[j]

        end
    end

    return V∇∇G

end

function Vandermontnn∇∇G(X, Y, n₁, n₂, ωₗ, ωᵣ)

    V∇∇G = zeros(length(X.r), length(Y.r))

    offset = (X.Singularity - Y.Singularity)* ( ! (X.Singularity ≈ Y.Singularity)  ) 

    for j = 1:length(Y)
        for i = 1:length(X)

            Δr = (X.r[i] - Y.r[j]) + offset
            V∇∇G[i,j] = ωₗ[i]*∇∇G(Δr, n₁[i], n₂[i])*ωᵣ[j]

        end
    end

    return V∇∇G

end





