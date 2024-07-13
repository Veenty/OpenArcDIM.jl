using Polynomials

function Quadrature_weights(n::Int)
    # Compute the weights for the quadrature rule using Cherbyshev points
    #Integral is from 0 to 1.
    #We want to be extremely accurate for the points near 0
    #We don't evaluate the function at 0.
    #Chebyshev points 

    θ = π*collect((2*(1:n) .- 1)/(2*n))
    x = (sin.(θ/2)).^2
    ŵ = zeros(n)
    ŵ[1] = 1.
    for i ∈ 3:n
        ŵ[i] = (cos(π*(i-1)) +1)/(2 - 2*(i-1)^2)
    end

    w = (vander(ChebyshevT, 1 .- 2*x , n-1))' \ ŵ

    return x, w

end


function quadIntegral(f,a,b;tol = 1e-8, levels = 10)
    x,w = Quadrature_weights(  -(Int∘floor∘log10)(tol)   )
    return RecursiveIntegral(f, a, b, x, w; tol = tol, levels = levels)
end



function RecursiveIntegral(f, a, b, x, w; tol = 1e-8, levels = 3)
    scale = (b-a)
    mid = (a+b)/2
    Iₜ = sum(scale*w .* f.(a .+scale*x))
    I₁ = sum((scale*0.5)*w .* f.(a .+ (scale*0.5)*x))
    I₂ = sum((scale*0.5)*w .* f.(mid .+ (scale*0.5)*x))

    # println("Iₜ: ", Iₜ)
    
    # println(I₁ + I₂)
    if abs(I₁ + I₂ - Iₜ) < tol || levels == 0
        return I₁ + I₂, abs(I₁ + I₂ - Iₜ) 
    else
        return RecursiveIntegral(f, a, mid, x, w; tol = tol, levels= levels -1) .+ RecursiveIntegral(f, mid, b, x, w; tol = tol, levels= levels -1)
    end


end


function RecursiveIntegralDetail(f, a, b, x, w; tol = 1e-10, levels = 10, subdivisions = Float64[a,b])
    scale = (b-a)
    mid = (a+b)/2
    append!(subdivisions, mid)
    Iₜ = sum(scale*w .* f.(a .+scale*x))
    I₁ = sum(scale/2*w .* f.(a .+ (scale/2)*x))
    I₂ = sum(scale/2*w .* f.(mid .+ (scale/2)*x))
    # println(levels) 
    if abs(I₁ + I₂ - Iₜ) < tol || levels == 0
        return I₁ + I₂, subdivisions
    else
        return RecursiveIntegralDetail(f, a, mid, x, w; tol = tol, levels = levels -1, subdivisions = subdivisions)[1] + RecursiveIntegralDetail(f, mid, b, x, w; tol = tol, levels = levels -1,subdivisions = subdivisions)[1], subdivisions
    end


end


# function RecursiveIntegral(f, a, b, x, w; tol = 1e-8, levels = 3)
#     scale = (b-a)
#     mid = (a+b)/2
#     Iₜ = sum(scale*w .* f.(a .+scale*x))
#     I₁ = sum(scale/2*w .* f.(a .+ (scale/2)*x))
#     I₂ = sum(scale/2*w .* f.(mid .+ (scale/2)*x))

#     # println("Iₜ: ", Iₜ)
    
#     # println(I₁ + I₂)
#     if abs(I₁ + I₂ - Iₜ) < tol || levels == 0
#         return I₁ + I₂, abs(I₁ + I₂ - Iₜ) 
#     else
#         return RecursiveIntegral(f, a, mid, x, w; tol = tol, levels= levels -1) .+ RecursiveIntegral(f, mid, b, x, w; tol = tol, levels= levels -1)
#     end


# end

function RecursiveIntegralMesh(f, a, b, x, w; tol = 1e-8, levels = 3, nodes = Float64[], weights = Float64[])
    scale = (b-a)
    mid = (a+b)/2
    Iₜ = sum(scale*w .* f.(a .+scale*x))
    I₁ = sum(scale/2*w .* f.(a .+ (scale/2)*x))
    I₂ = sum(scale/2*w .* f.(mid .+ (scale/2)*x))

    # println("Iₜ: ", Iₜ)
    
    # println(I₁ + I₂)
    if abs(I₁ + I₂ - Iₜ) < tol || levels == 0
        append!(nodes, a .+ (scale/2)*x)
        append!(weights, scale/2*w)
        append!(nodes, mid .+ (scale/2)*x)
        append!(weights, scale/2*w)
        return I₁ + I₂, nodes, weights
    else
        return (RecursiveIntegralMesh(f, a, mid, x, w; tol = tol, levels= levels -1, nodes = nodes, weights= weights)[1] 
        + RecursiveIntegralMesh(f, mid, b, x, w; tol = tol, levels= levels -1, nodes = nodes, weights= weights)[1]), nodes, weights
    end

end