using LinearAlgebra, Plots, Statistics, Random

function jacobi(A,b,ε=1e-12,MAX_ITER=1000)
    k=0
    x₀ = zeros(rank(A))
    r₀ = b - A*x₀
    A = Matrix(A)
    B = A-Diagonal(A)
    D⁻¹ = inv(Diagonal(A))
    
    while norm(r₀)>ε && MAX_ITER > k
        x₀ = D⁻¹*(b-B*x₀)
        r₀ = b - A*x₀
        k = k+1
    end
    println("$k iterations")
    return x₀
    
end
function cg(A, b, ε=1e-2, MAX_ITER = 1000)
    k  = 0
    x₀ = zeros(rank(A))
    r₀ = b - A*x₀
    p₀ = r₀
    xₖ₊₁=x₀
    while norm(r₀) > ε && MAX_ITER > k
        αₖ   = r₀⋅r₀/(p₀⋅(A*p₀))
        xₖ₊₁ = x₀ + αₖ*p₀
        rₖ₊₁ = r₀ - αₖ*A*p₀
        βₖ   = (rₖ₊₁⋅rₖ₊₁)/(r₀⋅r₀)
        pₖ₊₁ = rₖ₊₁ + βₖ*p₀
        #setting the previous set
        k  = k+1
        p₀ = pₖ₊₁
        x₀ = xₖ₊₁
        r₀ = rₖ₊₁
    end
    println("$k iterations")
    return xₖ₊₁
    #return k
end
function SSOL(A, ω)
 D = Diagonal(A)
 L = LowerTriangular(A)
 U = UpperTriangular(A)
 (D/ω + L)*ω/(2-ω)*inv(D)*(D/ω+transpose(L))
end
function GSp(A)
     D = Diagonal(A)
     L = LowerTriangular(A)
     U = UpperTriangular(A)
     (D+L)*inv(D)*(D+transpose(L))
end
function pcg(A,b,B=Diagonal(A),ε=1e-12,MAX_ITER=1000)
    k  = 0
    M= inv(B)
    x₀  = zeros(rank(A))
    r₀  = b - A*x₀
    z₀  = M*r₀
    v₀  = z₀
    xₖ₊₁= x₀
    while norm(A*x₀-b) > ε && MAX_ITER > k
        αₖ   = r₀⋅z₀/((A*v₀)⋅v₀)
        xₖ₊₁ = x₀ + αₖ*v₀
        rₖ₊₁ = r₀ - αₖ*A*v₀
        zₖ₊₁ = M*rₖ₊₁
        βₖ   = (rₖ₊₁⋅zₖ₊₁)/(r₀⋅z₀)
        vₖ₊₁ = zₖ₊₁ + βₖ*v₀
        #setting the previous set
        k  = k+1
        v₀ = vₖ₊₁
        z₀ = zₖ₊₁
        x₀ = xₖ₊₁
        r₀ = rₖ₊₁
    end
    println("$k iterations")
    return xₖ₊₁
    #return k
end