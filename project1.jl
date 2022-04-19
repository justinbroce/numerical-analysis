using LinearAlgebra, Symbolics
#=
Secant and newton methods, taken from:
Ökten, G. (2019). First Semester in Numerical Analysis with Julia. Retrieved from http://purl.flvc.org/fsu/fd/FSU_libsubv1_scholarship_submission_1556028278_15938059
=#
function secant(f::Function,pzero,pone,eps,N)
    n=1
    p=0. # to ensure the value of p carries out of the while loop
    while n<=N
        p=pone-f(pone)*(pone-pzero)/(f(pone)-f(pzero))
        if f(p)==0 || abs(p-pone)<eps
            return println("p is $p and the iteration number is $n")
        end
        pzero=pone
        pone=p
        n=n+1
    end
    y=f(p)
end

function newton(f::Function,fprime::Function,pin,eps,N)
    n=1
    p=0. # to ensure the value of p carries out of the while loop
    while n<=N
        p=pin-f(pin)/fprime(pin)
        if f(p)==0 || abs(p-pin)<eps
            return println("p is $p and the iteration number is $n")
        end
        pin=p
        n=n+1
    end
    y=f(p)
    println("Method did not converge. The last iteration gives$p with function value $y")
end
#lugauss adapted from our textbook:
function lugauss(A)
    n, m = size(A)
     if n != m
        print("A is not a square matrix")
        return -1.0
    end
    for k in 1:n-1 
        for i in k+1:n
            A[i,k] = A[i,k]/A[k,k]
            if isequal(A[k,k], 0)
               print("Null diagonal element")
            return -1.0
            end 
            j = k+1:n
            A[i,j] = A[i,j] - A[i,k] * A[k,j]
        end
    end
end