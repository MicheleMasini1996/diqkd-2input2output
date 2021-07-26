using Optim, ForwardDiff, Roots, Plots, LinearAlgebra
using DynamicPolynomials, SumOfSquares, MosekTools

function xlog2x(x)
    if x<0 || x>1
        0
    elseif x==0 || x==1
        0
    else
        -x*log2(x)
    end
end

function hshannon(prob) # Shannon entropy
    sum(xlog2x.(prob))
end

function hbin(x) # binary entropy
    hshannon((x,1-x))
end

function phi(x)
    hbin(1/2+x/2)    
end

function fq(x,q) # function f_q defined in the paper
    1+phi(sqrt((1-2*q)^2+4*q*(1-q)*x))-phi(sqrt(x))
end

function E(S) # analytical solution for E^2 when p=0.5
    if S<=2
        return 0
    else
        S=complex(S)
        x=(96 + 3*S^2 - (6144 + 384*S^2 + 9*S^4 + (384*(1024 - 352*S^2 + 49*S^4))/(262144 - 135168*S^2 + 9264*S^4 + 928*S^6 + 27*S^8 + 3*sqrt(3)*(S^4*(112 - 32*S^2 + S^4)^2*(-32768 + 3584*S^2 + 27*S^4))^0.5)^(1/3) + 96*(262144 - 135168*S^2 + 9264*S^4 + 928*S^6 + 27*S^8 + 3*sqrt(3)*(S^4*(112 - 32*S^2 + S^4)^2*(-32768 + 3584*S^2 + 27*S^4))^0.5)^(1/3))^0.5 + sqrt(6)*(-1024 - 64*S^2 + 3*(32 + S^2)^2 - (64*(1024 - 352*S^2 + 49*S^4))/(262144 - 135168*S^2 + 9264*S^4 + 928*S^6 + 27*S^8 + 3*sqrt(3)*(S^4*(112 - 32*S^2 + S^4)^2*(-32768 + 3584*S^2 + 27*S^4))^0.5)^(1/3) - 16*(262144 - 135168*S^2 + 9264*S^4 + 928*S^6 + 27*S^8 + 3*sqrt(3)*(S^4*(112 - 32*S^2 + S^4)^2*(-32768 + 3584*S^2 + 27*S^4))^0.5)^(1/3) - (3*sqrt(3)*S^2*(32 + S^2)^2)/(2048 + 128*S^2 + 3*S^4 + (128*(1024 - 352*S^2 + 49*S^4))/(262144 - 135168*S^2 + 9264*S^4 + 928*S^6 + 27*S^8 + 3*sqrt(3)*(S^4*(112 - 32*S^2 + S^4)^2*(-32768 + 3584*S^2 + 27*S^4))^0.5)^(1/3) + 32*(262144 - 135168*S^2 + 9264*S^4 + 928*S^6 + 27*S^8 + 3*sqrt(3)*(S^4*(112 - 32*S^2 + S^4)^2*(-32768 + 3584*S^2 + 27*S^4))^0.5)^(1/3))^0.5)^0.5)/96.
        E=(x^2+1)/(1-x)+S^2/4*(1+x)/(1-x)-S*(1+x)/(1-x)*((1+x)/2)^0.5
        return E.re
    end
end

function Ep(p,S) # finds E_p^2 using Lasserre
    @polyvar s c l m d 
    obj=s^2*l^2+c^2*m^2+2*(2*p-1)*s*c*l*m*d
    constr = @set c*l+s*m>=S/2 && l^2<=1 && m^2<=1 && (1-l^2)*(1-m^2)>=l^2*m^2*d^2 && c^2+s^2==1 && d^2<=1
    solver = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
    model = SOSModel(solver)
    @variable(model,α)
    @objective(model,Max,α)
    @constraint(model,c5, obj>=α, domain = constr, maxdegree = 6)
    optimize!(model)
    objective_value(model)
end

function rn(E,q) # key rate with noisy pre-processing and depolarizing noise
    dq=q+(1-2*q)*delta
    fq(E.real,q)-h(dq)
end

function shift(p,q)	# computes the entropy on shifted samples
    apts=10000
    pts=500
    if p==1
        samples=zeros(apts+1,2)
        S=2.
        samples[1,1]=S
        samples[1,2]=fq(S^2/4-1,q)
        for n = 2:apts+1
            samples[n,2]=fq(S^2/4-1,q)
            S=2+n*(2*(2^0.5)-2)/(apts+1)
            samples[n,1]=S
        end
        return samples
    elseif p==0.5
        samples=zeros(apts+1,2)
        S=2.
        samples[1,1]=S
        samples[1,2]=fq(E(S),q)
        for n = 2:apts+1
            samples[n,2]=fq(E(S),q)
            S=2+n*(2*(2^0.5)-2)/(apts+1)
            samples[n,1]=S
        end
        return samples
    else
        samples=zeros(pts+1,2)
        S=2.
        samples[1,1]=S
        samples[1,2]=fq(Ep(p,S),q)
        for n = 2:pts+1
            samples[n,2]=fq(Ep(p,S),q)
            S=2+n*(2*(2^0.5)-2)/(pts+1)
            samples[n,1]=S
        end
        return samples
    end
end

function chull(sam) # convexifies the function built on the samples
    X=sam[:,1]
    F=sam[:,2]
    N = min(length(X), length(F))

    stack = []
    k = 1

    while k <= N
        if length(stack) < 2
            push!(stack, k)
            k += 1
        else
            k1 = stack[end]
            k0 = stack[end-1]
            slope1 = (F[k] - F[k1]) / (X[k] - X[k1])
            slope0 = (F[k1] - F[k0]) / (X[k1] - X[k0])

            if slope1 > slope0
                push!(stack, k)
                k += 1
            else
                pop!(stack)
            end
        end
    end

    return [X[stack] F[stack]]
end

function centropy(s,p,q) # computes a lower bound on the conditional entropy on a single point
    if p==1
    	interp=shift(p,q)
    else
    	interp=chull(shift(p,q))
    end
    dim=size(interp[:,1])[1]
    ch=zeros(dim)
    for i=1:dim
        ch[i]=(abs(interp[i,1]-s))
    end
    a=argmin(ch)
    if a==1
    	mi=1
    	ma=2
    elseif a==dim
    	mi=dim-1
    	ma=dim
    elseif interp[a,1]-s>=0
        mi=a-1
        ma=a
    else
        mi=a
        ma=a+1
    end
    interp[mi,2]+(interp[ma,2]-interp[mi,2])/(interp[ma,1]-interp[mi,1])*(s-interp[mi,1])
end


function rate(delta,p,q) # computes the key rate as a function of the channel error rate for a single point
    if p==0.5
        pp=0.5
    elseif p==1
        pp=1
    else
        pp=(p - sqrt(p - p^2))/(-1 + 2*p)
    end
    deltaq=q+(1-2*q)*delta
    s=2*sqrt(2)*(1-2*delta)
    if p==1
        return centropy(s,p,q)-hbin(deltaq)
    else
        return (pp^2+(1-pp^2))*(centropy(s,p,q)-hbin(deltaq))
    end
end
	
function plot_centropy(p,q) # plots the conditional entropy between Alice and Eve as a function of S
    sam=shift(p,q)
    cv=chull(sam)
    gr(size=(700,500), html_output_format=:png)
    plot(cv[:,1],cv[:,2])
end

function plot_rate(p,q) # plots the key rate as a function of the channel error rate 
    sam=shift(p,q)
    cv=chull(sam)
    gr(size=(700,500), html_output_format=:png)
    if p==0.5
        pp=0.5
    elseif p==1
        pp=1
    else
        pp=(p - sqrt(p - p^2))/(-1 + 2*p)
    end
    dim=length(cv[:,1])
    prov=zeros(dim,2)
    count=0
    for i = 1:dim
        delta=0.5*(1-cv[dim-i+1,1]/(2*sqrt(2)))
        deltaq=q+(1-2*q)*delta
        prov[i,1]=delta
        if p==1
            prov[i,2]=(cv[dim-i+1,2]-hbin(deltaq))
        else
            prov[i,2]=(pp^2+(1-pp^2))*(cv[dim-i+1,2]-hbin(deltaq))
        end
        count=count+1
        if prov[i,2]<0
            break
        end
    end
    pl=zeros(count,2)
    for i=1:count
        pl[i,:]=prov[i,:]
    end
    plot(pl[:,1],pl[:,2])
end
