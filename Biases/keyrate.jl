function keyrate_corr(a1,s,Pa1b3,q,bound,param=nothing)
    #=
    compute the keyrate for given correlations using a given entropy bound
    return k,hae,hab
    - choices for bound = convexbound, sqbitbound_noa, sqbitbound
    =#
    (hae,)=bound(a1,s,q,param)
    Pb3=Pa1b3[1] .+ Pa1b3[2]
    Pa1b3 = (Pa1b3[1]...,Pa1b3[2]...) # we redefine Pa1b3 as a tuple of values
    hab=hcond(Pa1b3,Pb3)
    hae-hab,hae,hab
end

function keyrate(θ,α1,α2,β1,β2,β3,v,η,q,bound,param=nothing)
    #=
    compute the keyrate for given quantum strategy using a given entropy bound
    return k,hae,hab,a1,s
    - choices for bound = convexbound, sqbitbound_noa, sqbitbound
    =#
    a1=E1(θ,α1,v,η)
    s=S(θ,α1,α2,β1,β2,v,η)
    Pa1b3=Pab(θ,α1,β3,v,η,q)
    keyrate_corr(a1,s,Pa1b3,q,bound,param)...,a1,s
end

function keyrate(θ,α1,α2,β3,v,η,q,bound,param=nothing)
    #=
    as above but with optimal βs
    =#
    a1=E1(θ,α1,v,η)
    s=S(θ,α1,α2,v,η)
    Pa1b3=Pab(θ,α1,β3,v,η,q)
    keyrate_corr(a1,s,Pa1b3,q,bound,param)...,a1,s
end

function keyrate(θ,α2,v,η,q,bound,param=nothing)
    #=
    as above but with α1=β3=0
    =#
    keyrate(θ,0,α2,0,v,η,q,bound,param)
end

function keyrate_q(θ,α2,q,v,η,bound,param=nothing)
    #=
    as above but with α1=β3=0
    =#
    keyrate(θ,0,α2,0,v,η,q,bound,param)
end

function find_best_qstrat(x0,v,η,q,bound,param=nothing)
    #=
        example of starting point x0 assuming a1=b3=Z
        x0=[0.4 pi/2]
        example of starting point x0 for a1,b3 not equal
        x0=[0.4 0 pi/2 0.3]
        return k,hae,hab,a1,s,optimal x
    =#
        res=Optim.optimize(x->-keyrate(x...,v,η,q,bound,param)[1],x0)
        output=keyrate(Optim.minimizer(res)...,v,η,q,bound,param)
        output...,Optim.minimizer(res)
end

function find_best_qstrat_q(x0,v,η,bound,param=nothing)
    #=
        example of starting point x0 assuming a1=b3=Z
        x0=[0.4 pi/2]
        example of starting point x0 for a1,b3 not equal
        x0=[0.4 0 pi/2 0.3]
        return k,hae,hab,a1,s,optimal x
    =#
        lower=[0 0 0]
        upper=[2*pi 2*pi 0.499]
        res=Optim.optimize(x->-keyrate_q(x...,v,η,bound,param)[1],lower,upper,x0)
        output=keyrate_q(Optim.minimizer(res)...,v,η,bound,param)
        output...,Optim.minimizer(res)
end

function etathresh(x,v,q,bound,param=nothing)
    δ = 1e-6 # the value we consider a zero keyrate
    (k,)=keyrate(x...,v,1,q,bound,param)
    if k<δ
        return 1
    end
    tol=1e-16
    step=0.1
    eta=1
    while step>tol
        etalow=eta-step
        if keyrate(x...,v,etalow,q,bound,param)[1]>δ
            eta=etalow
        else
            step=step/2
        end
    end
    eta
end

function etathresh_n(x,v,q,δ,bound,param=nothing) #works better in the noisy cases
    #δ = 1e-11 # the value we consider a zero keyrate
    (k,)=keyrate(x...,v,1,q,bound,param)
    if k<δ
        return 1
    end
    tol=1e-16
    step=0.1
    eta=1
    while step>tol
        etalow=eta-step
        if keyrate(x...,v,etalow,q,bound,param)[1]>δ
            eta=etalow
        else
            step=step/2
        end
    end
    eta
end

function find_best_eta(x0,v,q,bound,param=nothing)
    res=Optim.optimize(x->etathresh(x,v,q,bound,param),x0)
    eta=Optim.minimum(res)
    output=keyrate(Optim.minimizer(res)...,v,eta,q,bound,param)
    eta,output...,Optim.minimizer(res)
end

function find_best_eta_n(x0,v,q,δ,bound,param=nothing)
    res=Optim.optimize(x->etathresh_n(x,v,q,δ,bound,param),x0)
    eta=Optim.minimum(res)
    output=keyrate(Optim.minimizer(res)...,v,eta,q,bound,param)
    eta,output...,Optim.minimizer(res)
end


function certifykeyrate(θ,α1,α2,β1,β2,β3,v,η,q,delta,kr)
    #=
    certify that the keyrate kr is indeed valid and return certificate w
    =#
    (k,hae,hab,a1,s)=keyrate(θ,α1,α2,β1,β2,β3,v,η,q,convexbound)
    certify(a1,s,q,delta,kr+hab)
end

function certifykeyrate(θ,α1,α2,β3,v,η,q,delta,kr)
    (k,hae,hab,a1,s)=keyrate(θ,α1,α2,β3,v,η,q,convexbound)
    certify(a1,s,q,delta,kr+hab)
end

function certifykeyrate(θ,α2,v,η,q,delta,kr)
    certifykeyrate(θ,0,α2,0,v,η,q,delta,kr)
end

