#=
Compute quantum correlations assuming
    - a state cos(θ)|00⟩+sin(θ)|11⟩
    - mesurements cos(α)Z+sin(α)X for Alice, cos(β)Z+sin(β)X for Bob 
    - visibility v
    - detector efficiency η
    - in case of E1, E2, S: binning of (-1,+1,∅) to (-1,+1) where ∅ → +1
    - in case of PAB: only binning of Alice's no detection outcomes, Bob's ones
      are kept 3-valued, and moreover Alice's resulting two outcomes are noisy preprocessed with prob q
=#

function E1(θ,γ)
    # 1-body correlator
    cos(γ)*cos(2θ)
end

function E1(θ,γ,v,η)
    # noisy 1-body correlator
    η*v*E1(θ,γ)+(1-η)
end

function E2(θ,α,β)
    # 2-body correlator
    cos(α)*cos(β)+sin(α)*sin(β)*sin(2θ)
end

function E2(θ,α,β,v,η)
    # noisy 2-body correlator
    η^2*v*E2(θ,α,β)+η*(1-η)*v*E1(θ,α)+η*(1-η)*v*E1(θ,β)+(1-η)^2
end

function S(θ,α1,α2,β1,β2,v,η)
    # noisy CHSH
    E2(θ,α1,β1,v,η)+E2(θ,α1,β2,v,η)+E2(θ,α2,β1,v,η)-E2(θ,α2,β2,v,η)
end

function S(θ,α1,α2,v,η)    
    # noisy CHSH with optimal βs
    p = η*(cos(α1)+cos(α2))+2*(1-η)*cos(2*θ)
    q = η*(sin(α1)+sin(α2))*sin(2*θ)
    r = cos(α1)-cos(α2)
    s = (sin(α1)-sin(α2))*sin(2*θ)
    v*η*sqrt(p^2+q^2)+v*η^2*sqrt(r^2+s^2)+v*2*η*(1-η)*cos(α1)*cos(2*θ)+2*(1-η)^2
end

function Pab(θ,α,β,v,η,q)
    # joint probabilities
    A = E1(θ,α,v,1)
    B = E1(θ,β,v,1)
    AB = E2(θ,α,β,v,1)
    # full 3x3 probability table
    P = ((η^2*(1+A+B+AB)/4, η^2*(1+(A-B)-AB)/4, η*(1-η)*(1+A)/2),
    (η^2*(1-(A-B)-AB)/4, η^2*(1-A-B+AB)/4, η*(1-η)*(1-A)/2),
    (η*(1-η)*(1+B)/2, η*(1-η)*(1-B)/2, (1-η)^2))
    # binning of Alice's outcomes
    P = (P[1] .+ P[3], P[2])
    # noisy preprocessing
    ((1-q).*P[1] .+ q.*P[2], (1-q).*P[2] .+ q.*P[1])
end