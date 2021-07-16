using LinearAlgebra, Optim, ForwardDiff

##############
# qubit bound
##############

function qbitbound(z,x²,q)
    # pay attention that the function takes as argument x²=x*x and not x
    rp=sqrt((1-2*q+z)^2+4*q*(1-q)*x²)
    rm=sqrt((1-2*q-z)^2+4*q*(1-q)*x²)
    φ(1/2*(rp+rm))+φ(1/2*(rp-rm))-φ(sqrt(z^2+x²))
end

function sqbitbound(a,s,q)
    qbitbound(a,stox²(s),q)
end

function stox²(s)
    s<=2 ? 0 : s^2/4-1
end

function sqbitbound_noa(a,s,q)
    qbitbound(0,stox²(s),q)
end


function smax(a)
    2*sqrt(2-a^2)
end


############################################
# convexified bound based on optimal attack
############################################

function convexbound(a,s,q)
    #=
    compute entropy bound h at the point (a,s) as a convex combination of the point (1,2) and the point ((1-x)*1+x*a,(1-x)*2+x*s) for the optimal value x 
    return (h,x)
    =#
    a=abs(a)
    # putting s<smax(a) instead of s<=smax(a) in the line below avoid that xmax<1 because of numerical accuracy, which would cause optimize to return an error
    if 2<s<smax(a)
        if a <= (2*sqrt(2)-s)/(2*sqrt(2)-2)
            xmax=1/(1-a)
        else
            xmax=4*(4-2*a-s)/((s-2)^2+(2a-2)^2)
        end
        res=Optim.optimize(x->convexboundpointx(x,a,s,q),1,xmax)
        Optim.minimum(res), Optim.minimizer(res)
    else
        0,1 # outputting 0 instead of squbitbound(a,s,q) helps the optimization over quantum strategies
    end
end

function convexboundpointx(x,a,s,q)
    (1-1/x)*sqbitbound(1,2,q)+1/x*sqbitbound((1-x)+x*a,(1-x)*2+x*s,q)
end

function convexboundprimal(a,s,q)
    #=
     like convexbound(a,s,q) but return the explicit convex decomposition yielding the entropy bound, i.e.,
     return (h,t,p1,p2) s.t. p = t[1].*p1 .+ t[2].*p2
     where p=(a,s,h)
    =# 
    (h,x)=convexbound(a,s,q)
    ax=1-x+x*a
    sx=(1-x)*2+x*s
    h,(1-1/x,1/x),[1; 2; sqbitbound(1,2,q)],[ax; sx; sqbitbound(ax,sx,q)]
end

function convexbounddual(a,s,q)
    #=
    return h,w where h is convexified entropy bound and w=(w1,w2,w3) is a vector that certifies the convex entropy bound :
    h(a',s')>=w1+w2*a'+w3*s' for all quantum (a',s')
    and
    h(a,s)=w1+w2*a+w3*s for the point (a,s) of interest
    =#
    (h,t,p1,p2)=convexboundprimal(a,s,q)
    grad = ForwardDiff.gradient(x -> sqbitbound(x...,q),p2[1:2])
    if t[2]<1-1e-7 # this choice is a bit arbitrary
        v = p1[1:2] - p2[1:2]
        v = v / norm(v)
        u=[-v[2]; v[1]]
        A=[1 p1[1:2]'; 1 p2[1:2]'; 0 u[1:2]']
        b=[p1[3]; p2[3]; grad'*u]
    else
        A=[1 p2[1:2]'; 0 1 0; 0 0 1]
        b=[p2[3]; grad[1:2]]
    end
        h,A \ b
end

function certify(a,s,q,delta,hbound)
    # certify that hbound is a valid bound on entropy and return certificate w
    (htrue,w)=convexbounddual(a,s,q)
    iswvalidonrect(w,(0,1),(2,2*sqrt(2)),htrue-hbound,delta,q),w
end

function certify(a,s,q,delta,hbound,w)
    # certify that hbound is a valid bound on entropy for given dual vector w
    htrue=w[1]+w[2]*a+w[3]*s
    iswvalidonrect(w,(0,1),(2,2*sqrt(2)),htrue-hbound,delta,q)
end


function iswvalidonrect(w,abounds,sbounds,tol,delta,q)
    #=
    check that h(a,s)>=w[1]+w[2]*a+w[3]*s-tol holds for all values of a,s in the rectangle defined by abounds=(amin,amax) and sbounds=(smin,smax)
        w=(w1,w2,w3) is a candidate entropy certificate 
        delta is the minimum subrectangle size
    =#
    if abounds[2]-abounds[1]<=delta || sbounds[2]-sbounds[1]<=delta
        false
    # ok if the rectangle is outside the quantum region
    elseif sbounds[1]>smax(abounds[1])
        true
    elseif w[1]+w[2]*abounds[2]+w[3]*sbounds[2]-sqbitbound(abounds[1],sbounds[1],q) <= tol
        true
    else
        aint=(abounds[1]+abounds[2])/2
        sint=(sbounds[1]+sbounds[2])/2
        iswvalidonrect(w,(aint,abounds[2]),(sint,sbounds[2]),tol,delta,q) &&
        iswvalidonrect(w,(aint,abounds[2]),(sbounds[1],sint),tol,delta,q) &&
        iswvalidonrect(w,(abounds[1],aint),(sint,sbounds[2]),tol,delta,q) &&
        iswvalidonrect(w,(abounds[1],aint),(sbounds[1],sint),tol,delta,q)
    end
end

function iswvalidonrect_cover!(cover,w,abounds,sbounds,tol,delta,q)
    #=
    as previous function but also return the rectangle covering lattice
    - cover is a n element vector of 5-element vectors, where each 5-element
        vector encodes a rectangle as [amin,amax,smin,smax,hdif]
        where hdif is how much certificate bound is above qubit bound
        if calling for the first time define cover=Vector{Float64}[]
    =#
    if abounds[2]-abounds[1]<=delta || sbounds[2]-sbounds[1]<=delta
        false
    # ok if the rectangle is outside the quantum region
    elseif sbounds[1]>smax(abounds[1])
        true
    else
        dif=w[1]+w[2]*abounds[2]+w[3]*sbounds[2]-sqbitbound(abounds[1],sbounds[1],q)
        if dif<= tol
            push!(cover,[abounds[1],abounds[2],sbounds[1],sbounds[2],dif])
            true
        else
            aint=(abounds[1]+abounds[2])/2
            sint=(sbounds[1]+sbounds[2])/2
            iswvalidonrect_cover!(cover,w,(aint,abounds[2]),(sint,sbounds[2]),tol,delta,q) &&
            iswvalidonrect_cover!(cover,w,(aint,abounds[2]),(sbounds[1],sint),tol,delta,q) &&
            iswvalidonrect_cover!(cover,w,(abounds[1],aint),(sint,sbounds[2]),tol,delta,q) &&
            iswvalidonrect_cover!(cover,w,(abounds[1],aint),(sbounds[1],sint),tol,delta,q)
        end
    end
end

# ####################################################
# # convexified bound based on a discrete grid and LP
# ####################################################
using JuMP, Gurobi
# using JuMP, GLPK

function gridpoints(w,q)
    # define non-rectangular grid with interpoint space w
    av=Float64[]
    sv=Float64[]
    hv=Float64[]
    a=0
    s=2
    while a<=1
        h=sqbitbound(a,s,q)
        if ~isnan(h)
            push!(av,a)
            push!(sv,s)
            push!(hv,h)
            s+=w
        else
            if last(sv)!=smax(a)
                push!(av,a)
                push!(sv,smax(a))
                push!(hv,sqbitbound(a,smax(a),q))
            end
            a+=w
            s=2
        end
    end
    hcat(av,sv,hv)
end

function extreme_grid_points(g)
    # v-representation based on gridpoints and ray [0 0 1]
    #p=polyhedron(vrep(g,[0 0 1]),CDDLib.Library()) 
    # v-representation based on gridpoints and no ray (not supported by QHull)
    p=polyhedron(vrep(g),QHull.Library()) 
    
    removevredundancy!(p)
    [i[j] for i in points(p), j=1:3]
end

function LP_model(g)
    # define LP model for JuMP based on grid points 
    model=Model(Gurobi.Optimizer)
    set_optimizer_attribute(model,"OutputFlag",0)
    # model=Model(GLPK.Optimizer)
    (n,)=size(g)
    @variable(model,t[1:n],lower_bound=0.0)
    @variable(model,corr[1:2])
    @constraint(model,norma,sum(t)==1)
    @constraint(model,rep,g[:,1:2]'*t .==corr)
    @objective(model,Min,g[:,3]'*t)
    model
end

function convexbound_LP(a,s,LPmodel)
    #= solve LP to determine convexified entropy for given p=[a,s]
        returns convex entropy value, the indices of extreme points contributing to convex decomposition, and weights of such extreme points
        requires JuMP "LPmodel" constructed with LP_model(g) specifying the grid points g
    =#
    fix(LPmodel[:corr][1],a)
    fix(LPmodel[:corr][2],s)
    optimize!(LPmodel)
    h=objective_value(LPmodel)
    ind=findall(value.(LPmodel[:t]).>=1e-6) 
    h,ind,value.(LPmodel[:t])[ind] 
end

function convexboundprimal_LP(a,s,LPmodel,g)
    (h,ind,t)=convexbound_LP(a,s,LPmodel)
    h,t,g[ind,:]
end

function gridpointsshifted(w,q)
    # define shifted non-rectangular grid with interpoint space w
    av=Float64[]
    sv=Float64[]
    hv=Float64[]
    a=0
    s=2
    while a<1
        a==0 ? h=sqbitbound(a,s-w,q) : h=sqbitbound(a-w,s-w,q)
        if ~isnan(h)
            push!(av,a)
            push!(sv,s)
            push!(hv,h)
            s+=w
        else
            a+=w
            s=2
        end
    end
    push!(av,1)
    push!(sv,2)
    push!(hv,sqbitbound(1,2,q))
    hcat(av,sv,hv)
end


##########################################
#=
These are a few method definitions so that we can call every bound as
    bound(a,s,q,param)
for instance in keyrate()
=#

function sqbitbound(a,s,q,param::Nothing)
    sqbitbound(a,s,q)
end

function sqbitbound_noa(a,s,q,param::Nothing)
    sqbitbound_noa(a,s,q)
end

function convexbound(a,s,q,param::Nothing)
    convexbound(a,s,q)
end

function convexbound_LP(a,s,q,LPmodel)
    convexbound_LP(a,s,LPmodel)
end