using LinearAlgebra, Optim, ForwardDiff

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
