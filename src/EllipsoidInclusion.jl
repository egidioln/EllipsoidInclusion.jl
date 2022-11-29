

function Base.:∉(elli1::Ellipsoid, elli2::Ellipsoid)
    !(elli1 ∈ elli2)
end


function Base.in(x::AbstractVecOrMat, elli::Ellipsoid)
    return (x-elli.c)'elli.P*(x-elli.c) ≤ 1
end

function Base.in(elli1::Ellipsoid, elli2::Ellipsoid)
    e_min = eigmin(elli1.P-elli2.P)
    if e_min<0
        return false
    elseif  elli1.c==elli2.c
        return e_min>=0
    elseif !(elli1.c ∈ elli2)
        return false
    else 
        L = cholesky((elli2.P+elli2.P')/2).L #TODO remove ' when Ellipsoid constructor fixed
        P = L\elli1.P/L';
        c = L'*(elli1.c -elli2.c)
        specDecomp = eigen(P)
        lb = specDecomp.values
        ct = specDecomp.vectors'*c
        α = 1/min(lb...)

        polPos(β) = -(1-β + sum((β*lb./(1 .- β*lb)).*(ct.^2)))
        dpolPos(β) = 1 - sum((lb./(1 .- β*lb).^2).*(ct.^2))
        ddpolPos(β) = 2*sum((lb.^2 ./(1 .- β*lb).^3).*(ct.^2))
        if(α==1)
            return polPos(1)<=0
        end
        if(α>1-norm(ct)^2)
            return false
        end
        (val, _) = dbisection(polPos, dpolPos, ddpolPos, interval=[α+1e-15, 1-norm(ct)^2], verbose=false, stopIfNegative=true)

        return val<=0
    end
end

function get_ℓ_ast(elli1::Ellipsoid, elli2::Ellipsoid)
    L = cholesky((elli2.P+elli2.P')/2).U #TODO remove ' when Ellipsoid constructor fixed
    P = L'\elli1.P/L;
    specDecomp = eigen(P)
    lb = specDecomp.values
    ct = specDecomp.vectors'*L*(elli1.c -elli2.c)
    α = 1/min(lb...)

    polPos(β) = -(1-β +sum((β*lb./(1 .- β*lb)).*(ct.^2)))
    dpolPos(β) = 1 - sum((lb./(1 .- β*lb).^2).*(ct.^2))
    ddpolPos(β) = -2*sum((lb.^2 ./(1 .- β*lb).^3).*(ct.^2))

    ub = α*2;
    while dpolPos(ub) < 0
        ub*=2
    end
    (val, β) = dbisection(polPos, dpolPos, ddpolPos, interval=[α+1e-15, ub], verbose=false, stopIfNegative=false)

    return -val-1, β
end


# minimise a convex function
function dbisection(f::Function,df::Function,ddf::Function; interval=[0, 1],  δ=1e-8, verbose=false, stopIfNegative=false)

    x = zeros(2)
    fval = zeros(2)
    dfval = zeros(2)
    L = 0

    x[1] = interval[1];
    x[2] = interval[2];
    
    fval[1] = f(x[1]);
    fval[2] = f(x[2]);
    dfval[1] = df(x[1]);
    dfval[2] = df(x[2]);
    L = ddf(x[1]);
    if dfval[2] < 0
        return (fval[2],x[2])
    end

    fβ, β = let β = (x[1] + x[2])/2, fβ = f(β)
        dfβ = df(β)
        if (verbose)
            println(string(dfval[1])*" "*string(dfval[2])*"    "*string(fβ))
        end
        while abs(x[1]-x[2])>δ && (!stopIfNegative || (fβ>0 && 2*fβ<-L*(x[2]-x[1])^2))
            
            if (verbose)
                println(string(x[1])*"\t"*string(x[2])*"\t"*string(fβ))
                println(string(dfval[1])*"\t"*string(dfval[2])*"\t"*string(fβ))
            end
            if(dfβ <0)
                x[1] = β;
                dfval[1] = dfβ;
                L = ddf(β);
            else
                x[2] = β;
                dfval[2] = dfβ;
            end

            β = (x[1] + x[2])/2
            fβ = f(β)
            dfβ = df(β)
        end
        (fβ, β)
    end
    (fβ, β)
end

