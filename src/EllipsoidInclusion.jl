module EllipsoidInclusion
using LinearAlgebra
using SpecialFunctions


struct Ellipsoid{T<:Real,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
    P::MT
    c::VT
    function Ellipsoid(P::MT,c::VT) where {T<:Real,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
        P_ = (P+P')
        if eigmin(P_) > 0
            new{T,MT,VT}(P_./2,c)
        else
            error("P must be a positive definite matrix")
        end

    end
end


function Base.:*(elli::Ellipsoid, r::Real)
    Ellipsoid(elli.P*(1/r), elli.c)
end

function Base.:*(r::Real, elli::Ellipsoid)
    elli*r
end

function Base.:∉(elli1::Ellipsoid, elli2::Ellipsoid)
    !(elli1 ∈ elli2)
end


function Base.in(elli1::Ellipsoid, elli2::Ellipsoid)
    e_max = eigmax(elli1.P-elli2.P)
    if e_max<0
        return false
    elseif  elli1.c==elli2.c
        return e_max>=0
    elseif !(elli1.c ∈ elli2)
        return false
    else 
        L = cholesky((elli2.P+elli2.P')/2).U #TODO remove ' when Ellipsoid constructor fixed
        P = L'\elli1.P/L;
        specDecomp = eigen(P)
        lb = specDecomp.values
        ct = specDecomp.vectors'*L*(elli1.c -elli2.c)
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
        #(val, _) = bisection(polPos, interval=[α+1e-15, 1-norm(ct)^2], verbose=false, stopIfNegative=true)
        (val, _) = dbisection(polPos, dpolPos, ddpolPos, interval=[α+1e-15, 1-norm(ct)^2], verbose=false, stopIfNegative=true)

        return val<=0
    end
end


function Base.in(x::AbstractVecOrMat, elli::Ellipsoid)
    return (x-elli.c)'elli.P*(x-elli.c) ≤ 1
end


function centerDistance(elli1::Ellipsoid,elli2::Ellipsoid)
    return norm(get_center(elli1)-get_center(elli2))
end

function pointCenterDistance(elli::Ellipsoid, x)
    return norm(get_center(elli)-x)
end

function volume(elli::Ellipsoid)
    N = size(elli.P,1)
    return pi^(N/2)/(gamma(N/2+1))*det(elli.P)^(-1/2)
end

function get_center(elli::Ellipsoid)
    return elli.c
end

function get_dims(elli::Ellipsoid)
    return length(elli.c)
end

function scale(elli::Ellipsoid, α)
    return Ellipsoid(elli.P*(1/α),elli.c*α)
end


function bisection(f::Function; interval=[0, 1],  δ=1e-8, verbose=false, stopIfNegative=false)

    x = zeros(4)
    fval = zeros(4)
    x[1] = interval[1];
    x[4] = interval[2];
    
    fval[1] = f(x[1]);
    fval[4] = f(x[4]);
     
    x[2] = x[4] - (x[4]-x[1])/3;
    x[3] = x[1] + (x[4]-x[1])/3;
    fval[2] = f(x[2]);
    fval[3] = f(x[3]);
    k=0;
    while abs(x[1]-x[4])>δ && (!stopIfNegative || all(i -> i >= 0, fval))
        k=k+1;
        
        if(fval[2] <fval[3])
            x[4] = x[3];
            fval[4] = fval[3];

            x[3] = x[2];
            fval[3] = fval[2];

            x[2] = x[4] - (x[4]-x[1])/phi;
            fval[2] =f(x[2]);
        else
            x[1] = x[2];
            fval[1] = fval[2];

            x[2] = x[3];
            fval[2] = fval[3];

            x[3] = x[1] + (x[4]-x[1])/phi;
            fval[3] = f(x[3]);
        end
        if (verbose)
            println(min(fval...))
        end
    end
    (_,im) = findmin(fval);
    (fval[im], x[im])
end



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

    β = (x[1] + x[2])/2
    fβ = f(β);
    dfβ = df(β);
    k=0;
    if (verbose)
        println(string(dfval[1])*" "*string(dfval[2])*"    "*string(fβ))
    end
    while abs(x[1]-x[2])>δ && (!stopIfNegative || (fβ>0 && 2*fβ<-L*(x[2]-x[1])^2))
        
        if (verbose)
            println(string(dfval[1])*" "*string(dfval[2])*"    "*string(fβ))
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
        fβ = f(β);
        dfβ = df(β);
    end
    (fβ, β)
end

export Ellipsoid

end # module
