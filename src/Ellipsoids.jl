module Ellipsoids
using LinearAlgebra
using SpecialFunctions
using Plots
using LazySets

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

function Base.:/(elli::Ellipsoid, r::Real)
    elli*(1/r)
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

function plot_ellipsoid!(elli::Ellipsoid; label="")
    E = LazySets.Ellipsoid(elli.c, inv(elli.P))
    plot!(E, 1e-3,label=label)
end


include("EllipsoidInclusion.jl")
include("EllipsoidIntersection.jl")

export Ellipsoid

end # module
