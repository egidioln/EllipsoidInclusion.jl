
using EllipsoidInclusion
using Mosek
using MosekTools
using SDPA
using JuMP
using LinearAlgebra
using Printf
using DataFrames
using CSV
using Random
Random.seed!(1)


opt_sdp_Mosek = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
opt_sdp_SDPA = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)

function sdpApproach(elli1, elli2, optimizer)
    P = elli1.P
    P0 = elli2.P
    c = elli1.c
    c0 = elli2.c

    model = Model(optimizer)
    @variable(model, β >= 0)

    t(x) = transpose(x);

    @constraint(model,
        [ P0         -P0*c0        
         t(-P0*c0)  c0'*P0*c0-1 ] <= 
     β.*[ P          -P*c        
         t(-P*c)    c'*P*c-1 ], PSDCone())

    optimize!(model)
    st = solution_summary(model).termination_status
    if st == MOI.OPTIMAL
        return 1
    elseif st == MOI.INFEASIBLE || st == MOI.INFEASIBLE_OR_UNBOUNDED || dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        return 0
    elseif st == MOI.SLOW_PROGRESS
        return 2
    elseif st == MOI.NUMERICAL_ERROR
        return 3
    else
        return -1
    end
end

function sdpApproachSDPA(elli1, elli2)
    sdpApproach(elli1, elli2, opt_sdp_SDPA)
end

function sdpApproachMosek(elli1, elli2)
    sdpApproach(elli1, elli2, opt_sdp_Mosek)
end

function ourApproach(elli1,elli2)
    elli1 ∈ elli2
end


ensureInside = false

ElSpan = Dict()
El0Span = Dict()
nSpan = [3, 10, 30,100] #, 100
K=100

df = DataFrame(n=[], our_time=[], sdpa_time=[], mosek_time=[],  
our_mem=[], sdpa_mem=[], mosek_mem=[],
our_res=[], sdpa_res=[], mosek_res=[])

for n = nSpan
    @printf "> Test for n = %d\n" n
    ElSpan[n] = []
    El0Span[n] = []
    for k=1:(K+1)
        aux = randn(n,n)
        P = aux'aux

        aux = randn(n,n)
        P0 = aux'aux
        l0 = eigmin(P)
        P0 = P0./eigmax(P0)*l0
        c0 = randn(n)
        
        
        El0 = Ellipsoid(P0, c0)
        
        c = randn(n)*0.1 + c0
        while c ∉ El0
            c = randn(n)*0.1 + c0
        end
        if ensureInside
            while Ellipsoid(P, c) ∉ El0 # ensure not in
                P = P*1.1             
            end
        else
            v = randn(n)*0.1
            while Ellipsoid(P, c) ∈ El0 # ensure not in
                c+= v
            end
        end
        El = Ellipsoid(P, c)
        push!(ElSpan[n],El)
        push!(El0Span[n],El0)
        
    end

    for k=1:(K+1)
        @printf "case %d / %d\n" k (K+1)
        tOur = @timed resOur = ourApproach(ElSpan[n][k], El0Span[n][k])
        tSDPA = @timed resSDPA = sdpApproachSDPA(ElSpan[n][k], El0Span[n][k])
        tMosek = @timed resMosek = sdpApproachMosek(ElSpan[n][k], El0Span[n][k])
        if k>1
            push!(df, [n, tOur.time, tSDPA.time,  tMosek.time,
             tOur.bytes, tSDPA.bytes,  tMosek.bytes,  
             resOur, resSDPA, resMosek])
        end
    end
end
filename = ensureInside ? "data_in.csv" :  "data_out.csv"
CSV.write(filename, df)

