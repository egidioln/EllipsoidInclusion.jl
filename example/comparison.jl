
using FileIO, JLD2
using EllipsoidInclusion
using Mosek
using MosekTools
using SDPA
using JuMP
using BenchmarkTools
using LinearAlgebra
using Printf


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

    solution_summary(model).termination_status == MOI.OPTIMAL
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

function batch!(counter,ElSpan,El0Span, approach)
    resOur = approach(ElSpan[counter[1]], El0Span[counter[1]])
    counter[1] = counter[1]+1
    resOur
end



ensureInside = false

ElSpan = Dict()
El0Span = Dict()
results = Dict()
nSpan = [3, 10, 30, 100]
K=100
res = []
for n = nSpan
    @printf "> Test for n = %d\n" n
    ElSpan[n] = []
    El0Span[n] = []
    results[n] = Dict()
    for k=1:(K+1)
        aux = randn(n,n)
        P = aux'aux

        aux = randn(n,n)
        P0 = aux'aux
        l0 = eigmin(P)
        P0 = P0./eigmax(P0)*l0
        c0 = randn(n)
        
        
        El0 = Ellipsoid(P0, c0)
        
        c = randn(n)*0.5 + c0
        while c ∉ El0
            c = randn(n)*0.2 + c0
        end
        if ensureInside
            while Ellipsoid(P, c) ∉ El0 # ensure not in
                P = P*1.01             
            end
        else
            v = randn(n)*0.01
            while Ellipsoid(P, c) ∈ El0 # ensure not in
                c+= v
            end
        end
        El = Ellipsoid(P, c)
        push!(ElSpan[n],El)
        push!(El0Span[n],El0)
    end
    counter = [1]
    ourBenchmark = @benchmarkable batch!($counter, $(ElSpan[n]), $(El0Span[n]), $ourApproach) 
    @printf "our Started\n"
    run(ourBenchmark, samples=1)
    results[n]["our"] = run(ourBenchmark, samples=K)

    counter = [1]
    SDPABenchmark = @benchmarkable batch!($counter, $(ElSpan[n]), $(El0Span[n]), $sdpApproachSDPA) 
    @printf "SDPA Started\n"
    run(SDPABenchmark, samples=1)
    results[n]["SDPA"] = run(SDPABenchmark, samples=K)

    counter = [1]
    MosekBenchmark = @benchmarkable batch!($counter, $(ElSpan[n]), $(El0Span[n]), $sdpApproachMosek) 
    @printf "Mosek Started\n"
    run(MosekBenchmark, samples=1)
    results[n]["Mosek"] = run(MosekBenchmark, samples=K)
    

    push!(res, [n, results[n]["our"].times, results[n]["SDPA"].times, results[n]["Mosek"].times,
    results[n]["our"].memory, results[n]["SDPA"].memory, results[n]["Mosek"].memory])
end

FileIO.save(ensureInside ? "data_el_in.jld2" :  "data_el_out.jld2", "res", res)