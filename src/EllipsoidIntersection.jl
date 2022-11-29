
function Base.:∩(elli1::Ellipsoid, elli2::Ellipsoid)
    intersect(elli1, elli2)
end

function intersect(elli1::Ellipsoid, elli2::Ellipsoid)
    if  elli1.c==elli2.c
        return true
    elseif (elli1.c ∈ elli2) || (elli2.c ∈ elli1)
        return true
    else 
        L = cholesky((elli2.P+elli2.P')/2).L #TODO remove ' when Ellipsoid constructor fixed
        P = L\elli1.P/L';
        c = L'*(elli1.c -elli2.c)
        specDecomp = eigen(P)
        lb = specDecomp.values
        ct = specDecomp.vectors'*c

        polPos(β) = -(-1-β + sum((β*lb./(1 .+ β*lb)).*(ct.^2)))
        dpolPos(β) = -(-1 + sum((lb./(1 .+ β*lb).^2).*(ct.^2)))
        ddpolPos(β) = -(-2*sum((lb.^2 ./(1 .+ β*lb).^3).*(ct.^2)))
        if(norm(ct)^2-1<0)
            return false
        end
        (val, _) = dbisection(polPos, dpolPos, ddpolPos, interval=[0, norm(ct)^2-1], verbose=false) #, stopIfPositive=true)

        return val>=0
    end
end


function get_ℓ_ast_intersect(elli1::Ellipsoid, elli2::Ellipsoid)
    L = cholesky((elli2.P+elli2.P')/2).U #TODO remove ' when Ellipsoid constructor fixed
    P = L'\elli1.P/L;
    specDecomp = eigen(P)
    lb = specDecomp.values
    ct = specDecomp.vectors'*L*(elli1.c -elli2.c)

    polPos(β) = -(-1-β + sum((β*lb./(1 .+ β*lb)).*(ct.^2)))
    dpolPos(β) = -(-1 + sum((lb./(1 .+ β*lb).^2).*(ct.^2)))
    ddpolPos(β) = -(-2*sum((lb.^2 ./(1 .+ β*lb).^3).*(ct.^2)))

    ub = norm(ct)^2-1;#α*2;
    #while dpolPos(ub) < 0
    #    ub*=2
    #end
    (val, β) = dbisection(polPos, dpolPos, ddpolPos, interval=[0, ub], verbose=false, stopIfNegative=false)

    return -val+1, β
end


