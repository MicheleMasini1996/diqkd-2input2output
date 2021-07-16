function xlog2x(x)
    if x<0 || x>1
        NaN
    elseif x==0 || x==1
        0
    else
        -x*log2(x)
    end
end

function hshannon(prob)
    sum(xlog2x.(prob))
end

function hbin(x)
    hshannon((x,1-x))
end

function φ(x)
    hbin(1/2+x/2)    
end

function hcond(pab,pb)
     hshannon(pab)-hshannon(pb)
end