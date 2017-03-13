function [F,G,ind]=OraclePG(qc,ind)
    
    q=q0+B*qc
    if ind == 2 then
        F = (1/3)*q'*(r.*q.*abs(q))+pr'*Ar*q;
        G = 0;
    end
    if ind == 3 then
        F = 0;
        G = B'*(r.*q.*abs(q)+Ar'*pr);
    end
    if ind == 4 then
        F=(1/3)*q'*(r.*q.*abs(q))+pr'*Ar*q;
        G=B'*(r.*q.*abs(q)+Ar'*pr);
    end
    
endfunction

