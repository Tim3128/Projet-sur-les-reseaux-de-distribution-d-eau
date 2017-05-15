function [F,G,H,ind]=OraclePH(qc,ind)

    q=q0+B*qc

    if ind == 2 then
        F = (1/3)*q'*(r.*q.*abs(q))+pr'*Ar*q;
        G = 0;
        H = 0;
    end
    if ind == 3 then
        F = 0;
        G = B'*(r.*q.*abs(q)+Ar'*pr);
        H = 0;
    end
    if ind == 4 then
        F = (1/3)*q'*(r.*q.*abs(q))+pr'*Ar*q;
        G = B'*(r.*q.*abs(q)+Ar'*pr);
        H = 0;
    end
    if ind == 5 then
        F = 0;
        G = 0;
        H = 2*B'*diag(r.*abs(q))*B;
    end
    if ind == 6 then
        F = 0;
        G = B'*(r.*q.*abs(q)+Ar'*pr);
        H = 2*B'*diag(r.*abs(q))*B;
    end
    if ind == 7 then
        F = (1/3)*q'*(r.*q.*abs(q))+pr'*Ar*q;
        G = B'*(r.*q.*abs(q)+Ar'*pr);
        H = 2*B'*diag(r.*abs(q))*B;
    end
    
endfunction




