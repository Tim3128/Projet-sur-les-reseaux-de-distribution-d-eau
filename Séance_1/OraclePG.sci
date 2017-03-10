function [F,G,ind]=OraclePG(qc,ind)
    q=q0+B*qc
    F=(1/3)*q'*(r.*q.*abs(q))+pr'*Ar*q
    G=B'*(r.*q.*abs(q)+Ar'*pr)
endfunction
