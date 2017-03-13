function [F,G,H,ind]=OraclePH(qc,ind)
    q=q0+B*qc
    F=(1/3)*q'*(r.*q.*abs(q))+pr'*Ar*q
    G=B'*(r.*q.*abs(q)+Ar'*pr)
    for j=1:n-md
        H(:,j)=B(:,j).*r.*(q0+B*qc)
    end
endfunction
