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
        for j=1:n-md
            H2(:,j)=2 * B'*(r.*abs(q).*B(:,j));
        end
        H = H2;
    end
    if ind == 6 then
        F = 0;
        G = B'*(r.*q.*abs(q)+Ar'*pr);
        for j=1:n-md
            H2(:,j)=2 * B'*(r.*abs(q).*B(:,j));
        end
        H = H2;
    end
    if ind == 7 then
        F = (1/3)*q'*(r.*q.*abs(q))+pr'*Ar*q;
        G = B'*(r.*q.*abs(q)+Ar'*pr);
        for j=1:n-md
            H2(:,j)=2 * B'*(r.*abs(q).*B(:,j));
        end
        H = H2;
    end
    
endfunction

