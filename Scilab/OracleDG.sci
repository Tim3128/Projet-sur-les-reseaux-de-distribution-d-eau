function [F,G,ind]=OracleDG(lambda,ind)

    z = Ar'*pr + Ad'*lambda;
    q = sqrt(abs(z)./r).*sign(-z);
    
    if ind == 2 then
        F = -((1/3)*q'*(r.*q.*abs(q)) + pr'*Ar*q + lambda'*(Ad*q-fd));
        G = 0;
    end
    if ind == 3 then
        F = 0;
        G = -(Ad*q - fd);
    end
    if ind == 4 then
        F = -((1/3)*q'*(r.*q.*abs(q)) + pr'*Ar*q + lambda'*(Ad*q-fd));
        G = -(Ad*q - fd);
    end
    
endfunction
