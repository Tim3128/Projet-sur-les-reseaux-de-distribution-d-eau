function [F,G,H,ind]=OracleDH(lambda,ind)

    z = Ar'*pr + Ad'*lambda;
    q = (abs(z)./r).^(1/2).*sign(-z);

    if ind == 2 then
        F = -((1/3)*q'*(r.*q.*abs(q)) + pr'*Ar*q + lambda'*(Ad*q-fd));
        G = 0;
        H = 0;
    end
    if ind == 3 then
        F = 0;
        G = -(Ad*q - fd);
        H = 0;
    end
    if ind == 4 then
        F = -((1/3)*q'*(r.*q.*abs(q)) + pr'*Ar*q + lambda'*(Ad*q-fd));
        G = -(Ad*q - fd);
        H = 0;
    end
    if ind == 5 then
        [n,k] = size(z);
        F = 0;
        G = 0;
        H = Ad*diag( ones(n)./( 2*sqrt(r).*sqrt(abs(z)) ) )*Ad';
    end
    if ind == 6 then
        F = 0;
        G = -(Ad*q - fd);
        H = Ad*diag( ones(n)./( 2*sqrt(r).*sqrt(abs(z)) ) )*Ad';
    end
    if ind == 7 then
        [n,k] = size(z);
        F = -((1/3)*q'*(r.*q.*abs(q)) + pr'*Ar*q + lambda'*(Ad*q-fd));
        G = -(Ad*q - fd);
        H = Ad*diag( ones(n)./( 2*sqrt(r).*sqrt(abs(z)) ) )*Ad';
    end
    
endfunction
