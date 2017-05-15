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
        F = 0;
        G = 0;
        //H = 2*B'*diag(r.*abs(q))*B;
        H=Ad*diag(1 ./ (2*abs(z).*abs(r)))*Ad';
    end
    if ind == 6 then
        F = 0;
        G = -(Ad*q - fd);
        H = Ad*diag(1 ./ (2*abs(z).*abs(r)))*Ad';
    end
    if ind == 7 then
        F = -((1/3)*q'*(r.*q.*abs(q)) + pr'*Ar*q + lambda'*(Ad*q-fd));
        G = -(Ad*q - fd);
        H = Ad*diag(1 ./ (2*sqrt(abs(z)).*abs(r)))*Ad';
    end
    
endfunction
