%Kepler solver function
function Q = kepler(X1, t, myu)
    %Unpack X1
    r1Vec = zeros(3,1);
    r1Vec(:) = X1(1:3);
    v1Vec = zeros(3,1);
    v1Vec(:) = X1(4:6);
    r1 = norm(r1Vec);
    v1 = norm(v1Vec);

    eps = 1e-13; %Epsillon value to check that delta x is sufficiently small
    eps1 = 0.01; %Epsillon value for C and S functions

    oneOver_a = 2/r1 - v1^2/myu;
    a = oneOver_a^(-1);
    hVec = cross(r1Vec,v1Vec);
    e_vec = cross(v1Vec, hVec)/myu - r1Vec/r1;
    e = norm(e_vec);
    
    if e > 1
        x_g = sign(t)*sqrt(-a)*log((-2*myu*t)/(a*(dot(r1Vec,v1Vec)+sign(t)*sqrt(-myu*a)*(1-r1/a))));
    elseif e == 1
        x_g = 0;
    else
        n = sqrt(myu/a^3);
        x_g = sqrt(a)*n*t;
    end

    x = x_g;
    check = true;
    i = 1;
    Q.x(i) = x;

    while check
        Q.iteration(i) = i;
        z = x^2/a;
        S = Sfunc(z,eps1);
        C = Cfunc(z,eps1);
        t_g = (x^3*S + dot(r1Vec, v1Vec)/sqrt(myu)*x^2*C + r1*(1 - z*S)*x)/sqrt(myu);
        r = x^2*C + dot(r1Vec, v1Vec)/sqrt(myu)*x*(1 - z*S) + r1*(1 - z*C);
        q = (t - t_g);
        Q.q(i) = q;
        dx = sqrt(myu)*q/r;
        x = x + dx;
        if abs(dx) < eps
            check = false;
        end
        i = i + 1;
        Q.x(i) = x;
    end

    z = x^2/a;
    Q.z = z;
    S = Sfunc(z, eps1);
    Q.S = S;
    C = Cfunc(z,eps1);
    Q.C = C;

    f = 1 - x^2/r1*C;
    Q.f = f;
    g = t - x^3/sqrt(myu)*S;
    Q.g = g;

    r2Vec = f*r1Vec + g*v1Vec;
    Q.r2Vec = r2Vec;
    r2 = norm(r2Vec);

    f_dot = sqrt(myu)/(r1*r2)*x*(z*S - 1);
    Q.f_dot = f_dot;
    g_dot = 1 - x^2/r2*C;
    Q.g_dot = g_dot;

    v2Vec = f_dot*r1Vec + g_dot*v1Vec;
    Q.v2Vec = v2Vec;

end