%Lambert solver function
function Q = lambert(r1Vec, r2Vec, miu, tof, direct)
    r1 = norm(r1Vec);
    Q.r1 = r1;
    r2 = norm(r2Vec);
    Q.r2 = r2;
    h = cross(r1Vec, r2Vec);
    if dot(h,[0 0 1])
        kstar = 1;
    else
        kstar = -1;
    end
    if direct
        k = kstar;
    else
        k = -kstar;
    end
    Q.k = k;

    theta = acos(dot(r1Vec, r2Vec)/(r1*r2));
    Q.theta = theta;
    zmin = -9e9;
    zmax = (2*pi)^2;
    z = 0;
    eps = 0.01;
    eps1 = 1e-9;
    iteration = 0;

    while true
        iteration = iteration + 1;
        C = Cfunc(z, eps);
        dC = dCfunc(z, eps);
        S = Sfunc(z, eps);
        dS = dSfunc(z, eps);

        [tz,dtdz,~,~,iflag] = getTimeLambert(z, miu, k, theta, r1, r2, C, S, dC, dS);

        if iflag == -1
            deltaz = deltaz/2;
            z = z - deltaz;
            continue
        end
        deltaz = (tof - tz)/dtdz;
        ztry = z + deltaz;
        while true
            if ztry < zmin || ztry > zmax
                deltaz = deltaz/2;
                ztry = z + deltaz;
            else
                break
            end
        end
        fz = abs(ztry - z);
        if fz < eps1
            break
        end
        z = ztry;
    end
    Q.fz = fz;
    Q.C = C;
    Q.dC = dC;
    Q.S = S;
    Q.dS = dS;
    Q.iter = iteration;

    A = k*sqrt((1+cos(theta))*r1*r2);
    Q.A = A;
    y = r1 + r2 - A*(1 - z*S)/sqrt(C);
    Q.y = y;
    f = 1 - y/r1;
    Q.f = f;
    g = A*sqrt(y/miu);
    Q.g = g;
    gdot = 1 - y/r2;
    Q.gdot = gdot;
    v1Vec = (r2Vec - f*r1Vec)/g;
    Q.v1Vec = v1Vec;
    v2Vec = (gdot*r2Vec - r1Vec)/g;
    Q.v2Vec = v2Vec;
    v1 = norm(v1Vec);
    v2 = norm(v2Vec);
    energy1 = v1^2/2 - miu/r1;
    energy2 = v2^2/2 - miu/r2;
    Q.energyErr = abs(energy2 - energy1);
    Q.z = z;
    
end