%Derivative of C(z)
function dC = dCfunc(z, eps)
    if z >= eps || z <= -eps
        dC = (1 - 2*Cfunc(z, eps) - z*Sfunc(z, eps))/(2*z);
    else
        dC = -1/factorial(4) + 2*z/factorial(6) - 3*z^2/factorial(8) + 4*z^3/factorial(10) - 5*z^4/factorial(12) + 6*z^5/factorial(14);
    end
end