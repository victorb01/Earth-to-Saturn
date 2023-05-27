%Derivative of S(z)
function dS = dSfunc(z, eps)
    if z >= eps || z <= -eps
        dS = (Cfunc(z, eps) - 3*Sfunc(z, eps))/(2*z);
    else
        dS = -1/factorial(5) + 2*z/factorial(7) - 3*z^2/factorial(9) + 4*z^3/factorial(11) - 5*z^4/factorial(13) + 6*z^5/factorial(15);
    end
end