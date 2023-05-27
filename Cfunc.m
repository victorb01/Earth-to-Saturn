%Computes C(z)
function C = Cfunc(z, eps)
    if z >= eps
        C = (1 - cos(sqrt(z)))/z;
    elseif z <= -eps
        C = (1 - cosh(sqrt(-z)))/z;
    else
        C = 1/2 - z/factorial(4) + z^2/factorial(6) - z^3/factorial(8) + z^4/factorial(10) - z^5/factorial(12);
    end
end