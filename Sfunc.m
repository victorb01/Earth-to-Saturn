%Computes S(z)
function S = Sfunc(z, eps)
    if z >= eps
        S = (sqrt(z) - sin(sqrt(z)))/z^(3/2);
    elseif z <= -eps
        S = (sinh(sqrt(-z)) - sqrt(-z))/(-z)^(3/2);
    else
        S = 1/factorial(3) - z/factorial(5) + z^2/factorial(7) - z^3/factorial(9) + z^4/factorial(11) - z^5/factorial(13);
    end
end