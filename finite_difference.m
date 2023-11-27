% Padé approximation of the first derivative stencil xi−2, xi−1, xi



% the second derivative xi−1, xi, xi+1, xi+2, xi+3

syms h;
f = exp(h);
T = taylor(f, 'Order', 8)