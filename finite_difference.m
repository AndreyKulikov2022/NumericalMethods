%% Pade Approximation
syms h;
syms D;
syms a;
syms b;
syms c;
syms d;
syms e;

f_2=taylor(exp(-2*h*D),D,'Order',8);
f_1=taylor(exp(-h*D),D,'Order',8);
f=1;

expr= h*D*(a*f_2+b*f_1+1)-(c*f_2 + d*f_1 + e);
coefficients=coeffs(expr,D);
eqs=coefficients(1:5)==[0,0,0,0,0]
[A,B] = equationsToMatrix(eqs, [a, b, c, d, e]);
X = linsolve(A,B)
subs(coefficients(6), [a, b, c, d, e], X')


%% Finite difference
syms h;
syms D;
syms a;
syms b;
syms c;
syms d;
syms e;

f_1=taylor(exp(-h*D),D,'Order',8);
f=1;
f1=taylor(exp(h*D),D,'Order',8);
f2=taylor(exp(2*h*D),D,'Order',8);
f3=taylor(exp(3*h*D),D,'Order',8);

p=a*f_1+b*f+c*f1+d*f2+e*f3;
coefficients=coeffs(p,D);
eqs=coefficients(1:5)==[0,0,1,0,0];
[A,B] = equationsToMatrix(eqs, [a, b, c, d, e]);
X = linsolve(A,B)
subs(coefficients(6), [a, b, c, d, e], X')