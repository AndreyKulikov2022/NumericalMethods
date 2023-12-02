Export=false;

syms x;
df = diff(1/((x-1)^2 + 0.002) + 1/((x-0.2)^2 + 0.005) - 5);

N=[8,16,32,64,128,256];
err_trap=zeros(length(N),1);
err_trap_end=zeros(length(N),1);
err_simp=zeros(length(N),1);
I_analitical=(1/sqrt(0.002)*atan(1/sqrt(0.002)))+(1/sqrt(0.005)*atan(0.8/sqrt(0.005)))+(1/sqrt(0.005)*atan(0.2/sqrt(0.005)))-5;
for i=1:length(N)
   X=linspace(0,1,N(i)+1);
   I_trap=0;
   I_simp=0;
   for j=1:length(X)-1
        I_trap=I_trap + trapesoidal(X(j),X(j+1));
   end
   I_trap_end=I_trap-(1/N(i)^2)/12*(eval(subs(df,x,1))-eval(subs(df,x,0)));
    
   for j=1:2:length(X)-2
        I_simp=I_simp + simpson(X(j),X(j+2));
   end

   err_trap(i)=abs((I_analitical-I_trap)/I_analitical);
   err_trap_end(i)=abs((I_analitical-I_trap_end)/I_analitical);
   err_simp(i)=abs((I_analitical-I_simp)/I_analitical);
end

plot_error(N,err_trap);
ax = gca;
if Export exportgraphics(ax,['Integrals/','Trapesoidal','.png']); end

plot_error(N,err_trap_end);
ax = gca;
if Export exportgraphics(ax,['Integrals/','Trapesoidal_End','.png']); end


plot_error(N,err_simp);
ax = gca;
if Export exportgraphics(ax,['Integrals/','Simpson','.png']); end

tols=[1e-1,1e-3,1e-5,1e-7];
n=1000;
err=zeros(length(tols),1);
I_analitical=(1/sqrt(0.002)*atan(1/sqrt(0.002)))+(1/sqrt(0.005)*atan(0.8/sqrt(0.005)))+(1/sqrt(0.005)*atan(0.2/sqrt(0.005)))-5;

number_of_calls=zeros(length(tols),1);


for i=1:length(tols)

    [I,distrib]=I_mid(0,1,tols(i));
    figure
    hold on
    plot(linspace(0,1,n),f(linspace(0,1,n)));
    plot(distrib,zeros(length(distrib),1),'o');
    ax=gca;
    if Export exportgraphics(ax,['Integrals/','Adaptive_',num2str(i),'.png']); end
    err(i)=abs(I_analitical-I);
    number_of_calls(i)=length(distrib);

end

plot_error(number_of_calls,tols);
ax = gca;
ylabel('tolerance','Interpreter','latex');
if Export exportgraphics(ax,['Integrals/','Adaptive','.png']); end

figure
hold on
[dens,xi] = ksdensity(distrib,'Support','positive');
histogram(distrib,100);
%plot(xi,dens);
xi=xi(xi<1.01);
fd2=zeros(length(xi),1);
for i=1:length(xi)
    fd2(i)=abs(eval(subs((diff(df)),x,xi(i))));
end

plot(xi,fd2/100);
ylabel('$|\frac{d^2 f}{dx^2}|10^{-2}$','Interpreter','latex');
xlabel('x');
ax = gca;
if Export exportgraphics(ax,['Integrals/','Adaptive_distrib','.png']); end


function I = trapesoidal(x0,x1)

I=(x1-x0)/2*(f(x0)+f(x1));

end

function I = simpson(x0,x1)

I =(x1-x0)/6 * (f(x0)+4*f((x0+x1)/2)+f(x1));

end


function res=f(x)
res = 1./((x-1).^2 + 0.002) + 1./((x-0.2).^2 + 0.005) - 5;
end 

function [I,distrib]=I_mid(a,b,tol)
I1=(b-a)*f((a+b)/2);
I2=(b-a)/2*f((a+(a+b)/2)/2)+(b-a)/2*f(((a+b)/2+b)/2);
if abs(I2-I1)<tol
    distrib=(a+b)/2;
    I = I1;
else
    [Ia,distrib_a]=I_mid(a,(a+b)/2,tol/2);
    [Ib,distrib_b]=I_mid((a+b)/2,b,tol/2);
    I=Ia+Ib;
    distrib=[distrib_a, distrib_b];
end
end

function plot_error(N,err)
    figure
    hold on
    ord = polyfit(log(N), log(err), 1);
    plot(N,err,'b');
    plot(N,(N.^ord(1)),'r');
    text(N(2),2*N(2).^ord(1),['~h^{',num2str(-ord(1)),'}'],'Color','r');
    xlabel('Number of points.');
    ylabel('$\frac{|I-I_{num}|}{|I|}$','Interpreter','latex');
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log');
end