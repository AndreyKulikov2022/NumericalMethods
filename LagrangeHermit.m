Export=true;
Lagrange=true;
Hermit=true;

n=10000;
x=linspace(0,1,n);

% F_check= @(x) (x^5+4*x^2-1);
% D_check= @(x) (5*x^4+8*x);
% N_check=6;
% 
% P=zeros(1,n);
% for i=1:n
%     P(i)=F_check(x(i));
% end
% %Lagr
% X=gen_eq(N_check-1);
% Y=zeros(1,N_check);
% for i=1:N_check
%     Y(i)=F_check(X(i));
% end
% 
% Interpolant = zeros(1,n);
% for i=1:n
%     Interpolant(i)=Lagr(X,Y,x(i));
% end
% error=max(abs(Interpolant-P));
% disp(error);
% 
% %Herm
% X=gen_eq(N_check/2-1);
% Y=zeros(1,N_check/2);
% dY=zeros(1,N_check/2);
% for i=1:N_check/2
%     Y(i)=F_check(X(i));
%     dY(i)=D_check(X(i));
% end
% 
% Interpolant = zeros(1,n);
% for i=1:n
%     Interpolant(i)=Herm(X,Y,dY,x(i));
% end
% error=max(abs(Interpolant-P));
% disp(error);
% 




F1=@(x) 1/(1+x^2);
F2=@(x)(x-1/2)^2*sign(x-1/2);
F3=@(x)abs(x-1/2);
F4=@(x)sqrt(1-x^2);

D1=@(x) -2*x/((1+x^2)^2);
D2=@(x) 2*(x-1/2)*sign(x-1/2);
D3=@(x) sign(x-1/2);
D4=@(x) -x/sqrt(1-x^2);

Fs={F1,F2,F3,F4};
Ds={D1,D2,D3,D4};

for i_f=1:length(Fs)
    F=zeros(1,n);
    for i=1:n
        F(i)=Fs{i_f}(x(i));
    end
    % Lagrange
    if Lagrange
    N=[10,20,40,80];
    err_eq=zeros(1,length(N));
    err_cheb=zeros(1,length(N));
    err_asin=zeros(1,length(N));
    for iter=1:length(N)
        X=gen_eq(N(iter));
        Y=zeros(1,N(iter)+1);
        for i=1:N(iter)+1
            Y(i)=Fs{i_f}(X(i));
        end

        Interpolant = zeros(1,n);
        for i=1:n
            Interpolant(i)=Lagr(X,Y,x(i));
        end
        error=max(abs(Interpolant-F));
        err_eq(iter)=error;
        plot_origin(x,F);
        %plot(X,Y,'o');
        plot(x,Interpolant,'r');
        text(0.5,0.03,['Equispaced, ','$N= $',num2str(N(iter)),', ','$\|f(x)-y(x)\|_{\infty}= $',num2str(error,6),'.'],...
            'Units','normalized','HorizontalAlignment','center','Interpreter','latex','FontSize',8);
        ax = gca;
        if Export exportgraphics(ax,['LagrangeHermit/','Lagr_Equi','_F',num2str(i_f),'_N',num2str(N(iter)),'.png']); end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X=gen_cheb(N(iter));
        Y=zeros(1,N(iter)+1);
        for i=1:N(iter)+1
            Y(i)=Fs{i_f}(X(i));
        end
        Interpolant = zeros(1,n);
        for i=1:n
            Interpolant(i)=Lagr(X,Y,x(i));
        end
        error=max(abs(Interpolant-F));
        err_cheb(iter)=error;
        plot_origin(x,F);
        %plot(X,Y,'o');
        plot(x,Interpolant,'r');
        text(0.5,0.03,['Chebyshev, ','$N= $',num2str(N(iter)),', ','$\|f(x)-y(x)\|_{\infty}= $',num2str(error,6),'.'],...
            'Units','normalized','HorizontalAlignment','center','Interpreter','latex','FontSize',8);
        ax = gca;
        if Export exportgraphics(ax,['LagrangeHermit/','Lagr_Cheb','_F',num2str(i_f),'_N',num2str(N(iter)),'.png']); end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X=gen_asin(N(iter));
        Y=zeros(1,N(iter)+1);
        for i=1:N(iter)+1
            Y(i)=Fs{i_f}(X(i));
        end
        Interpolant = zeros(1,n);
        for i=1:n
            Interpolant(i)=Lagr(X,Y,x(i));
        end
        error=max(abs(Interpolant-F));
        err_asin(iter)=error;
        plot_origin(x,F);
        %plot(X,Y,'o');
        plot(x,Interpolant,'r');
        text(0.5,0.03,['Asin, ','$N= $',num2str(N(iter)),', ','$\|f(x)-y(x)\|_{\infty}= $',num2str(error,6),'.'],...
            'Units','normalized','HorizontalAlignment','center','Interpreter','latex','FontSize',8);
        ax = gca;
        if Export exportgraphics(ax,['LagrangeHermit/','Lagr_Asin','_F',num2str(i_f),'_N',num2str(N(iter)),'.png']); end
    end
   
    plot_error(N,err_eq);
    ax = gca;
    if Export exportgraphics(ax,['LagrangeHermit/','Lagr_eq_err','_F',num2str(i_f),'.png']); end
    
    plot_error(N,err_cheb);
    ax = gca;
    if Export exportgraphics(ax,['LagrangeHermit/','Lagr_cheb_err','_F',num2str(i_f),'.png']); end
    
    plot_error(N,err_asin);
    ax = gca;
    if Export exportgraphics(ax,['LagrangeHermit/','Lagr_asin_err','_F',num2str(i_f),'.png']); end
    end
    % Hermit
    if Hermit
    N=[5,10,20,40];
    err_eq=zeros(1,length(N));
    err_cheb=zeros(1,length(N));
    err_asin=zeros(1,length(N));
    for iter=1:length(N)
        X=gen_eq(N(iter));
        Y=zeros(1,N(iter)+1);
        dY=zeros(1,N(iter)+1);
        for i=1:N(iter)+1
            Y(i)=Fs{i_f}(X(i));
            dY(i)=Ds{i_f}(X(i));
        end

        Interpolant = zeros(1,n);
        for i=1:n
            Interpolant(i)=Herm(X,Y,dY,x(i));
        end
        error=max(abs(Interpolant-F));
        err_eq(iter)=error;
        plot_origin(x,F);
        %plot(X,Y,'o');
        plot(x,Interpolant,'r');
        text(0.5,0.03,['Equispaced, ','$N= $',num2str(N(iter)),', ','$\|f(x)-y(x)\|_{\infty}= $',num2str(error,6),'.'],...
            'Units','normalized','HorizontalAlignment','center','Interpreter','latex','FontSize',8);
        ax = gca;
        if Export exportgraphics(ax,['LagrangeHermit/','Herm_Equi','_F',num2str(i_f),'_N',num2str(N(iter)),'.png']); end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X=gen_cheb(N(iter));
        Y=zeros(1,N(iter)+1);
        dY=zeros(1,N(iter)+1);
        for i=1:N(iter)+1
            Y(i)=Fs{i_f}(X(i));
            dY(i)=Ds{i_f}(X(i));
        end
        Interpolant = zeros(1,n);
        for i=1:n
            Interpolant(i)=Herm(X,Y,dY,x(i));
        end
        error=max(abs(Interpolant-F));
        err_cheb(iter)=error;
        plot_origin(x,F);
        %plot(X,Y,'o');
        plot(x,Interpolant,'r');
        text(0.5,0.03,['Chebyshev, ','$N= $',num2str(N(iter)),', ','$\|f(x)-y(x)\|_{\infty}= $',num2str(error,6),'.'],...
            'Units','normalized','HorizontalAlignment','center','Interpreter','latex','FontSize',8);
        ax = gca;
        if Export exportgraphics(ax,['LagrangeHermit/','Herm_Cheb','_F',num2str(i_f),'_N',num2str(N(iter)),'.png']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X=gen_asin(N(iter));
        Y=zeros(1,N(iter)+1);
        dY=zeros(1,N(iter)+1);
        for i=1:N(iter)+1
            Y(i)=Fs{i_f}(X(i));
            dY(i)=Ds{i_f}(X(i));
        end
        Interpolant = zeros(1,n);
        for i=1:n
            Interpolant(i)=Herm(X,Y,dY,x(i));
        end
        error=max(abs(Interpolant-F));
        err_asin(iter)=error;
        plot_origin(x,F);
        %plot(X,Y,'o');
        plot(x,Interpolant,'r');
        text(0.5,0.03,['Asin, ','$N= $',num2str(N(iter)),', ','$\|f(x)-y(x)\|_{\infty}= $',num2str(error,6),'.'],...
            'Units','normalized','HorizontalAlignment','center','Interpreter','latex','FontSize',8);
        ax = gca;
        if Export exportgraphics(ax,['LagrangeHermit/','Herm_Asin','_F',num2str(i_f),'_N',num2str(N(iter)),'.png']); end
    end
    
    plot_error(N,err_eq);
    ax = gca;
    if Export exportgraphics(ax,['LagrangeHermit/','Herm_eq_err','_F',num2str(i_f),'.png']); end
 
    plot_error(N,err_cheb);
    ax = gca;
    if Export exportgraphics(ax,['LagrangeHermit/','Herm_cheb_err','_F',num2str(i_f),'.png']); end
    
    plot_error(N,err_asin);
    ax = gca;
    if Export exportgraphics(ax,['LagrangeHermit/','Herm_asin_err','_F',num2str(i_f),'.png']); end
    end
end

function X=gen_eq(N)
X=linspace(0,1,N+1);
end

function X=gen_cheb(N)
X=1/2 - cos(linspace(0,1,N+1)*pi)/2;
end

function X=gen_asin(N)
X=1/2 + asin(2*linspace(0,1,N+1)-1)/pi;
end

function res=Lagr(X,Y,x)
res=0;
for i=1:length(X)
    res=res+Y(i)*Li(X,x,i);
end
end

function L=Li(X,x,i)
L=1;
for j=1:length(X)
    if j~=i
        L=L*(x-X(j))/(X(i)-X(j));
    end
end
end

function dL=dLi(X,i)
dL=0;
for j=1:length(X)
    if j~=i
        dL=dL+(1/(X(i)-X(j)));
    end
end
end

function p=Herm(X,Y,dY,x)
p=0;
for i=1:length(X)
    Hi=(1-2*dLi(X,i)*(x-X(i)))*(Li(X,x,i)^2);
    Hdi=(x-X(i))*(Li(X,x,i)^2);
    p=p+Hi*Y(i)+Hdi*dY(i);
end

end


function plot_origin(x,F)
figure
hold on
plot(x,F,'b');
xlim([min(x),max(x)]);
ylim([min(F)-(max(F)-min(F))/15,max(F)]);
ylabel('$f(x)$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
end

function plot_error(N,err)
    figure
    hold on
    ord = polyfit(log(N), log(err), 1);
    plot(N,err);
    xlabel('Number of data points.');
    ylabel('Error');
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log');
end