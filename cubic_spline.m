Export=true;
n=1000;
N=[8,12,16,20];
for i=1:length(N)
    t=linspace(0,2*pi,N(i));
    h=t(2)-t(1);
    t=[t,2*pi+t(2)];%provides two equal data points
    x=cos(t)';
    y=sqrt(2)*sin(t)';
    figure
    hold on;
    plot(x,y,'.')
    plot(cos(linspace(0,2*pi,n)),sqrt(2)*sin(linspace(0,2*pi,n)),'b');
    %A *fd2 = B* f
    A=diag(2*h/3*ones(N(i)+1,1))+diag(h/6*ones(N(i),1),1)+diag(h/6*ones(N(i),1),-1);
    B=diag(-2*ones(N(i)+1,1))+diag(ones(N(i),1),1)+diag(ones(N(i),1),-1);
    % Apply periodic conditions
    A(1,:)=0;
    A(end,:)=0;
    A(1,1)=1;
    A(1,end-1)=-1;
    A(end,2)=1;
    A(end,end)=-1;
    B(1,:)=0;
    B(end,:)=0;
    xd2=linsolve(A,B*x);
    yd2=linsolve(A,B*y);

    T=linspace(0,2*pi,n);
    X=zeros(n,1);
    Y=zeros(n,1);
    spline_n=sum(T>=t(1:end-1)',1);
    for j=1:n
        X(j)=spline_el(h,x(spline_n(j)),x(spline_n(j)+1),xd2(spline_n(j)),xd2(spline_n(j)+1),t(spline_n(j)),t(spline_n(j)+1),T(j));
        Y(j)=spline_el(h,y(spline_n(j)),y(spline_n(j)+1),yd2(spline_n(j)),yd2(spline_n(j)+1),t(spline_n(j)),t(spline_n(j)+1),T(j));
    end
    plot(X,Y,'r');
    ax=gca;
    if Export exportgraphics(ax,['CubicSpline/','cubic_spline','_N',num2str(N(i)),'.png']); end
end

function X=spline_el(h,x,x1,xd2,x1d2,t,t1,T)
X=((t1-T)^3)/(6*h)*xd2 + ((T-t)^3)/(6*h)*x1d2 + (t1-T)*(x/h - h/6*xd2)+(T-t)*(x1/h - h/6*x1d2);

end