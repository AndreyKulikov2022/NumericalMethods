Export = false;

I1_analitical=1/101;
I2_analitical=pi/4;
I3_analitical=1/3;
I4_analitical=sqrt(pi)/sqrt(sqrt(exp(1)));

N=[2,3,4,5];

err_1=zeros(length(N),1);
err_2=zeros(length(N),1);
err_3=zeros(length(N),1);
err_4=zeros(length(N),1);

%Semi-Infinite
x_l={[0.58578, 3.41421],[0.415774556783 2.294280360279 6.289945082937],...
    [0.322547689619 1.745761101158 4.536620296921 9.395070912301],...
    [0.263560319718 1.413403059107 3.596425771041 7.085810005859 12.640800844276]};
w_l={[0.853553, 0.146446],[0.711093009929 0.278517733569 0.0103892565016],...
    [0.603154104342 0.357418692438 0.0388879085150 0.000539294705561], ...
    [0.521755610583 0.398666811083 0.0759424496817 0.00361175867992 0.0000233699723858]};

for i=1:length(N)

    I1=0;
    I2=0;

    for j=1:length(x_l{i})

        I1=I1+w_l{i}(j)*f1(x_l{i}(j));
        I2=I2+w_l{i}(j)*f2(x_l{i}(j));

    end
    err_1(i)=abs(I1-I1_analitical);
    err_2(i)=abs(I2-I2_analitical);


end


%Infinite

x_h={[-0.707106781186548, 0.707106781186548],[-1.224744871391589, 0, 1.224744871391589],...
    [-1.65060123885785, -0.524647623275290 ,1.65060123885785, 0.524647623275290], ...
    [-2.02018270456086, -0.958572464613819,0, 2.02018270456086, 0.958572464613819]};
w_h={[0.8862269254528, 0.8862269254528],[0.2954089751509,1.181635900604,0.2954089751509],...
    [0.08131283544725,0.8049140900055,0.08131283544725,0.8049140900055], ...
    [0.01995324205905,0.3936193231522,0.9453087204829,0.01995324205905,0.3936193231522]};

for i=1:length(N)

    I3=0;
    I4=0;

    for j=1:length(x_h{i})

        I3=I3+w_h{i}(j)*f3(x_h{i}(j));
        I4=I4+w_h{i}(j)*f4(x_h{i}(j));

    end
    err_3(i)=abs(I3-I3_analitical);
    err_4(i)=abs(I4-I4_analitical);


end

plot_error(N,err_1);
ax = gca;
if Export exportgraphics(ax,['Integrals/','Improper1','.png']); end

plot_error(N,err_2);
ax = gca;
if Export exportgraphics(ax,['Integrals/','Improper2','.png']); end


plot_error(N,err_3);
ax = gca;
if Export exportgraphics(ax,['Integrals/','Improper3','.png']); end


plot_error(N,err_4);
ax = gca;
if Export exportgraphics(ax,['Integrals/','Improper4','.png']); end


function res=f1(x)

res=sin(x/10)/10;

end

function res=f2(x)

res=1/(1+exp(-2*x));

end

function res=f3(x)

res=abs(x)/3;

end


function res=f4(x)

res=cos(x);

end

function plot_error(N,err)
    figure
    hold on
    ord = polyfit(log(N), log(err), 1);
    plot(N,err,'b');
    plot(N,(N.^ord(1)),'r');
    text(N(2),2*N(2).^ord(1),['~$\frac{1}{N^{',num2str(-ord(1)),'}}$'],'Interpreter','latex','Color','r');
    xlabel('Number of points.');
    ylabel('error','Interpreter','latex');
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log');
end