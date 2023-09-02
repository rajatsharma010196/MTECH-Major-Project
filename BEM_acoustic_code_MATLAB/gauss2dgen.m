function [gauss,w]=gauss2dgen(N,a,b)
xi=zeros(N,1);
eta=xi;
w=xi;
[xi1,wxi1]=lgwt(N-1,a,b);
count=1;
for i=1:(N-1)
    [xi2,wxi2]=lgwt((N+1-i),a,b);
    for j=(1:N+1-i)
        xi(count)=(1+xi1(i))/2;
        eta(count)=(1-xi1(i))*(1+xi2(j))/4;
        w(count)=(1-xi1(i))*wxi1(i)*wxi2(j)/8;
        count=count+1;
    end
end
gauss=[xi eta];

