clear;
%******Properties of acoustic medium******
%[angular frequency, sound velocity, density]
k_ana=2;
ep=[1 1 1];
rad=1;
ndiv=9;
%ndiv=max(8,round((4*k_ana(m)*rad)/3)+1);
phi=pi/(2*(ndiv));
theta=pi/(2*(ndiv));
sk=size(k_ana,1);
p_kvar=zeros(sk,1);
    ep(1)=k_ana*ep(2);

mu=1i*k_ana/(k_ana^2+(pi)^2);
bcele=zeros(8,2);
bcele(1,:)=[1,1+0i];
ele_coord=zeros(8,3);
ele_coord(1,:)=[rad,0,0];
count=2;
%first row
    for j=1:(ndiv+1)
        ele_coord(count,1)=rad*cos(phi);
        ele_coord(count,2)=rad*sin(phi)*cos((j-1)*theta);
        ele_coord(count,3)=rad*sin(phi)*sin((j-1)*theta);
        bcele(count,1)=count;
        bcele(count,2)=ele_coord(count,1)/rad+sqrt(1-ele_coord(count,1)^2)*0i/rad;
        count=count+1;
    end
    %rest of the rows
    for i=2:ndiv
        theta=pi/(2*(ndiv+i-1));
        for j=1:(ndiv+i)   
        ele_coord(count,1)=rad*cos(i*phi);
        ele_coord(count,2)=rad*sin(i*phi)*cos((j-1)*theta);
        ele_coord(count,3)=rad*sin(i*phi)*sin((j-1)*theta);
        bcele(count,1)=count;
        bcele(count,2)=ele_coord(count,1)/rad+sqrt(1-ele_coord(count,1)^2)*0i/rad;
        count=count+1;
        end
     end
    for i=1:size(ele_coord,1)
        if(real(bcele(i,2))<1/100000)
            bcele(i,2)=0;
        end
        for j=1:3
            if ele_coord(i,j)<rad/100000
                ele_coord(i,j)=0;
            end
        end
    end
    edof=zeros(8,4);
    ex8=zeros(8,3);
    ey8=ex8;
    ez8=ex8;
    %first row
    count=1;
   for i=1:ndiv
       edof(count,1)=count;
       edof(count,2)=1;
       edof(count,3)=i+1;
       edof(count,4)=i+2;
       ex8(count,1)=ele_coord(1,1);
       ex8(count,2)=ele_coord(i+1,1);
       ex8(count,3)=ele_coord(i+2,1);
       ey8(count,1)=ele_coord(1,2);
       ey8(count,2)=ele_coord(i+1,2);
       ey8(count,3)=ele_coord(i+2,2);
       ez8(count,1)=ele_coord(1,3);
       ez8(count,2)=ele_coord(i+1,3);
       ez8(count,3)=ele_coord(i+2,3);
       count=count+1;
   end
   nct=2;
   ncb=2+ndiv+1;
for i=2:(ndiv)
    for j=1:(2*ndiv+2*i-3)
       if rem((count-ndiv+i),2)==1
       edof(count,1)=count;
       edof(count,2)=nct;
       edof(count,3)=ncb;
       edof(count,4)=ncb+1;
       ex8(count,1)=ele_coord(nct,1);
       ex8(count,2)=ele_coord(ncb,1);
       ex8(count,3)=ele_coord(ncb+1,1);
       ey8(count,1)=ele_coord(nct,2);
       ey8(count,2)=ele_coord(ncb,2);
       ey8(count,3)=ele_coord(ncb+1,2);
       ez8(count,1)=ele_coord(nct,3);
       ez8(count,2)=ele_coord(ncb,3);
       ez8(count,3)=ele_coord(ncb+1,3);
       count=count+1;
       else
       edof(count,1)=count;
       edof(count,2)=nct;
       edof(count,3)=ncb+1;
       edof(count,4)=nct+1;
       ex8(count,1)=ele_coord(nct,1);
       ex8(count,2)=ele_coord(ncb+1,1);
       ex8(count,3)=ele_coord(nct+1,1);
       ey8(count,1)=ele_coord(nct,2);
       ey8(count,2)=ele_coord(ncb+1,2);
       ey8(count,3)=ele_coord(nct+1,2);
       ez8(count,1)=ele_coord(nct,3);
       ez8(count,2)=ele_coord(ncb+1,3);
       ez8(count,3)=ele_coord(nct+1,3);
       count=count+1;
       ncb=ncb+1;
       nct=nct+1;
       end
    end
    ncb=ncb+2;
    nct=nct+1;
end
ez4=[ez8;ez8]; ez2=[ez4;ez4]; ez=[ez2;-ez2];
ex4=[ex8;ex8]; ex2=[ex4;-ex4]; ex=[ex2;ex2];
ey4=[ey8;-ey8]; ey2=[ey4;ey4]; ey=[ey2;ey2];
%*******Node Coordinates*******
n_ele=size(edof,1);
node_coord=zeros(n_ele,3);
bcnv=zeros(n_ele,2);
for i=1:n_ele
    id1=edof(i,2);
    id2=edof(i,3);
    id3=edof(i,4);
    node_coord(i,:)=(ele_coord(id1,:)+ele_coord(id2,:)+ele_coord(id3,:))/3;
    bcnv(i,1)=i;
    %bcnv(i,2)=(bcele(id1,2)+bcele(id2,2)+bcele(id3,2))/3;
    bcnv(i,2)=node_coord(i,3)/(norm(node_coord(i,:)))-sqrt(1-(node_coord(i,1)/norm(node_coord(i,:)))^2)*0i/rad;
    %bcnv(i,2)=node_coord(i,1)/cos(phi)+sqrt(100-(node_coord(i,1)))^2*0i/rad;
end
n_node=n_ele;
%*****Reversed element normal direction******
n=ones(8*n_ele,1); n((n_ele+1):3*n_ele)=-1; n((4*n_ele+1):(5*n_ele))=-1; n((7*n_ele+1):(8*n_ele))=-1;
%******Boundary condition matrices******
bcpr=[]; bcim=[];bcpi=[];bcvi=[];
oct_sign=[1,1,1,1,-1,-1,-1,-1];