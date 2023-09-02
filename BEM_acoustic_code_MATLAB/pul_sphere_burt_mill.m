%******Properties of acoustic medium******
% wavenumber k_ana
% [angular frequency, sound velocity, density] ep
% radius rad
% number of divisions along generator arc of sphere for an octant ndiv
k_ana=0.2;
ep=[k_ana*341 341 1.2];
rad=10;
ndiv=22;
%ndiv=max(8,round((4*k_ana*rad)/3)+1);
phi=pi/(2*(ndiv));
theta=pi/(2*(ndiv));
sk=size(k_ana,1);
p_kvar=zeros(sk,1);
for m=1:sk
    lambda=2*pi/(k_ana(m)*rad);
    ep(1)=k_ana(m)*ep(2);
    mu=1i/(1i*(k_ana*rad)+lambda);
%ele_coord stores the coordinates of each point on the element.    
ele_coord=zeros(8,3);
ele_coord(1,:)=[rad,0,0];
count=2;
%first row of sphere
    for j=1:(ndiv+1)
        ele_coord(count,1)=rad*cos(phi);
        ele_coord(count,2)=rad*sin(phi)*cos((j-1)*theta);
        ele_coord(count,3)=rad*sin(phi)*sin((j-1)*theta);
        count=count+1;
    end
    %rest of the rows
    for i=2:ndiv
        theta=pi/(2*(ndiv+i-1));
        for j=1:(ndiv+i)   
        ele_coord(count,1)=rad*cos(i*phi);
        ele_coord(count,2)=rad*sin(i*phi)*cos((j-1)*theta);
        ele_coord(count,3)=rad*sin(i*phi)*sin((j-1)*theta);
        count=count+1;
        end
     end
    for i=1:size(ele_coord,1)
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
%mirroring the coordinates to all the other octants
ez4=[ez8;ez8]; ez2=[ez4;ez4]; ez=[ez2;-ez2];
ex4=[ex8;ex8]; ex2=[ex4;-ex4]; ex=[ex2;ex2];
ey4=[ey8;-ey8]; ey2=[ey4;ey4]; ey=[ey2;ey2];
%*******Node Coordinates*******
n_ele=size(edof,1);
node_coord=zeros(n_ele,3);
for i=1:n_ele
    id1=edof(i,2);
    id2=edof(i,3);
    id3=edof(i,4);
    node_coord(i,:)=(ele_coord(id1,:)+ele_coord(id2,:)+ele_coord(id3,:))/3;
end
n_node=n_ele;
%*****Reversed element normal direction to take care of symmetry******
n=ones(8*n_ele,1); n((n_ele+1):3*n_ele)=-1; n((4*n_ele+1):(5*n_ele))=-1; n((7*n_ele+1):(8*n_ele))=-1;
%******Assemble influence matrices******
H=zeros(n_node,n_node); G=zeros(n_node,n_node); Mt=G; N=H;
%Quadrature rules
%The number of gauss points is gloop
%XQON,WQON when the node lies on element
%XQOFF, WQOFF when the node does not lie on element
%Alternatively, N_count can be initialized separately for the two cases
N_count=4;
[XQOFF,WQOFF]=gauss2dgen(N_count,-1,1);
gloop=(N_count*(N_count+1)/2)-1;
XQON=zeros(3*gloop,2);
WQON=zeros(3*gloop,1);
N_gauss=3*gloop;
for j=1:gloop
        XQON(j,1)=XQOFF(j,1)/3+XQOFF(j,2);
        XQON(j,2)=XQOFF(j,1)/3;
        WQON(j)=WQOFF(j)/3;
        XQON(j+gloop,1)=XQOFF(j,1)/3;
        XQON(j+gloop,2)=XQOFF(j,1)/3+XQOFF(j,2);
        WQON(j+gloop)=WQOFF(j)/3;
        XQON(j+2*gloop,1)=(1/3)*(1+2*XQOFF(j,1)-XQOFF(j,2));
        XQON(j+2*gloop,2)=(1/3)*(1-XQOFF(j,1)+2*XQOFF(j,2));
        WQON(j+2*gloop)=WQOFF(j)/3;
end


%k will traverse therough the node coordinate array
for k=1:(n_node)
r=1;
Ler=zeros(1,n_ele);
Mer=zeros(1,n_ele);
Mter=zeros(1,n_ele);
Ner=zeros(1,n_ele);
%j will traverse the element coordinate array. We access all the elements
%for the entire sphere, but the nodes are accessed only for the first
%octant.
for j=1:(8*n_ele)
if r==(n_ele+1)
r=1;
end
if k==r
gauss=XQON;
weight=WQON;
g_loop=N_gauss;
else
gauss=XQOFF;
weight=WQOFF;
g_loop=gloop;
end
%calculate the coefficients for velocity and pressure matrix
[Lk,Mk,Mkt,Nk]=bem_infl3q_trial(node_coord(k,:),ex(j,:),ey(j,:),ez(j,:),ep,n(j),edof,k,ele_coord,mu,g_loop,gauss,weight);
%We calculate a row first and then assemble it in the main coefficient matrix.
%This way, parallel for loop can be used to reduce time, in systems with
%more processors
Ler(r)=Ler(r)+Lk;
Mer(r)=Mer(r)+Mk;
Mter(r)=Mter(r)+Mkt;
Ner(r)=Ner(r)+Nk;
%For oscillating sphere:
% for j=1:(4*n_ele)
% add the coefficient rows
% for j=(4*n_ele+1:8*n_ele)
% subtract the coefficient rows
r=r+1;
end
%Coefficient matrices from conventional Helmhotz BEM
H(k,:)=H(k,:)+Mer;
G(k,:)=G(k,:)+Ler;
%Coefficient matrices from Derievative of conventinal Helmholtz BEM
Mt(k,:)=Mt(k,:)+Mter;
N(k,:)=N(k,:)+Ner;
end
cdia=(0.5).*ones(n_node,1);
C=diag(cdia);
%The burton miller formula
H(1:n_node,1:n_node)=H(1:n_node,1:n_node)-1.*C;
H=H+mu*N;
G=G+mu*(Mt+C);
G=-1i*ep(1)*ep(3)*G;
cond(H)
end
%******Boundary condition matrices******
bcnv=zeros(n_node,2);
for t=1:(n_node)
bcnv(t,1)=t;
bcnv(t,2)=1;
end

bcpr=[];%if pressure bc was given
bcim=[];%if impedance bc was given
%*****Solve the BEM model******
lvec=H*bcnv(:,2);
rvec=G*bcnv(:,2);
pr=H\rvec;
nv=bcnv(:,2);

%*****Post-processing the BEM model******
r_ana=[20 50];
 s_ana=size(r_ana,2);
 coordans=zeros(s_ana,3);
 for i=1:s_ana
     coordans(i,2)=r_ana(i);
 end
 p_sol=zeros(s_ana,1);
 for j=1:s_ana
p_sol(j)=bem_acouspost_burt_mill(node_coord,coordans(j,:),ex,ey,ez,ep,pr,nv,edof,n,ele_coord,mu);
end
