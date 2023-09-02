input_acoustic;
for m=1:sk
ep(1)=k_ana(m)*ep(2);
%ele_coord stores the coordinates of each point on the element.    
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
%*****Solve the BEM model******
[pr,nv]=bem_solveq(G,H,bcpr,bcnv,bcim);
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
