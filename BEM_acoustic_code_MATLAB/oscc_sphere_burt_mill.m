input_acoustic;
%******Assemble influence matrices******
H=zeros(n_node,n_node); G=zeros(n_node,n_node); Mt=H; N=H;
%Quadrature rules
%For the hypersingular kernel, the standard gauss quadrature converges very
%slowly, so we need to take significantly more points.
N_count=4;
[gauss2,weight2]=gauss2dgen(N_count,-1,1);
gloop2=(N_count*(N_count+1)/2)-1;
%For the other three kernels we can use fewer gauss points
N_count=4;
[XQOFF,WQOFF]=gauss2dgen(N_count,-1,1);
gloop=(N_count*(N_count+1)/2)-1;
%Traversing from node to node
for k=1:(n_node)
r=1;
%To take for octant symmetry
for j=1:(size(oct_sign,2))
    for l=((j-1)*n_node+1):j*n_node
    if r==(n_node+1)
    r=1;
    end
    gauss=XQOFF;
    weight=WQOFF;
    g_loop=gloop;
    [Lk,Mk,Mkt,Nk]=bem_infl3q_trial(node_coord(k,:),ex(l,:),ey(l,:),ez(l,:),ep,n(l),edof,k,ele_coord,mu,g_loop,gauss,weight,gloop2,gauss2,weight2);
    H(k,r)=H(k,r)+oct_sign(j)*Mk;
    G(k,r)=G(k,r)+oct_sign(j)*Lk;
    Mt(k,r)=Mt(k,r)+oct_sign(j)*Mkt;
    N(k,r)=N(k,r)+oct_sign(j)*Nk;
    r=r+1;
    end
end
end
cdia=(0.5).*ones(n_node,1);
C=diag(cdia);
H(1:n_node,1:n_node)=H(1:n_node,1:n_node)-1.*C;
H=H+mu*N;
G=G+mu*(Mt+C);
G=-1i*ep(1)*ep(3)*G;
%*****Solve the BEM model******
[pr,nv]=bem_solveq(G,H,bcpr,bcnv,bcim,bcpi,bcvi,mu);
%*****Post-processing the BEM model******
coordans=zeros(2,3);
s_ana=size(coordans,1);
coordans(1,:)=node_coord(2,:);
coordans(2,:)=[0,-50,0];
 p_sol=zeros(s_ana,1);
 for j=1:s_ana
     if (norm(coordans(j,:))==norm(node_coord(1,:)))
         out_on=2;
     else
         out_on=1;
     end
%If coordans lies on the boundary but on one of the edges or corners of the traingle then out_on=2 will no longer
%be valid. It is preferable that a point within the node_triangle is taken to calculate the surface pressure. As the element discretization will increase,
%so will the angles at the corners and edges approach pi, and out_on=2 will get closer to the answer       
p_sol(j)=bem_acouspost_burt_mill(node_coord,coordans(j,:),ex,ey,ez,ep,pr,nv,edof,n,ele_coord,mu,oct_sign,out_on);
 end