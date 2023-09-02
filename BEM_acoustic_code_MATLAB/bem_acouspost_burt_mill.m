function p=bem_acouspost_burt_mill(nodes,coord,ex,ey,ez,ep,pr,nv,edof,n,ele_coord,mu,oct_sign,out_on)
%****Gauss points****
N_count=4;
[gauss,weight]=gauss2dgen(N_count,-1,1);
g_loop=(N_count*(N_count+1)/2)-1;
%Size of the matrices required
[dof,~]=size(pr);
nel=dof;
G=zeros(1,dof); H=zeros(1,dof);
for s=1:size(coord,1)
r=1;
%To take care of octants
for j=1:size(oct_sign,2)
 for l=((j-1)*nel+1):j*nel
if r==(nel+1)
r=1;
end
[Lk,Mk,~,~]=bem_infl3q_trial(coord(s,:),ex(l,:),ey(l,:),ez(l,:),ep,n(l),edof,s,ele_coord,mu,g_loop,gauss,weight,g_loop,gauss,weight);
H(s,r)=H(s,r)+oct_sign(j)*Mk;
G(s,r)=G(s,r)+oct_sign(j)*Lk;
r=r+1;
end
end
%For oscillating sphere:
% for j=1:(4*n_ele)
% add the coefficients
% for j=(4*n_ele+1:8*n_ele)
% subtract the coefficients
G=-1i*ep(1)*ep(3)*G;
p=out_on*(-G*nv+H*pr);
%if the point lies on boundary use
% p=2*(-G*nv+H*pr)
end
%-----------------------------------end--------------------------------