%******Properties of acoustic medium******
%[angular frequency, sound velocity, density]
k_ana=0.1;
kz=k_ana;
ep=[1 341 1.2];
rad=10;
ndiv=8;
%ndiv=max(8,round((4*k_ana(m)*rad)/3)+1);
phi=pi/(2*(ndiv));
theta=pi/(2*(ndiv));
sk=size(k_ana,1);
%p_kvar=zeros(sk,1);
for m=1:sk
    ep(1)=k_ana(m)*ep(2);
%     mu=min(1i/(k_ana+1),(k_ana*1i/(1+k_ana)));
%     if rem(k_ana,pi)==0
%         mu=k_ana*1i/(1+k_ana);
%     end
    mu=-0i/(k_ana+1);
    %mu=-1i/(k_ana+1);
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
%node_coord=zeros(2*n_ele,3);
node_coord=[(ex4(:,1)+ex4(:,2)+ex4(:,3))./3,(ey4(:,1)+ey4(:,2)+ey4(:,3))./3,(ez4(:,1)+ez4(:,2)+ez4(:,3))./3];
bcnv=zeros(2*n_ele,2);
bcpi=zeros(2*n_ele,2);
for i=1:2*n_ele
%     id1=edof(i,2);
%     id2=edof(i,3);
%     id3=edof(i,4);
%     %node_coord(i,:)=(ele_coord(id1,:)+ele_coord(id2,:)+ele_coord(id3,:))/3;
    bcnv(i,1)=i;
    bcpi(i,1)=i;
    %bcpi(i,2)=-exp(1i*kz*norm([0,0,500]-node_coord(i,:)))/(4*pi*norm([0,0,500]-node_coord(i,:)));
    %if i<=n_ele
    bcpi(i,2)=-exp(1i*kz*(node_coord(i,2)))*rad/(norm(node_coord(i,:)));
    %else
     %   bcpi(i,2)=0;
    %end;
    %bcpi(i,2)=1;
%     rdis=[node_coord(i,1),node_coord(i,2),node_coord(i,3)];
%     ndis=node_coord(i,:)/(norm(node_coord(i,:)));
%     drndis=rdis*ndis'/(norm(rdis));
%     bcnv(i,2)=0*bcpi(i,2)*node_coord(i,3)/(norm(node_coord(i,:)))/(ep(1)*ep(3));
    
end
n_node=2*n_ele;
%*****Reversed element normal direction******
n=ones(8*n_ele,1); n((n_ele+1):3*n_ele)=-1; n((4*n_ele+1):(5*n_ele))=-1; n((7*n_ele+1):(8*n_ele))=-1;
%******Assemble influence matrices******
H=zeros(n_node,n_node); G=zeros(n_node,n_node); Mt=H; N=H;
%Quadrature rules
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

for k=1:(n_node)
r=1;
sgn=1;
for j=1:(2*n_ele)
% if r==(n_ele+1)
% r=1;
%end
% if r==n_ele+1
%     sgn=-1;
% end
if k==r
gauss=XQON;
weight=WQON;
g_loop=N_gauss;
else
gauss=XQOFF;
weight=WQOFF;
g_loop=gloop;
end
[Lk,Mk,Mkt,Nk]=bem_infl3q_trial(node_coord(k,:),ex(j,:),ey(j,:),ez(j,:),ep,n(j),[edof;edof],k,ele_coord,mu,g_loop,gauss,weight);
H(k,r)=H(k,r)+sgn*Mk;
G(k,r)=G(k,r)+sgn*Lk;
Mt(k,r)=Mt(k,r)+Mkt;
N(k,r)=N(k,r)+Nk;
% G=bem_assem_burtmill(edof,G,Lk,k,r);
% Mt=bem_assem_burtmill(edof,Mt,Mkt,k,r);
% N=bem_assem_burtmill(edof,N,Nk,k,r);
r=r+1;
end

gauss=XQOFF;
weight=WQOFF;
g_loop=gloop;
sgn=1;
for j=(2*n_ele+1):(4*n_ele)
 if r==(2*n_ele+1)
 r=1;
 end
%  if r==n_ele+1
%      sgn=-1;
%  end
[Lk,Mk,Mkt,Nk]=bem_infl3q_trial(node_coord(k,:),ex(j,:),ey(j,:),ez(j,:),ep,n(j),[edof;edof],k,ele_coord,mu,g_loop,gauss,weight);
H(k,r)=H(k,r)+sgn*Mk;
G(k,r)=G(k,r)+sgn*Lk;
Mt(k,r)=Mt(k,r)+Mkt;
N(k,r)=N(k,r)+Nk;
r=r+1;
end
sgn=1;
for j=(4*n_ele+1:8*n_ele)
if r==(2*n_ele+1)
r=1;
end

%  if r==n_ele+1
%      sgn=-1;
%  end
[Lk,Mk,Mkt,Nk]=bem_infl3q_trial(node_coord(k,:),ex(j,:),ey(j,:),ez(j,:),ep,n(j),[edof;edof],k,ele_coord,mu,g_loop,gauss,weight);
H(k,r)=H(k,r)+sgn*Mk;
G(k,r)=G(k,r)+sgn*Lk;
Mt(k,r)=Mt(k,r)+Mkt;
N(k,r)=N(k,r)+Nk;
r=r+1;
end

end
cdia=(0.5).*ones(n_node,1);
C=diag(cdia);
% for i=1:n_node
%     for j=1:i
%         H(i,j)=H(j,i);
%     end
% end
H(1:n_node,1:n_node)=H(1:n_node,1:n_node)-1.*C;
%for i=1:n_node 
%G(1:n_node,1:n_node)=1i*ep(1)*ep(3).*(G(1:n_node,1:n_node)-mu.*C);
%G(i,i)=mu*1i*ep(1)*ep(3)*C(i,i);
%end
cond(H)
%if cond(H)>60
H=H+mu*N;
G=G+mu*(Mt+C);
end
G=-1i*ep(1)*ep(3)*G;
%end
%******Boundary condition matrices******
%pr=H\(G*bcnv(:,2));
%nv=bcnv(:,2);
%p_kvar(m)=pr(1);
% for i=(n_node-2*ndiv):n_node
%     pr(i,1)=0;
% end
%bcnv(n_node+1,1)=n_node+1;
bcpr=[]; bcim=[];
%*****Solve the BEM model******
%[pr,nv]=bem_solveq(G,H,bcpr,bcnv,bcim);
% H_new=H'*H;
% G_new=H'*G;
lvec=H*bcnv(:,2);
rvec=G*bcnv(:,2)+bcpi(:,2);
%rvec=bcpi(:,2);
pr=H\(rvec);
nv=bcnv(:,2);
p_kvar=(pr(end));
fid = fopen('scattering_pressure.txt', 'a+');
fprintf(fid, 'k = %f pr = %f + %fi \n',k_ana,real(p_kvar),imag(p_kvar));
fclose(fid);

%*****Post-processing the BEM model******
 %p_point=bem_acouspost2(coord,[0 0 50],ex,ey,ez,ep,pr,nv,edof,n);
 %p_kvar(m)=p_point;
% %  
   r_ana=-1*rad;
%  xcoord=zeros(101);
%  zcoord=zeros(101);
%  pr_ext=zeros(101);
%  for i=1:101
%      for j=1:101
%          xcoord(i,j)=1*(j-51);
%          zcoord(i,j)=1*(i-51);
%          coordans=[xcoord(i,j),zcoord(i,j),0];
%          if abs(sqrt(xcoord(i,j)^2+zcoord(i,j)^2))>(1.05*rad)
%          pr_ext(i,j)=bem_acouspost_rigid_infl3q(node_coord,coordans,ex,ey,ez,ep,pr,nv,edof,n,ele_coord,mu);
%          
%          else
%          pr_ext(i,j)=NaN+1i*NaN;
%          end
%      end
%  end
%  
 s_ana=size(r_ana,2);
 coordans=[0,-10,0];
%  for i=1:n_node
%      coordans(i,:)=(node_coord(i,:)/norm(node_coord(i,:)))*10;
%  end
  s_ana=size(coordans,1);
%  for i=1:s_ana
%      coordans(i,2)=r_ana(i);
%  end
 p_sol=zeros(s_ana,1);
 for j=1:s_ana
p_sol(j)=bem_acouspost_rigid_infl3q(node_coord,coordans(j,:),ex,ey,ez,ep,pr,nv,edof,n,ele_coord,mu);
 end
 %p_sol=1*exp(1i*kz*(r_ana))-p_sol;
 irat=real(p_sol)/real(exp(1i*kz*(r_ana)));
% fid = fopen('scat_farfield.txt', 'a+');
% fprintf(fid, 'k = %f r=50 pr = %f + %fi \n',k_ana,real(p_sol),imag(p_sol));
% fclose(fid);
%end
writematrix(node_coord,'coordinates.txt');
writematrix(pr,'psurf.txt');