%======================================================================
%Boundary element code to find displacements in a 2D body
%----------------------------------------------------------------------
%You can modify the Elastic constants, number of elements, and the node
%coordinates according to requirement. Constant shape function
%approximation is used. 
%======================================================================
 run('cook2.m');   
g_size=0;
 node=0;
 %number of nodes
if(p_type==0)
    node=ele;
else
    node=p_type*ele;
end    
%Initialize the coefficient matrices
A=zeros(node,ele*(p_type+1));%Equation 4.24
B=zeros(node,ele*(p_type+1));%Equation 4.25
C=zeros(node);%Equation 4.20
D=zeros(node);%Equation 4.21
E=zeros(node,ele*(p_type+1));%Equation 4.22
F=zeros(node,ele*(p_type+1));%Equation 4.23
G=zeros(node);%Equation 4.26
H=zeros(node);%Equation 4.27
%Creating dominant diagonal matrix
half=0.5*zeros(ele,1);
Major=diag(half);
c_size=size(corner_node,1);
%if nodes are at some angled corners
if(c_size~=0)
    for m=1:c_size
        index=corner_node(m,1);
        Major(index,index)=0*corner_node(m,2);
    end
end

%Gaussian Quadrature
[w,xi]=get_gauss(p_type);
g_size=size(w,1);
%for the distance parameter
rx=zeros(ele,1);
ry=zeros(ele,1);
r=zeros(ele,1);
%for normals
nx=zeros(ele,1);
ny=zeros(ele,1);
%Calculating normals for all elements
for i=1:ele
    p=i+1;
    if i==ele
        p=1;
    end
    rx(i)=X_ele(p,1)-X_ele(i,1);
    ry(i)=X_ele(p,2)-X_ele(i,2);
    r(i)=sqrt(rx(i)^2 + ry(i)^2);
    nx(i)=ry(i)/r(i);
    ny(i)=-rx(i)/r(i);
end
%Calculating the coefficient matrices
for i=1:node
for k=1:g_size
a=zeros(node,ele*(p_type+1));
b=zeros(node,ele*(p_type+1));
c=zeros(node);
d=zeros(node);
e=zeros(node,ele*(p_type+1));
f=zeros(node,ele*(p_type+1));
g=zeros(node);
h=zeros(node);
N=get_N(p_type,xi(k)); 
n_size=size(N,2);
%for i=1:node
    %Calculating distances from node i to gauss point xi(k) on element j 
    for j=1:ele
        p=j+1;
        if j==ele
            p=1;
        end
        rx(j)=((X_ele(p,1)-X_ele(j,1))*xi(k))/2+((X_ele(p,1)+X_ele(j,1))*1)/2-X_node(i,1);
        ry(j)=((X_ele(p,2)-X_ele(j,2))*xi(k))/2+((X_ele(p,2)+X_ele(j,2))*1)/2-X_node(i,2);
        r(j)=sqrt(rx(j)^2+ry(j)^2);

       %jacobian is length of element by 2
       jacob = sqrt((X_ele(p,1)-X_ele(j,1))^2 + (X_ele(p,2)-X_ele(j,2))^2)/2;
       %to accomodate for the increase in nodes per element due to p
       %increase
       for m=1:n_size
        sum_pos=max(0,(p_type-1))*(j-1)+j+m-1;
        if sum_pos>node
           sum_pos=sum_pos-node;
        end
       %Coefficients obtained by using Kelvin's solution of Navier
       %relations
       
       a(i,(p_type+1)*(j-1)+m)=C1*((3-4*nu)*log(r(j))-(rx(j)/r(j))^2);% eqn 4.24
       b(i,(p_type+1)*(j-1)+m)=C1*(-((rx(j)*ry(j))/(r(j))^2));% eqn 4.25
       e(i,(p_type+1)*(j-1)+m)=C1*(-((rx(j)*ry(j))/(r(j))^2));% eqn 4.26
       f(i,(p_type+1)*(j-1)+m)=C1*((3-4*nu)*log(r(j))-(ry(j)/r(j))^2);% eqn 4.27
       c(i,sum_pos)=C2*(1/r(j))*(((1-2*nu)+2*(rx(j)/r(j))^2)*((rx(j)*nx(j)+ry(j)*ny(j))/(r(j))));% eqn 4.20
       d(i,sum_pos)=C2*(1/r(j))*(((2*rx(j)*ry(j))/(r(j)^2))*((rx(j)*nx(j)+ry(j)*ny(j))/(r(j)))-((1-2*nu)*((rx(j)*ny(j)-ry(j)*nx(j))/(r(j)))));% eqn 4.21
       g(i,sum_pos)=C2*(1/r(j))*(((2*rx(j)*ry(j))/(r(j)^2))*((rx(j)*nx(j)+ry(j)*ny(j))/(r(j)))-((1-2*nu)*((ry(j)*nx(j)-rx(j)*ny(j))/(r(j)))));% eqn 4.22
       h(i,sum_pos)=C2*(1/r(j))*(((1-2*nu)+2*(ry(j)/r(j))^2)*((rx(j)*nx(j)+ry(j)*ny(j))/(r(j))));% eqn 4.23
       
       
       %%===========Integration Using Gaussian Quadrature==========%%
       
       %For traction the nodes for each element are divided between two
       %sharing elements. This is because traction on the same node can
       %change depending on element
       A(i,(p_type+1)*(j-1)+m)=A(i,(p_type+1)*(j-1)+m)+w(k)*jacob*N(m)*a(i,(p_type+1)*(j-1)+m);
       B(i,(p_type+1)*(j-1)+m)=B(i,(p_type+1)*(j-1)+m)+w(k)*jacob*N(m)*b(i,(p_type+1)*(j-1)+m);
       E(i,(p_type+1)*(j-1)+m)=E(i,(p_type+1)*(j-1)+m)+w(k)*N(m)*jacob*e(i,(p_type+1)*(j-1)+m);
       F(i,(p_type+1)*(j-1)+m)=F(i,(p_type+1)*(j-1)+m)+w(k)*N(m)*jacob*f(i,(p_type+1)*(j-1)+m);
       %For displacements the coefficients of node for sharing elements are
       %added. The displacement on node does not change due to element
       C(i,sum_pos)=C(i,sum_pos)+w(k)*N(m)*jacob*c(i,sum_pos);
       D(i,sum_pos)=D(i,sum_pos)+w(k)*N(m)*jacob*d(i,sum_pos);
       G(i,sum_pos)=G(i,sum_pos)+w(k)*N(m)*jacob*g(i,sum_pos);
       H(i,sum_pos)=H(i,sum_pos)+w(k)*N(m)*jacob*h(i,sum_pos);
       end
end
end
end
%%To set coefficient using corner shapes
for i=1:node
C(i,i)=C(i,i)+Major(i);
H(i,i)=H(i,i)+Major(i);
end;
%%======Coefficients of displacements and tractions====%%
U_mat=zeros(2*node);
T_mat=zeros(2*node,2*ele*(p_type+1));
U_mat(1:node,1:node)=C;
U_mat(1:node,node+1:2*node)=D;
U_mat(node+1:2*node,1:node)=G;
U_mat(node+1:2*node,node+1:2*node)=H;
T_mat(1:node,1:(ele*(p_type+1)))=A;
T_mat(1:node,(ele*(p_type+1))+1:2*(ele*(p_type+1)))=B;
T_mat(node+1:2*node,1:(ele*(p_type+1)))=E;
T_mat(node+1:2*node,(ele*(p_type+1))+1:2*(ele*(p_type+1)))=F;

%%====Using the fact that there is no net force on body in static
%%analysis. Sum along a row will be 0 and Uii=-sum(Uij)==============%%
for i=1:node
            U_mat(i,i)=0;
            U_mat(i+node,i)=0;
            U_mat(i,i+node)=0;
            U_mat(i+node,i+node)=0;
        
    end
    
for i=1:node
    for(j=1:node)
        if(i~=j)
            U_mat(i+node,i)=U_mat(i+node,i)-U_mat(i+node,j);
            U_mat(i,i)=U_mat(i,i)-U_mat(i,j);
            U_mat(i+node,i+node)=U_mat(i+node,i+node)-U_mat(i+node,j+node);
            U_mat(i,i+node)=U_mat(i,i+node)-U_mat(i,j+node);
        end
    end
    
end
%%==========Arranging in the form AX=B for solution======================%%
%LHS matrix
%Arranged as u v tx ty
size1=size(u_unknown,1);
size2=size(v_unknown,1);
size3=size(tx_unknown,1);
size4=size(ty_unknown,1);
size5=size(tx_eq_no,1);
size6=size(ty_eq_no,1);
size_net=size1+size2+size3+size4;
LHS=zeros(size_net);
k=1;
%initialize LHS from unknowns
if size1~=0
for i=1:size1
    index=u_unknown(i);
    LHS(1:2*node,k)=1.*U_mat(:,index);
    k=k+1;
end
end
if size2~=0
for i=1:size2
    index=v_unknown(i);
    LHS(1:2*node,k)=1.*U_mat(:,index+node);
    k=k+1;
end
end
if size3~=0
for i=1:size3
    index=(tx_unknown(i,1)-1)*(p_type+1)+tx_unknown(i,2);
    LHS(1:2*node,k)=-1.*T_mat(:,index);
    k=k+1;
end
end
if (size4~=0)
for i=1:size4
    index=(ele*(p_type+1))+(ty_unknown(i,1)-1)*(p_type+1)+ty_unknown(i,2);
    LHS(1:2*node,k)=-1.*T_mat(:,index);
    k=k+1;
end
k;
end
%equality equations
k=k-size5-size6;
k;
if size5~=0
    for i=1:size5
        m=0;
        while m<=2
        index1=tx_eq_no(i,1+m);
        index2=tx_eq_no(i,2+m);
        for j=1:size4
            if (tx_unknown(j,1)==index1)&&(tx_unknown(j,2)==index2)
                LHS(k,size_net-size3-size4+j)=(-1)^(m/2);
                m=m+2;
            end
        end
        end
        k=k+1;
    end
end
if size6~=0
    for i=1:size6
        m=0;
        while m<=2
        index1=ty_eq_no(i,1+m);
        index2=ty_eq_no(i,2+m);
        for j=1:size4
            if (ty_unknown(j,1)==index1)&&(ty_unknown(j,2)==index2)
                LHS(k,size_net-size4+j)=(-1)^(m/2);
                m=m+2;
            end
        end
        end
        k=k+1;
    end
end
%Right Hand side
size1=size(u_known,1);
size2=size(v_known,1);
size3=size(tx_known,1);
size4=size(ty_known,1);
RHS=zeros(size_net,size1+size2+size3+size4);
k=1;
if size1~=0
for i=1:size1
    index=u_known(i,1);
    RHS(1:2*node,k)=-1.*U_mat(:,index);
    k=k+1;
end
end
if size2~=0
    for i=1:size2
        index=v_known(i,1);
        RHS(1:2*node,k)=-1.*U_mat(:,index+node);
        k=k+1;
    end
end
if size3~=0
    for i=1:size3
        index=(tx_known(i,1)-1)*(p_type+1)+tx_known(i,2);
        RHS(1:2*node,k)=1.*T_mat(:,index);
        k=k+1;
    end
end
if (size4~=0)
    for i=1:size4
    index=(ele*(p_type+1))+(ty_known(i,1)-1)*(p_type+1)+ty_known(i,2);
    RHS(1:2*node,k)=1.*T_mat(:,index);
    k=k+1;
    end
end
%Rvector
R_vector=zeros(size1+size2+size3+size4,1);
k=1;
if size1~=0
for i=1:size1
    R_vector(k)=u_known(i,2);
    k=k+1;
end
end
if size2~=0
for i=1:size2
    R_vector(k)=v_known(i,2);
    k=k+1;
end
end
if size3~=0
for i=1:size3
    R_vector(k)=tx_known(i,3);
    k=k+1;
end
end
if size4~=0
for i=1:size4
    R_vector(k)=ty_known(i,3);
    k=k+1;
end
end
Right_v=RHS*R_vector;
L_vector=LHS\Right_v;
%Preparing the displacement and traction vectors for post processing
u_node=zeros(node,1);
v_node=zeros(node,1);
tx_ele=zeros(ele*(p_type+1),1);
ty_ele=zeros(ele*(p_type+1),1);
%First collect all u, v, tx, ty from unknowns
size1=size(u_unknown,1);
size2=size(v_unknown,1);
size3=size(tx_unknown,1);
size4=size(ty_unknown,1);
k=1;
%u
if size1~=0
    for i=1:size1
        index=u_unknown(i);
        u_node(index)=L_vector(k);
        k=k+1;
    end;
end;
%v
if size2~=0
    for i=1:size2
        index=v_unknown(i);
        v_node(index)=L_vector(k);
        k=k+1;
    end;
end;
%tx
if size3~=0
    for i=1:size3
        index1=tx_unknown(i,1);
        index2=tx_unknown(i,2);
        index=(p_type+1)*(index1-1)+index2;
        tx_ele(index)=L_vector(k);
        k=k+1;
    end;
end;
%ty
if size4~=0
    for i=1:size4
        index1=ty_unknown(i,1);
        index2=ty_unknown(i,2);
        index=(p_type+1)*(index1-1)+index2;
        ty_ele(index)=L_vector(k);
        k=k+1;
    end;
end;
%Now get known u v tx ty
size1=size(u_known,1);
size2=size(v_known,1);
size3=size(tx_known,1);
size4=size(ty_known,1);
%u from known
if size1~=0
    for i=1:size1
        index=u_known(i,1);
        u_node(index)=u_known(i,2);
        
    end;
end;
%v
if size2~=0
    for i=1:size2
        index=v_known(i,1);
        v_node(index)=v_known(i,2);
        
    end;
end;
%tx
if size3~=0
    for i=1:size3
        index1=tx_known(i,1);
        index2=tx_known(i,2);
        index=(p_type+1)*(index1-1)+index2;
        tx_ele(index)=tx_known(i,3);
        
    end;
end;
%ty
if size4~=0
    for i=1:size4
        index1=ty_known(i,1);
        index2=ty_known(i,2);
        index=(p_type+1)*(index1-1)+index2;
        ty_ele(index)=ty_known(i,3);
        
    end;
end;
%final shape
X_final=X_node(:,1)+u_node;
Y_final=X_node(:,2)+v_node;
f1=figure();
hold on;
plot([X_node(:,1);X_node(1,1)],[X_node(:,2);X_node(1,2)],'-s');
plot([X_final(:);X_final(1)],[Y_final(:);Y_final(1)],'-o');
legend('Initial Shape','Final Shape','FontSize',14,'Location','southeast');
grid on;
hold off;

%%===================================================================================================%%
%%Post Processing for displacement
%%=====================================================================================================%%
if postp==1
pts=size(xcoord,1);
xcoord_matrix=zeros(dats);
ycoord_matrix=zeros(dats);
ucoord_matrix=zeros(dats);
vcoord_matrix=zeros(dats);
for i=1:dats
    for j=1:dats
        xcoord(i,1)=i/reconi;
        ycoord(i,1)=j/reconj;
    end
end
for i=1:ele
    p=i+1;
    if i==ele
        p=1;
    end
    rx(i)=X_ele(p,1)-X_ele(i,1);
    ry(i)=X_ele(p,2)-X_ele(i,2);
    r(i)=sqrt(rx(i)^2 + ry(i)^2);
    nx(i)=-ry(i)/r(i);
    ny(i)=rx(i)/r(i);
end

for i=1:dats
    for j=1:dats
    
        [in,on] = inpolygon(xcoord(i,1),ycoord(i,1),X_ele(:,1),X_ele(:,2));
        inon = in | on;                                            
                                                
        
Ap=zeros(1,ele*(p_type+1));
Bp=zeros(1,ele*(p_type+1));
Cp=zeros(1,node);
Dp=zeros(1,node);
Ep=zeros(1,ele*(p_type+1));
Fp=zeros(1,ele*(p_type+1));
Gp=zeros(1,node);
Hp=zeros(1,node);    

ucoord=zeros(1,1);
vcoord=zeros(1,1);
    for k=1:g_size
a=zeros(1,ele*(p_type+1));
b=zeros(1,ele*(p_type+1));
c=zeros(1,node);
d=zeros(1,node);
e=zeros(1,ele*(p_type+1));
f=zeros(1,ele*(p_type+1));
g=zeros(1,node);
h=zeros(1,node);
N=get_N(p_type,xi(k)); 
n_size=size(N,2);
    %Calculating distances from node i to gauss point xi(k) on element j 
    for q=1:ele
        p=q+1;
        if q==ele
            p=1;
        end
        rx(q)=((X_ele(p,1)-X_ele(q,1))*xi(k))/2+((X_ele(p,1)+X_ele(q,1))*1)/2-xcoord(i,1);
        ry(q)=((X_ele(p,2)-X_ele(q,2))*xi(k))/2+((X_ele(p,2)+X_ele(q,2))*1)/2-ycoord(i,1);
        r(q)=sqrt(rx(q)^2+ry(q)^2);

    end
    
    
  for q=1:ele
       p=q+1;
       if q==ele
            p=1;
       end
       %jacobian is length of element by 2
       jacob = sqrt((X_ele(p,1)-X_ele(q,1))^2 + (X_ele(p,2)-X_ele(q,2))^2)/2;
       %to accomodate for the increase in nodes per element due to p
       %increase
       for m=1:n_size
        sum_pos=q+m-1;
        if sum_pos>ele
           sum_pos=sum_pos-ele;
        end
       %Coefficients obtained by using Kelvin's solution of Navier
       %relations
       a(1,(p_type+1)*(q-1)+m)=C1*((3-4*nu)*log(r(q))-(rx(q)/r(q))^2);%coeff of tx in u eqn
       b(1,(p_type+1)*(q-1)+m)=C1*(-((rx(q)*ry(q))/(r(q))^2));%coeff of ty in u eqn
       e(1,(p_type+1)*(q-1)+m)=C1*(-((rx(q)*ry(q))/(r(q))^2));%coeff of u in u eqn
       f(1,(p_type+1)*(q-1)+m)=C1*((3-4*nu)*log(r(q))-(ry(q)/r(q))^2);%coeff of v in u eqn
       c(1,sum_pos)=C2*(((1-2*nu)+2*(rx(q)/r(q))^2)*((rx(q)*nx(q)+ry(q)*ny(q))/(r(q))^2));
       d(1,sum_pos)=C2*(((2*rx(q)*ry(q))/(r(q)^2))*((rx(q)*nx(q)+ry(q)*ny(q))/(r(q))^2)-((1-2*nu)*((rx(q)*ny(q)-ry(q)*nx(q))/(r(q))^2)));
       g(1,sum_pos)=C2*(((2*rx(q)*ry(q))/(r(q)^2))*((rx(q)*nx(q)+ry(q)*ny(q))/(r(q))^2)-((1-2*nu)*((ry(q)*nx(q)-rx(q)*ny(q))/(r(q))^2)));
       h(1,sum_pos)=C2*(((1-2*nu)+2*(ry(q)/r(q))^2)*((rx(q)*nx(q)+ry(q)*ny(q))/(r(q))^2));
       %for traction the nodes for each element are divided between two
       %sharing elements. This is because traction on the same node can
       %change depending on element
       Ap(1,(p_type+1)*(q-1)+m)=Ap(1,(p_type+1)*(q-1)+m)+w(k)*jacob*N(m)*a(1,(p_type+1)*(q-1)+m);
       Bp(1,(p_type+1)*(q-1)+m)=Bp(1,(p_type+1)*(q-1)+m)+w(k)*jacob*N(m)*b(1,(p_type+1)*(q-1)+m);
       %For displacements the coefficients of node for sharing elements are
       %added. The displacement on node does not change due to element
       Cp(1,sum_pos)=Cp(1,sum_pos)+w(k)*N(m)*jacob*c(1,sum_pos);
       Dp(1,sum_pos)=Dp(1,sum_pos)+w(k)*N(m)*jacob*d(1,sum_pos);
       Ep(1,(p_type+1)*(q-1)+m)=Ep(1,(p_type+1)*(q-1)+m)+w(k)*N(m)*jacob*e(1,(p_type+1)*(q-1)+m);
       Fp(1,(p_type+1)*(q-1)+m)=Fp(1,(p_type+1)*(q-1)+m)+w(k)*N(m)*jacob*f(1,(p_type+1)*(q-1)+m);
       Gp(1,sum_pos)=Gp(1,sum_pos)+w(k)*N(m)*jacob*g(1,sum_pos);
       Hp(1,sum_pos)=Hp(1,sum_pos)+w(k)*N(m)*jacob*h(1,sum_pos);
       end
end
    end
 if inon==1
ucoord_matrix(i,j)=Ap*tx_ele+Bp*ty_ele-Cp*u_node-Dp*v_node;
vcoord_matrix(i,j)=Ep*tx_ele+Fp*ty_ele-Gp*u_node-Hp*v_node;
  else
 ucoord_matrix(i,j)=NaN;
  vcoord_matrix(i,j)=NaN;
  end
    end
    
end
end;
%%===================================================================================================%%
%%Post Processing for stress
%%=====================================================================================================%%

if postp==2
pts=size(xcoord,1);
sigmaxx=zeros(pts,1);
sigmaxy=zeros(pts,1);
sigmayy=zeros(pts,1);
for i=1:ele
    p=i+1;
    if i==ele
        p=1;
    end
    rx(i)=X_ele(p,1)-X_ele(i,1);
    ry(i)=X_ele(p,2)-X_ele(i,2);
    r(i)=sqrt(rx(i)^2 + ry(i)^2);
    nx(i)=ry(i)/r(i);
    ny(i)=-rx(i)/r(i);
end

for i=1:pts
                                     
        
Dxxx=zeros(1,ele*(p_type+1));
Dyxx=zeros(1,ele*(p_type+1));
Dxxy=zeros(1,ele*(p_type+1));
Dyxy=zeros(1,ele*(p_type+1));
Dxyy=zeros(1,ele*(p_type+1));
Dyyy=zeros(1,ele*(p_type+1));
Sxxx=zeros(1,node);
Syxx=zeros(1,node);
Sxxy=zeros(1,node);
Syxy=zeros(1,node);
Sxyy=zeros(1,node);
Syyy=zeros(1,node);

ucoord=zeros(1,1);
vcoord=zeros(1,1);
    for k=1:g_size

dxxx=zeros(1,ele*(p_type+1));
dyxx=zeros(1,ele*(p_type+1));
dxxy=zeros(1,ele*(p_type+1));
dyxy=zeros(1,ele*(p_type+1));
dxyy=zeros(1,ele*(p_type+1));
dyyy=zeros(1,ele*(p_type+1));
sxxx=zeros(1,node);
syxx=zeros(1,node);
sxxy=zeros(1,node);
syxy=zeros(1,node);
sxyy=zeros(1,node);
syyy=zeros(1,node);
N=get_N(p_type,xi(k)); 
n_size=size(N,2);
dist=zeros(ele,1);
for q=1:ele
    p=q+1;
    if q==ele
        p=1;
    end
    
    if rx~=0
        ta=ry(q)/rx(q);
       dist(q)=abs(ta*xcoord(i,1)-ycoord(i,1)+X_ele(q,2)-ta*X_ele(q,1))*(rx(q)/r(q))^2;
    else
        dist(q)=abs(xcoord(i,1)-X_ele(q,1));
    end
    sig=(X_ele(q,1)-xcoord(i,1))*(X_ele(p,2)-ycoord(i,1))-(X_ele(p,1)-xcoord(i,1))*(X_ele(q,2)-ycoord(i,1));
    if sig<0
        dist(q)=-dist(q);
    end

end

    %Calculating distances from node i to gauss point xi(k) on element j 
    for q=1:ele
        p=q+1;
        if q==ele
            p=1;
        end
        rx(q)=((X_ele(p,1)-X_ele(q,1))*xi(k))/2+((X_ele(p,1)+X_ele(q,1))*1)/2-xcoord(i,1);
        ry(q)=((X_ele(p,2)-X_ele(q,2))*xi(k))/2+((X_ele(p,2)+X_ele(q,2))*1)/2-ycoord(i,1);
        r(q)=sqrt(rx(q)^2+ry(q)^2);
    end
    
    
  for q=1:ele
       p=q+1;
       if q==ele
            p=1;
       end
       %jacobian is length of element by 2
       jacob = sqrt((X_ele(p,1)-X_ele(q,1))^2 + (X_ele(p,2)-X_ele(q,2))^2)/2;
       %to accomodate for the increase in nodes per element due to p
       %increase
       for m=1:n_size
        sum_pos=q+m-1;
        if sum_pos>ele
           sum_pos=sum_pos-ele;
        end
       %Coefficients obtained by using Kelvin's solution of Navier
       %relations
       C1=1/(4*pi*(1-nu)*r(q));
       C2=GS/(2*pi*(1-nu)*r(q)^2);
       %drn=(rx(q)*nx(q)+ry(q)*ny(q))/r(q);
       %drn=dist(q)/r(q);
       drn=(nx(q)*rx(q)+ny(q)*ry(q))/r(q);
       drx=rx(q)/r(q);
       dry=ry(q)/r(q);
       
       dxxx(1,(p_type+1)*(q-1)+m)=C1*((1-2*nu)*drx+2*drx^3);%coeff of tx in u eqn
       dyxx(1,(p_type+1)*(q-1)+m)=C1*((1-2*nu)*(-dry)+2*drx*drx*dry);%coeff of ty in u eqn
       dxxy(1,(p_type+1)*(q-1)+m)=C1*((1-2*nu)*dry+2*dry*drx*drx);%coeff of u in u eqn
       dyxy(1,(p_type+1)*(q-1)+m)=C1*((1-2*nu)*(rx(q)/r(q))+2*(rx(q)/r(q))*(ry(q)/r(q))^2);%coeff of v in u eqn
       dxyy(1,(p_type+1)*(q-1)+m)=C1*((1-2*nu)*(-rx(q)/r(q))+2*(rx(q)/r(q))*(ry(q)/r(q))^2);%coeff of u in u eqn
       dyyy(1,(p_type+1)*(q-1)+m)=C1*((1-2*nu)*(ry(q)/r(q))+2*(ry(q)/r(q))^3);%coeff of v in u eqn
       sxxx(1,sum_pos)=C2*(2*drn*((1-2*nu)*drx+nu*2*drx-4*drx^3)+4*nu*nx(q)*drx^2+(1-2*nu)*(2*nx(q)*drx^2+2*nx(q))-(1-4*nu)*nx(q));
       syxx(1,sum_pos)=C2*(2*drn*((1-2*nu)*dry-4*dry*drx^2)+4*nu*nx(q)*drx*dry+(1-2*nu)*2*ny(q)*drx^2-(1-4*nu)*ny(q));
       sxxy(1,sum_pos)=C2*(2*drn*(nu*dry-4*dry*drx^2)+2*nu*(nx(q)*drx*dry+ny(q)*drx*drx)+(1-2*nu)*(2*nx(q)*drx*dry+ny(q)));
       syxy(1,sum_pos)=C2*(2*drn*(nu*drx-4*drx*dry^2)+2*nu*(ny(q)*drx*dry+nx(q)*dry*dry)+(1-2*nu)*(2*ny(q)*drx*dry+nx(q)));
       sxyy(1,sum_pos)=C2*(2*drn*((1-2*nu)*drx-4*drx*dry^2)+4*nu*ny(q)*dry*drx+(1-2*nu)*2*nx(q)*dry^2-(1-4*nu)*nx(q));
       syyy(1,sum_pos)=C2*(2*drn*((1-2*nu)*dry+nu*2*dry-4*dry^3)+4*nu*ny(q)*dry^2+(1-2*nu)*(2*ny(q)*dry^2+2*ny(q))-(1-4*nu)*ny(q));
       %for traction the nodes for each element are divided between two
       %sharing elements. This is because traction on the same node can
       %change depending on element
       Dxxx(1,(p_type+1)*(q-1)+m)=Dxxx(1,(p_type+1)*(q-1)+m)+w(k)*N(m)*jacob*dxxx(1,(p_type+1)*(q-1)+m);
       Dyxx(1,(p_type+1)*(q-1)+m)=Dyxx(1,(p_type+1)*(q-1)+m)+w(k)*N(m)*jacob*dyxx(1,(p_type+1)*(q-1)+m);
       Dxxy(1,(p_type+1)*(q-1)+m)=Dxxy(1,(p_type+1)*(q-1)+m)+w(k)*N(m)*jacob*dxxy(1,(p_type+1)*(q-1)+m);
       Dyxy(1,(p_type+1)*(q-1)+m)=Dyxy(1,(p_type+1)*(q-1)+m)+w(k)*N(m)*jacob*dyxy(1,(p_type+1)*(q-1)+m);
       Dxyy(1,(p_type+1)*(q-1)+m)=Dxyy(1,(p_type+1)*(q-1)+m)+w(k)*N(m)*jacob*dxyy(1,(p_type+1)*(q-1)+m);
       Dyyy(1,(p_type+1)*(q-1)+m)=Dyyy(1,(p_type+1)*(q-1)+m)+w(k)*N(m)*jacob*dyyy(1,(p_type+1)*(q-1)+m);

       %For displacements the coefficients of node for sharing elements are
       %added. The displacement on node does not change due to element
       Sxxx(1,sum_pos)=Sxxx(1,sum_pos)+w(k)*N(m)*jacob*sxxx(1,sum_pos);
       Syxx(1,sum_pos)=Syxx(1,sum_pos)+w(k)*N(m)*jacob*syxx(1,sum_pos);
       Sxxy(1,sum_pos)=Sxxy(1,sum_pos)+w(k)*N(m)*jacob*sxxy(1,sum_pos);
       Syxy(1,sum_pos)=Syxy(1,sum_pos)+w(k)*N(m)*jacob*syxy(1,sum_pos);
       Sxyy(1,sum_pos)=Sxyy(1,sum_pos)+w(k)*N(m)*jacob*sxyy(1,sum_pos);
       Syyy(1,sum_pos)=Syyy(1,sum_pos)+w(k)*N(m)*jacob*syyy(1,sum_pos);
       end
end
    end

sigmaxx(i)=Dxxx*tx_ele+Dyxx*ty_ele-Sxxx*u_node-Syxx*v_node;
sigmaxy(i)=Dxxy*tx_ele+Dyxy*ty_ele-Sxxy*u_node-Syxy*v_node;
sigmayy(i)=Dxyy*tx_ele+Dyyy*ty_ele-Sxyy*u_node-Syyy*v_node;
   end
    
stressp1=0.5.*(sigmaxx+sigmayy)+sqrt((0.5.*(sigmaxx-sigmayy)).^2+sigmaxy.^2);
stressp2=0.5.*(sigmaxx+sigmayy)-sqrt((0.5.*(sigmaxx-sigmayy)).^2+sigmaxy.^2);
Z=zeros(pts,1);
for i=1:pts
    
        Z(i)=max(abs(stressp1(i)),abs(stressp2(i)));
    end

[xi, yi] = meshgrid(...
    linspace(min(xcoord),max(xcoord)),...
    linspace(min(ycoord),max(ycoord)));
zi = griddata(xcoord,ycoord,Z, xi,yi);
Z_x=flip(zi,1);
Z_y=flip(zi,2);
Z_xy=flip(Z_x,2);
f3=figure(3);
hold on;
title('Stress Plot Title');
contourf(xi,yi,zi,[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24],'Linecolor','none');
c=colorbar;
c.Ticks=[0 4 8 12 16 20 24];
c.Label.String = 'Stress';
c.Label.FontSize=16;

hold off;
end;