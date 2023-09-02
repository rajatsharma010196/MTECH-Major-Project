%====================================================================%
%Enter the Elastic constants
nu=0.3/1.3;  Elast=1*10^5;
GS=(Elast)/(2*(1+nu));
kay=3-4*nu;
%Constants for calculations
C1=-1/(8*pi*(1-nu)*GS);
C2=-1/(4*pi*(1-nu));

%===================================================================%
%Enter the element coordinates
X_ele=zeros(60,2);
for i=1:17
    X_ele(i,1)=0;
    X_ele(i,2)=1+(i-1)*3/16;
end
for i=1:6
    X_ele(i+17,1)=-0-(i)*4/6;
    X_ele(i+17,2)=4;
end
for i=1:6
    X_ele(i+17+6,2)=4-(i)*4/6;
    X_ele(i+17+6,1)=-4;
end
for i=1:16
    X_ele(i+17+6*2,1)=-4+(i)*3/16;
    X_ele(i+17+12,2)=0;
end
for i=1:15
    X_ele(i+17+28,1)=-cos(i*pi/32);
    X_ele(i+17+28,2)=sin(i*pi/32);
end

%Coordinates of nodes
ele=size(X_ele,1);
X_node=zeros(ele,2);
for i=1:ele
    p=i+1;
    if(i==ele)
        p=1;
    end
    X_node(i,:)=0.5.*(X_ele(i,:)+X_ele(p,:));
end

%======================================================================%
%Known Parameters
%displacement=global node number, value
u_known=zeros(16,2);
for i=1:16
    u_known(i,1)=i;
    u_known(i,2)=0;
end
v_known=zeros(15,2);
for i=30:44
    v_known(i-29,1)=i;
    v_known(i-29,2)=0;
end
%traction= Element Number, Local Node Number, Value
tx_known=zeros(44,3);
ty_known=zeros(44,3);

for i=17:22
    tx_known(i-16,1)=i;
    tx_known(i-16,2)=1;
    tx_known(i-16,3)=0;
end
for i=23:28
    tx_known(i-16,1)=i;
    tx_known(i-16,2)=1;
    tx_known(i-16,3)=-20;
end
for i=29:60
    tx_known(i-16,1)=i;
    tx_known(i-16,2)=1;
    tx_known(i-16,3)=0;
end
for i=1:16
    ty_known(i,1)=i;
    ty_known(i,2)=1;
    ty_known(i,3)=0;
end

for i=17:22
    ty_known(i,1)=i;
    ty_known(i,2)=1;
    ty_known(i,3)=20;
end

for i=23:28
    ty_known(i,1)=i;
    ty_known(i,2)=1;
    ty_known(i,3)=0;
end

for i=45:60
    ty_known(i-16,1)=i;
    ty_known(i-16,2)=1;
    ty_known(i-16,3)=0;
end

%========================================================================%
%Unknown parameters 
%Displacemnt=Global Node Numbers
u_unknown1=[17:1:60];
u_unknown=u_unknown1';
v_unknown1=[1:28,45:60];
v_unknown=v_unknown1';
%Traction=Element number, Local Node Number
tx_unknown=zeros(16,2);
ty_unknown=zeros(16,2);
for i=1:16
    tx_unknown(i,1)=i;
    tx_unknown(i,2)=1;
end
for i=29:44
    ty_unknown(i-28,1)=i;
    ty_unknown(i-28,2)=1;
end

      tx_eq_no=[];
      ty_eq_no=[
                ];
node=ele;
corner_node=zeros(node,2);
for i=1:ele
    corner_node(i,:)=[i,0.5];
end;

%========================================================================%
%Problem Type
%0-constant
%1-Linear
p_type=0;

%========================================================================%
%Post Processing
%1-displacement
%2-stress
%%Internal nodes for post processing
   postp=2;  
X_bound=zeros(100,2);
 X_bound(1:5,:)= [0.05,1.05;
           0.05,3.95;
           -3.95,3.95;
           -3.95,0.05;
           -1.05,0.05;];
  for i=6:104
      X_bound(i,1)=-1.05*cos((i-5)*pi/200);
      X_bound(i,2)=1.05*sin((i-5)*pi/200);
  end;
k=1;
X_in=zeros(2401,2);
for i=1:49
    for j=1:49
        X_in(k,1)=-4*i/40;
        X_in(k,2)=4*j/40;
        k=k+1;
    end
end
reconi=-40/4;
reconj=40/4;
dats=39;
        [in,on] = inpolygon(X_in(:,1),X_in(:,2), X_bound(:,1), X_bound(:,2));
        inon = in;                                            
        idx = find(inon(:));                                        
        xcoord = X_in(idx,1);                                           
        ycoord = X_in(idx,2);    
         