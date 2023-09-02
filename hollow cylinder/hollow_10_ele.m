%Enter the Material constants
nu=0.25;  Elast=200000; GS=(Elast)/(2*(1+nu));
kay=3-4*nu
%Constants for calculations
C1=-1/(8*pi*(1-nu)*GS);
C2=-1/(4*pi*(1-nu));
%Enter the number of elements and their coordinates
X_ele=zeros(32,2);
for i=1:2
    X_ele(i,1)=10+15*(i-1);
    X_ele(i,2)=0;
end
for i=1:10
    X_ele(2+i,1)=25*cos((i)*pi/20);
    X_ele(2+i,2)=25*sin((i)*pi/20);
end
for i=1:1
    X_ele(12+i,1)=0;
    X_ele(12+i,2)=25-15*(i);
end
for i=1:19
    X_ele(12+1+i,1)=10*cos(pi/2-(i)*pi/40);
    X_ele(12+1+i,2)=10*sin(pi/2-(i)*pi/40);
end
%Coordinates of nodes
ele=size(X_ele,1);
X_node=X_ele;

%%=========Known Values============%%
%%=Displacement = Global Node Number, Value======%
%%=Traction = Element Number, Local Node Number, Value======% 
u_known=zeros(2,2);
for i=1:2;
    u_known(i,1)=12+i-1;
end
v_known=zeros(2,2);
for i=1:2;
    v_known(i,1)=1+i-1;
end
      
for i=1:11
    tx_known(2*i-1,1)=i;
    tx_known(2*i,1)=i;
    tx_known(2*i-1,2)=1;
    tx_known(2*i,2)=2;
    tx_known(2*i-1,3)=0;
    tx_known(2*i,3)=0;
end;
tx_known(1,3)=0;
for i=12:31
    tx_known(2*i-1,1)=i+1;
    tx_known(2*i,1)=i+1;
    tx_known(2*i-1,2)=1;
    tx_known(2*i,2)=2;
    tx_known(2*i-1,3)=100*cos(pi/2-(i-12)*pi/40);
    tx_known(2*i,3)=100*cos(pi/2-(i-11)*pi/40);
end;
for i=3:11
    ty_known(2*i-1-4,1)=i;
    ty_known(2*i-4,1)=i;
    ty_known(2*i-1-4,2)=1;
    ty_known(2*i-4,2)=2;
    ty_known(2*i-1-4,3)=0;
    ty_known(2*i-4,3)=0;
end;

for i=12:31
    ty_known(2*i-1-4,1)=i+1;
    ty_known(2*i-4,1)=i+1;
    ty_known(2*i-1-4,2)=1;
    ty_known(2*i-4,2)=2;
    ty_known(2*i-1-4,3)=100*sin(pi/2-(i-12)*pi/40);
    ty_known(2*i-4,3)=100*sin(pi/2-(i-11)*pi/40);
end;
%%=========UnKnown Values============%%
%%=Displacement = Global Node Number======%
%%=Traction = Element Number, Local Node Number======% 
u_unknown=[1:11,14:32
          ]';      
v_unknown=[3:32
          ]'; 
for i=1:1
    tx_unknown(2*i-1,1)=11+i;
    tx_unknown(2*i,1)=11+i;
    tx_unknown(2*i-1,2)=1;
    tx_unknown(2*i,2)=2; 
end
    
for i=1:1
    ty_unknown(2*i-1,1)=i;
    ty_unknown(2*i,1)=i;
    ty_unknown(2*i-1,2)=1;
    ty_unknown(2*i,2)=2; 
end
      tx_eq_no=[];
      ty_eq_no=[];
node=ele;
corner_node=zeros(node,2);
for i=1:ele
    corner_node(i,:)=[i,0.5];
end;

%%====Problem Type=1(Linear),0(Constant)=====%
p_type=1;

%%==============Post Processing=1(displacement),2(stress)=============%%
postp=0;
X_bound=X_ele;
k=1;
X_in=zeros(1521,2);
X_in_2=X_in;
for i=1:39
    for j=1:39
        X_in(k,1)=25*i/40;
        X_in(k,2)=25*j/40;
        X_in_2(k,2)=25*i/200;
        X_in_2(k,1)=25*j/200;
        k=k+1;
    end
end
reconi=40/25;
reconj=40/25;
dats=39;

        [in,on] = inpolygon(X_in(:,1),X_in(:,2), X_bound(:,1), X_bound(:,2));%Checks whether in bounds or not
        inon = in;                                            
        idx = find(inon(:));                                        
        xcoord = X_in(idx,1);                                           
        ycoord = X_in(idx,2);
        
            