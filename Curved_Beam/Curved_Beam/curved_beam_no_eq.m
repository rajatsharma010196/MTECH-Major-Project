%========Material constants===============%
nu=0;  Elast=1000; GS=(Elast)/(2*(1+nu));
kay=3-4*nu;
%========Constants for calculations
C1=-1/(8*pi*(1-nu)*GS);
C2=-1/(4*pi*(1-nu));
%Enter the number of elements and their coordinates
X_bound= [9.95,0;
        10.05,0;
        10.05*cos(pi/64),10.05*sin(pi/64);
        10.05*cos(2*pi/64),10.05*sin(2*pi/64);
        10.05*cos(3*pi/64),10.05*sin(3*pi/64);
        10.05*cos(4*pi/64),10.05*sin(4*pi/64);
        10.05*cos(5*pi/64),10.05*sin(5*pi/64);
        10.05*cos(6*pi/64),10.05*sin(6*pi/64);
        10.05*cos(7*pi/64),10.05*sin(7*pi/64);
        10.05*cos(8*pi/64),10.05*sin(8*pi/64);
        10.05*cos(9*pi/64),10.05*sin(9*pi/64);
        10.05*cos(10*pi/64),10.05*sin(10*pi/64);
        10.05*cos(11*pi/64),10.05*sin(11*pi/64);
        10.05*cos(12*pi/64),10.05*sin(12*pi/64);
        10.05*cos(13*pi/64),10.05*sin(13*pi/64);
        10.05*cos(14*pi/64),10.05*sin(14*pi/64);
        10.05*cos(15*pi/64),10.05*sin(15*pi/64);
        10.05*cos(16*pi/64),10.05*sin(16*pi/64);
        10.05*cos(17*pi/64),10.05*sin(17*pi/64);
        10.05*cos(18*pi/64),10.05*sin(18*pi/64);
        10.05*cos(19*pi/64),10.05*sin(19*pi/64);
        10.05*cos(20*pi/64),10.05*sin(20*pi/64);
        10.05*cos(21*pi/64),10.05*sin(21*pi/64);
        10.05*cos(22*pi/64),10.05*sin(22*pi/64);
        10.05*cos(23*pi/64),10.05*sin(23*pi/64);
        10.05*cos(24*pi/64),10.05*sin(24*pi/64);
        10.05*cos(25*pi/64),10.05*sin(25*pi/64);
        10.05*cos(26*pi/64),10.05*sin(26*pi/64);
        10.05*cos(27*pi/64),10.05*sin(27*pi/64);
        10.05*cos(28*pi/64),10.05*sin(28*pi/64);
        10.05*cos(29*pi/64),10.05*sin(29*pi/64);
        10.05*cos(30*pi/64),10.05*sin(30*pi/64);
        10.05*cos(31*pi/64),10.05*sin(31*pi/64);
       0,10.05;
       0,9.95;
       9.95*cos(31*pi/64),9.95*sin(31*pi/64);
       9.95*cos(30*pi/64),9.95*sin(30*pi/64);
        9.95*cos(29*pi/64),9.95*sin(29*pi/64);
        9.95*cos(28*pi/64),9.95*sin(28*pi/64);
       9.95*cos(27*pi/64),9.95*sin(27*pi/64);
        9.95*cos(26*pi/64),9.95*sin(26*pi/64);
        9.95*cos(25*pi/64),9.95*sin(25*pi/64);
       9.95*cos(24*pi/64),9.95*sin(24*pi/64);
       9.95*cos(23*pi/64),9.95*sin(23*pi/64);
        9.95*cos(22*pi/64),9.95*sin(22*pi/64);
        9.95*cos(21*pi/64),9.95*sin(21*pi/64);
       9.95*cos(20*pi/64),9.95*sin(20*pi/64);
        9.95*cos(19*pi/64),9.95*sin(19*pi/64);
        9.95*cos(18*pi/64),9.95*sin(18*pi/64);
       9.95*cos(17*pi/64),9.95*sin(17*pi/64);
       9.95*cos(16*pi/64),9.95*sin(16*pi/64);
       9.95*cos(15*pi/64),9.95*sin(15*pi/64);
        9.95*cos(14*pi/64),9.95*sin(14*pi/64);
        9.95*cos(13*pi/64),9.95*sin(13*pi/64);
       9.95*cos(12*pi/64),9.95*sin(12*pi/64);
        9.95*cos(11*pi/64),9.95*sin(11*pi/64);
        9.95*cos(10*pi/64),9.95*sin(10*pi/64);
       9.95*cos(9*pi/64),9.95*sin(9*pi/64);
       9.95*cos(8*pi/64),9.95*sin(8*pi/64);
        9.95*cos(7*pi/64),9.95*sin(7*pi/64);
        9.95*cos(6*pi/64),9.95*sin(6*pi/64);
       9.95*cos(5*pi/64),9.95*sin(5*pi/64);
        9.95*cos(4*pi/64),9.95*sin(4*pi/64);
        9.95*cos(3*pi/64),9.95*sin(3*pi/64);
       9.95*cos(2*pi/64),9.95*sin(2*pi/64);
       9.95*cos(1*pi/64),9.95*sin(1*pi/64);];
%Coordinates of nodes
for i=1:2
    X_ele(i,1)=9.95+(i-1)*0.1;
    X_ele(i,2)=0;
end
for i=1:2398
    X_ele(i+2,1)=10.05*cos(i*pi/4798);
    X_ele(i+2,2)=10.05*sin(i*pi/4798);
end
for i=1:2
    X_ele(i+2+2398,2)=10.05-(i-1)*0.1;
    X_ele(i+2+2398,1)=0;
end
for i=1:2398
    X_ele(i+2+2398+2,1)=9.95*cos((2399-i)*pi/4798);
    X_ele(i+2+2398+2,2)=9.95*sin((2399-i)*pi/4798);
end

ele=size(X_ele,1);
X_node=X_ele;
%%======Known values========%%
%%======displacement=node,value======%%
%%======traction=element, local node number, value=========%%
u_known=[2401,0;
    2402,0;
          ];
      
v_known=u_known;
tx_known=zeros(9598,2);
ty_known=tx_known;
j=1;
for i=1:9600
    if ((i<4801)||(i>4802))
        k=0;
        if rem(i,2)==1
            k=(i+1)/2;
        else
            k=i/2;
        end;
    tx_known(j,1)=k;
    tx_known(j,2)=2-rem(i,2);
    tx_known(j,3)=0;
    ty_known(j,1)=k;
    ty_known(j,2)=2-rem(i,2);
    ty_known(j,3)=0;
    j=j+1;
    end;
end;
tx_known(1:2,3)=-0.0001*ones(2,1);
      

%%=======Unknown Parameters=========%%
%%=======displacement=node==========%%
%%=======traction=element, local node number=======%%
u_unknown1=[1:2400,2403:4800];
             
v_unknown1=[1:2400,2403:4800];
      u_unknown=u_unknown1';
      v_unknown=v_unknown1';
      tx_unknown=[2401,1;
                  2401,2;
          ];
      ty_unknown=tx_unknown;
      tx_eq_no=[];
      ty_eq_no=tx_eq_no;
node=ele;
corner_node=zeros(node,2);
for i=1:ele
    corner_node(i,:)=[i,0.5];
end;
%%======Problem Type==========%%
p_type=1;

%%======Post-Processing========%%
%%=====1-displacement, 2-stress, anything else-no postprocessing===%%
postp=0;
k=1;
X_in=zeros(2401,2);
X_in_2=X_in;
for i=1:49
    for j=1:49
        X_in(k,1)=10.48*i/50;
        X_in(k,2)=10.48*j/50;
        X_in_2(k,2)=10.48*i/50;
        X_in_2(k,1)=10.48*j/50;
        k=k+1;
    end
end
reconi=50/10.48;
reconj=50/10.48;
dats=49;
        [in,on] = inpolygon(X_in(:,1),X_in(:,2), X_bound(:,1), X_bound(:,2));
        inon = in ;                                            
        idx = find(inon(:));                                        
        xcoord = X_in(idx,1);                                           
        ycoord = X_in(idx,2);
        [in,on] = inpolygon(X_in_2(:,1),X_in_2(:,2), X_bound(:,1), X_bound(:,2));
        inon = in | on;                                            
        idx = find(inon(:));                                        
        xcoord_2 = X_in_2(idx,1);                                           
        ycoord_2 = X_in_2(idx,2);    
        
        