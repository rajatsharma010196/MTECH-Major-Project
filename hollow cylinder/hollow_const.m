%Enter the Elastic constants
nu=0.25;  Elast=200000; GS=(Elast)/(2*(1+nu));
kay=3-4*nu;
%Constants for calculations
C1=-1/(8*pi*(1-nu)*GS);
C2=-1/(4*pi*(1-nu));
%%====Element Coordinates===========%%
X_ele=zeros(160,2);
for i=1:41
    X_ele(i,1)=10+15*(i-1)/40;
    X_ele(i,2)=0;
end
for i=1:40
    X_ele(41+i,1)=25*cos((i)*pi/80);
    X_ele(41+i,2)=25*sin((i)*pi/80);
end
for i=1:40
    X_ele(81+i,1)=0;
    X_ele(81+i,2)=25-15*(i)/40;
end
for i=1:39
    X_ele(121+i,1)=10*cos(pi/2-(i)*pi/80);
    X_ele(121+i,2)=10*sin(pi/2-(i)*pi/80);
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

%%==========Known Parameters=============%%
u_known=zeros(40,2);
for i=1:40;
    u_known(i,1)=81+i-1;
end
v_known=zeros(40,2);
for i=1:40;
    v_known(i,1)=1+i-1;
end
      
for i=1:80
    tx_known(i,1)=i;
    tx_known(i,2)=1;
    tx_known(i,3)=0;
end;
tx_known(1,3)=0
for i=81:120
    tx_known(i,1)=i+40;
    tx_known(i,2)=1;
    tx_known(i,3)=100*cos(pi/2-(2*i-161)*pi/160);
end;
for i=41:120
    ty_known(i-40,1)=i;
    ty_known(i-40,2)=1;
    ty_known(i-40,3)=0;
   end;

for i=121:160
    ty_known(i-40,1)=i;
    ty_known(i-40,2)=1;
    ty_known(i-40,3)=100*sin(pi/2-(2*i-241)*pi/160);
  
end;

%%============Unknown Parameters==============================%%
u_unknown=[1:80,121:160
          ]';      
v_unknown=[41:160
          ]'; 
for i=1:40
    tx_unknown(i,1)=80+i;
    tx_unknown(i,2)=1;
    
end
    
for i=1:40
    ty_unknown(i,1)=i;
    ty_unknown(i,2)=1;
     
end
      tx_eq_no=[];
      ty_eq_no=[];

node=ele;
corner_node=zeros(node,2);
for i=1:ele
    corner_node(i,:)=[i,0.5];
end;
%%=======Problem type===============%%
     p_type=0;%0 for constant, 1 for linear
     
%%========Post Processing===========%%
  X_bound= [10,0;
        17.5,0;
        24,0;
        24*cos(pi/64),24*sin(pi/64);
        24*cos(2*pi/64),24*sin(2*pi/64);
        24*cos(3*pi/64),24*sin(3*pi/64);
        24*cos(4*pi/64),24*sin(4*pi/64);
        24*cos(5*pi/64),24*sin(5*pi/64);
        24*cos(6*pi/64),24*sin(6*pi/64);
        24*cos(7*pi/64),24*sin(7*pi/64);
        24*cos(8*pi/64),24*sin(8*pi/64);
        24*cos(9*pi/64),24*sin(9*pi/64);
        24*cos(10*pi/64),24*sin(10*pi/64);
        24*cos(11*pi/64),24*sin(11*pi/64);
        24*cos(12*pi/64),24*sin(12*pi/64);
        24*cos(13*pi/64),24*sin(13*pi/64);
        24*cos(14*pi/64),24*sin(14*pi/64);
        24*cos(15*pi/64),24*sin(15*pi/64);
        24*cos(16*pi/64),24*sin(16*pi/64);
        24*cos(17*pi/64),24*sin(17*pi/64);
        24*cos(18*pi/64),24*sin(18*pi/64);
        24*cos(19*pi/64),24*sin(19*pi/64);
        24*cos(20*pi/64),24*sin(20*pi/64);
        24*cos(21*pi/64),24*sin(21*pi/64);
        24*cos(22*pi/64),24*sin(22*pi/64);
        24*cos(23*pi/64),24*sin(23*pi/64);
        24*cos(24*pi/64),24*sin(24*pi/64);
        24*cos(25*pi/64),24*sin(25*pi/64);
        24*cos(26*pi/64),24*sin(26*pi/64);
        24*cos(27*pi/64),24*sin(27*pi/64);
        24*cos(28*pi/64),24*sin(28*pi/64);
        24*cos(29*pi/64),24*sin(29*pi/64);
        24*cos(30*pi/64),24*sin(30*pi/64);
        24*cos(31*pi/64),24*sin(31*pi/64);
       0,24;
       0,17.5;
       0,11;
       11*cos(31*pi/64),11*sin(31*pi/64);
       11*cos(30*pi/64),11*sin(30*pi/64);
        11*cos(29*pi/64),11*sin(29*pi/64);
        11*cos(28*pi/64),11*sin(28*pi/64);
       11*cos(27*pi/64),11*sin(27*pi/64);
        11*cos(26*pi/64),11*sin(26*pi/64);
        11*cos(25*pi/64),11*sin(25*pi/64);
       11*cos(24*pi/64),11*sin(24*pi/64);
       11*cos(23*pi/64),11*sin(23*pi/64);
        11*cos(22*pi/64),11*sin(22*pi/64);
        11*cos(21*pi/64),11*sin(21*pi/64);
       11*cos(20*pi/64),11*sin(20*pi/64);
        11*cos(19*pi/64),11*sin(19*pi/64);
        11*cos(18*pi/64),11*sin(18*pi/64);
       11*cos(17*pi/64),11*sin(17*pi/64);
       11*cos(16*pi/64),11*sin(16*pi/64);
       11*cos(15*pi/64),11*sin(15*pi/64);
        11*cos(14*pi/64),11*sin(14*pi/64);
        11*cos(13*pi/64),11*sin(13*pi/64);
       11*cos(12*pi/64),11*sin(12*pi/64);
        11*cos(11*pi/64),11*sin(11*pi/64);
        11*cos(10*pi/64),11*sin(10*pi/64);
       11*cos(9*pi/64),11*sin(9*pi/64);
       11*cos(8*pi/64),11*sin(8*pi/64);
        11*cos(7*pi/64),11*sin(7*pi/64);
        11*cos(6*pi/64),11*sin(6*pi/64);
       11*cos(5*pi/64),11*sin(5*pi/64);
        11*cos(4*pi/64),11*sin(4*pi/64);
        11*cos(3*pi/64),11*sin(3*pi/64);
       11*cos(2*pi/64),11*sin(2*pi/64);
       11*cos(1*pi/64),11*sin(1*pi/64);];
k=1;
X_in=zeros(361,2);
X_in_2=X_in;
for i=1:19
    for j=1:19
        X_in(k,1)=25*i/20;
        X_in(k,2)=25*j/20;
        X_in_2(k,2)=25*i/70;
        X_in_2(k,1)=25*j/70;
        k=k+1;
    end
end
reconi=20/25;
reconj=20/25;
dats=19;
        [in,on] = inpolygon(X_in(:,1),X_in(:,2), X_bound(:,1), X_bound(:,2));
        inon = in | on;                                            
        idx = find(inon(:));                                        
        xcoord = X_in(idx,1);                                           
        ycoord = X_in(idx,2);
        [in,on] = inpolygon(X_in_2(:,1),X_in_2(:,2), X_bound(:,1), X_bound(:,2));
        inon = in | on;                                            
        idx = find(inon(:));                                        
        xcoord_2 = X_in_2(idx,1);                                           
        ycoord_2 = X_in_2(idx,2);    
        
            postp=2;%1 for displacement, 2 for stress