%Enter the Material constants
nu=0.3/1.3;%Plane stress
Elast=1*10^5;
GS=(Elast)/(2*(1+nu));
kay=3-4*nu;
%Constants for calculations
C1=-1/(8*pi*(1-nu)*GS);
C2=-1/(4*pi*(1-nu));

%======================================================================%
%Enter the element coordinates
X_ele= [0,1;
        0,1.5;
        0,2;
        0,2.5;
        0,3;
        0,3.5;
        0,4;
        -0.5,4;
        -1,4;
        -1.5,4;
        -2,4;
        -2.5,4;
        -3,4;
        -3.5,4;
       -4,4;
       -4,3.5;
       -4,3;
       -4,2.5;
       -4,2;
       -4,1.5;
       -4,1;
       -4,0.5;
       -4,0;
       -3.5,0;
       -3,0;
       -2.5,0;
       -2,0;
       -1.5,0;
       -1,0;
    -cos(pi/32),sin(pi/32);
    -cos(2*pi/32),sin(2*pi/32);
    -cos(3*pi/32),sin(3*pi/32);
    -cos(4*pi/32),sin(4*pi/32);
    -cos(5*pi/32),sin(5*pi/32);
    -cos(6*pi/32),sin(6*pi/32);
    -cos(7*pi/32),sin(7*pi/32);
    -cos(8*pi/32),sin(8*pi/32);
    -cos(9*pi/32),sin(9*pi/32);
    -cos(10*pi/32),sin(10*pi/32);
    -cos(11*pi/32),sin(11*pi/32);
    -cos(12*pi/32),sin(12*pi/32);
    -cos(13*pi/32),sin(13*pi/32);
    -cos(14*pi/32),sin(14*pi/32);
    -cos(15*pi/32),sin(15*pi/32);
    
           ];

%Coordinates of nodes
ele=size(X_ele,1);
X_node=X_ele;

%========================================================================%
%========Known Parameters============%
u_known=[1,0;%global node, value
        2,0;
        3,0;
        4,0;
        5,0;
        6,0;
        7,0;
          ];
      
v_known=[ 
    23,0;%global node, value
    24,0;
    25,0;
    26,0;
    27,0;
    28,0;
    29,0;
          ];
      
      tx_known=[ 
           7,1,0;%element, local node, value
           7,2,0;
           8,1,0;
           8,2,0;
           9,1,0;
           9,2,0;
           10,1,0;
           10,2,0;
           11,1,0;
           11,2,0;
           12,1,0;
           12,2,0;
           13,1,0;
           13,2,0;
           14,1,0;
           14,2,0;
           15,1,-20;
           15,2,-20;
           16,1,-20;
           16,2,-20;
           17,1,-20;
           17,2,-20;
           18,1,-20;
           18,2,-20;
           19,1,-20;
           19,2,-20;
           20,1,-20;
           20,2,-20;
           21,1,-20;
           21,2,-20;
           22,1,-20;
           22,2,-20;
           23,1,0;
           23,2,0;
           24,1,0;
           24,2,0;
           25,1,0;
           25,2,0;
           26,1,0;
           26,2,0;
           27,1,0;
           27,2,0;
           28,1,0;
           28,2,0;
           29,1,0;
           29,2,0;
           30,1,0;
           30,2,0;
           31,1,0;
           31,2,0;
           32,1,0;
           32,2,0;
           33,1,0;
           33,2,0;
           34,1,0;
           34,2,0;
           35,1,0;
           35,2,0;
           36,1,0;
           36,2,0;
           37,1,0;
           37,2,0;
           38,1,0;
           38,2,0;
           39,1,0;
           39,2,0;
           40,1,0;
           40,2,0;
           41,1,0;
           41,2,0;
           42,1,0;
           42,2,0;
           43,1,0;
            ];
   
      
      ty_known=[ 1,1,0;%elemnt, local node, value
            1,2,0;
            2,1,0;
            2,2,0;
            3,1,0;
            3,2,0;
          4,1,0;
            4,2,0;
          5,1,0;
          5,2,0
           6,1,0;
           6,2,0
           7,1,20;
           7,2,20;
           8,1,20;
           8,2,20;
           9,1,20;
           9,2,20
           10,1,20;
           10,2,20;
           11,1,20;
           11,2,20;
           12,1,20;
           12,2,20;
           13,1,20;
           13,2,20;
           14,1,20;
           14,2,20;
           15,1,0;
           15,2,0;
           16,1,0;
           16,2,0;
           17,1,0;
           17,2,0;
           18,1,0;
           18,2,0;
           19,1,0;
           19,2,0;
           20,1,0;
           20,2,0;
           21,1,0;
           21,2,0;
           22,1,0;
           22,2,0;
           29,1,0;
           29,2,0;
           30,1,0;
           30,2,0;
           31,1,0;
           31,2,0;
           32,1,0;
           32,2,0;
           33,1,0;
           33,2,0;
           34,1,0;
           34,2,0;
           35,1,0;
           35,2,0;
           36,1,0;
           36,2,0;
           37,1,0;
           37,2,0;
           38,1,0;
           38,2,0;
           39,1,0;
           39,2,0;
           40,1,0;
           40,2,0;
           41,1,0;
           41,2,0;
           42,1,0;
           42,2,0;
           43,1,0;
            ];
 %===============================================================%       
%Unknown parameters

u_unknown1=[8:1:44;%global nodes
          ];
      u_unknown=u_unknown1';
v_unknown1=[1:1:22,30:1:44;%global nodes
          ]; 
      v_unknown=v_unknown1';
      tx_unknown=[1,1;1,2;%element, local node number
          2,1;2,2;
          3,1;3,2;
          4,1;4,2;
          5,1;5,2;
          6,1;6,2;
                  
                  ];
      ty_unknown=[
          23,1;%element, local node number
          23,2;
          24,1;
          24,2;
          25,1;
          25,2;
          26,1;
          26,2;
          27,1;
          27,2;
          28,1;
          28,2;
          ];
      %common nodes(unknown traction) between two elements
      tx_eq_no=[1,2,2,1;
          2,2,3,1;
          3,2,4,1;
          4,2,5,1;
          5,2,6,1;];
      ty_eq_no=[23,2,24,1;
          24,2,25,1;
          25,2,26,1;
          26,2,27,1;
          27,2,28,1;
                ];
node=ele;
corner_node=zeros(node,2);
for i=1:ele
    corner_node(i,:)=[i,0.5];
end;
%======================================================================%
%Problem type
%0-constant
%1-linear
p_type=1;

%=====================================================================%
 %%Post processing
 %1-displacement
 %2-stress
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
X_in=zeros(841,2);
for i=1:49
    for j=1:49
        X_in(k,1)=-4*i/50;
        X_in(k,2)=4*j/50;
        k=k+1;
    end
end
reconi=-50/4;
reconj=50/4;
dats=49;
        [in,on] = inpolygon(X_in(:,1),X_in(:,2), X_bound(:,1), X_bound(:,2));
        inon = in;                                            
        idx = find(inon(:));                                        
        xcoord = X_in(idx,1);                                           
        ycoord = X_in(idx,2);    
    