%Enter the Elastic constants
nu=0;  Elast=1000; GS=(Elast)/(2*(1+nu));
kay=3-4*nu;
%Constants for calculations
C1=-1/(8*pi*(1-nu)*GS);
C2=-1/(4*pi*(1-nu));
%Enter the number of elements and their coordinates
X_bound= [9.5,0;
        10.5,0;
        10.5*cos(pi/64),10.5*sin(pi/64);
        10.5*cos(2*pi/64),10.5*sin(2*pi/64);
        10.5*cos(3*pi/64),10.5*sin(3*pi/64);
        10.5*cos(4*pi/64),10.5*sin(4*pi/64);
        10.5*cos(5*pi/64),10.5*sin(5*pi/64);
        10.5*cos(6*pi/64),10.5*sin(6*pi/64);
        10.5*cos(7*pi/64),10.5*sin(7*pi/64);
        10.5*cos(8*pi/64),10.5*sin(8*pi/64);
        10.5*cos(9*pi/64),10.5*sin(9*pi/64);
        10.5*cos(10*pi/64),10.5*sin(10*pi/64);
        10.5*cos(11*pi/64),10.5*sin(11*pi/64);
        10.5*cos(12*pi/64),10.5*sin(12*pi/64);
        10.5*cos(13*pi/64),10.5*sin(13*pi/64);
        10.5*cos(14*pi/64),10.5*sin(14*pi/64);
        10.5*cos(15*pi/64),10.5*sin(15*pi/64);
        10.5*cos(16*pi/64),10.5*sin(16*pi/64);
        10.5*cos(17*pi/64),10.5*sin(17*pi/64);
        10.5*cos(18*pi/64),10.5*sin(18*pi/64);
        10.5*cos(19*pi/64),10.5*sin(19*pi/64);
        10.5*cos(20*pi/64),10.5*sin(20*pi/64);
        10.5*cos(21*pi/64),10.5*sin(21*pi/64);
        10.5*cos(22*pi/64),10.5*sin(22*pi/64);
        10.5*cos(23*pi/64),10.5*sin(23*pi/64);
        10.5*cos(24*pi/64),10.5*sin(24*pi/64);
        10.5*cos(25*pi/64),10.5*sin(25*pi/64);
        10.5*cos(26*pi/64),10.5*sin(26*pi/64);
        10.5*cos(27*pi/64),10.5*sin(27*pi/64);
        10.5*cos(28*pi/64),10.5*sin(28*pi/64);
        10.5*cos(29*pi/64),10.5*sin(29*pi/64);
        10.5*cos(30*pi/64),10.5*sin(30*pi/64);
        10.5*cos(31*pi/64),10.5*sin(31*pi/64);
       0,10.5;
       0,9.5;
       9.5*cos(31*pi/64),9.5*sin(31*pi/64);
       9.5*cos(30*pi/64),9.5*sin(30*pi/64);
        9.5*cos(29*pi/64),9.5*sin(29*pi/64);
        9.5*cos(28*pi/64),9.5*sin(28*pi/64);
       9.5*cos(27*pi/64),9.5*sin(27*pi/64);
        9.5*cos(26*pi/64),9.5*sin(26*pi/64);
        9.5*cos(25*pi/64),9.5*sin(25*pi/64);
       9.5*cos(24*pi/64),9.5*sin(24*pi/64);
       9.5*cos(23*pi/64),9.5*sin(23*pi/64);
        9.5*cos(22*pi/64),9.5*sin(22*pi/64);
        9.5*cos(21*pi/64),9.5*sin(21*pi/64);
       9.5*cos(20*pi/64),9.5*sin(20*pi/64);
        9.5*cos(19*pi/64),9.5*sin(19*pi/64);
        9.5*cos(18*pi/64),9.5*sin(18*pi/64);
       9.5*cos(17*pi/64),9.5*sin(17*pi/64);
       9.5*cos(16*pi/64),9.5*sin(16*pi/64);
       9.5*cos(15*pi/64),9.5*sin(15*pi/64);
        9.5*cos(14*pi/64),9.5*sin(14*pi/64);
        9.5*cos(13*pi/64),9.5*sin(13*pi/64);
       9.5*cos(12*pi/64),9.5*sin(12*pi/64);
        9.5*cos(11*pi/64),9.5*sin(11*pi/64);
        9.5*cos(10*pi/64),9.5*sin(10*pi/64);
       9.5*cos(9*pi/64),9.5*sin(9*pi/64);
       9.5*cos(8*pi/64),9.5*sin(8*pi/64);
        9.5*cos(7*pi/64),9.5*sin(7*pi/64);
        9.5*cos(6*pi/64),9.5*sin(6*pi/64);
       9.5*cos(5*pi/64),9.5*sin(5*pi/64);
        9.5*cos(4*pi/64),9.5*sin(4*pi/64);
        9.5*cos(3*pi/64),9.5*sin(3*pi/64);
       9.5*cos(2*pi/64),9.5*sin(2*pi/64);
       9.5*cos(1*pi/64),9.5*sin(1*pi/64);];
%Coordinates of nodes
for i=1:17
    X_ele(i,1)=9.5+(i-1)*1/16;
    X_ele(i,2)=0;
end
for i=1:127
    X_ele(i+17,1)=10.5*cos(i*pi/256);
    X_ele(i+17,2)=10.5*sin(i*pi/256);
end
for i=1:17
    X_ele(i+17+127,2)=10.5-(i-1)*1/16;
    X_ele(i+17+127,1)=0;
end
for i=1:127
    X_ele(i+17+127+17,1)=9.5*cos((128-i)*pi/256);
    X_ele(i+17+127+17,2)=9.5*sin((128-i)*pi/256);
end

ele=size(X_ele,1);
X_node=X_ele;
phi_ele=zeros(ele,1);
for i=1:ele
    phi_ele(i)=atan2(X_node(i,2),X_node(i,1));
end
u_ana=zeros(ele,1);
v_ana=zeros(ele,1);
for i=1:ele
    r=sqrt(X_ele(i,1)^2+X_ele(i,2)^2);
   u_ana(i,1)=1/(8*GS)*((r*(kay+1)*cos(phi_ele(i)))+2*((1+kay)*cos(phi_ele(i))+cos(3*phi_ele(i)))/r-2*cos(3*phi_ele(i))/r^3);
   v_ana(i,1)=1/(8*GS)*((r*(kay-1)*sin(phi_ele(i)))+2*((1-kay)*sin(phi_ele(i))+sin(3*phi_ele(i)))/r-2*sin(3*phi_ele(i))/r^3);
end

u_known=zeros(17,2);
u_known(:,1)=[145:161]';
      
v_known=u_known;
tx_known=zeros(540,2);
ty_known=tx_known;
j=1;
for i=1:576
    if ((i<289)||(i>320))
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
tx_known(1:32,3)=-0.1*ones(32,1);
      
%     for i=3:4
%     r=sqrt(X_ele(i,1)^2+X_ele(i,2)^2);
%     r1=sqrt(X_ele(i+1,1)^2+X_ele(i+1,2)^2);
%     tx_known(2*i-1-4,3)=-1*(1-(1/r^2)*(1.5*cos(2*phi_ele(i))+cos(4*phi_ele(i)))+(1.5*(cos(4*phi_ele(i)))/(r^4)));
%     tx_known(2*i-4,3)=-1*(1-(1/r1^2)*(1.5*cos(2*phi_ele(i+1))+cos(4*phi_ele(i+1)))+1.5*(cos(4*phi_ele(i+1)))/(r1^4));
%     ty_known(2*i-1,3)=-1*(-(1/r^2)*(0.5*sin(2*phi_ele(i))+sin(4*phi_ele(i)))+1.5*(cos(4*phi_ele(i)))/(r^4));
%     ty_known(2*i,3)=-1*(-(1/r1^2)*(0.5*sin(2*phi_ele(i+1))+sin(4*phi_ele(i+1)))+1.5*(cos(4*phi_ele(i+1)))/(r1^4));
%     end
% 
% %  ty_known(16:21,3)=ty_known(16:21,3)*1;
% %  tx_known(10:15,3)=tx_known(10:15,3)*1;
% for i=3:4
%        r=sqrt(X_ele(i,1)^2+X_ele(i,2)^2);
%     r1=sqrt(X_ele(i+1,1)^2+X_ele(i+1,2)^2);
%     ty_known(2*i-1,3)=1*(-(1/r^2)*(0.5*cos(2*phi_ele(i))-cos(4*phi_ele(i)))-1.5*(cos(4*phi_ele(i)))/(r^4));
%     ty_known(2*i,3)=1*(-(1/r1^2)*(0.5*cos(2*phi_ele(i+1))-cos(4*phi_ele(i+1)))-1.5*(cos(4*phi_ele(i+1)))/(r1^4));
%     tx_known(2*i-1-4,3)=1*(-(1/r^2)*(0.5*sin(2*phi_ele(i))+sin(4*phi_ele(i)))+1.5*(cos(4*phi_ele(i)))/(r^4));
%     tx_known(2*i-4,3)=1*(-(1/r1^2)*(0.5*sin(2*phi_ele(i+1))+sin(4*phi_ele(i+1)))+1.5*(cos(4*phi_ele(i+1)))/(r1^4));
%  
% end

% ty_known(8:13,3)=ty_known(8:13,3)*1;
% tx_known(2:7,3)=tx_known(2:7,3)*1;
%Unknown parameters (Only nodes)       
u_unknown1=[1:144,162:288];
             
v_unknown1=[1:144,162:288];
      u_unknown=u_unknown1';
      v_unknown=v_unknown1';
      tx_unknown=zeros(32,2);
      for i=1:32
          tx_unknown(i,1)=round((144*2+i)/2);
          tx_unknown(i,2)=2-rem(i,2);
      end;
%       tx_unknown=[37,1;
%           37,2;
%           38,1;
%           38,2;
%           39,1;
%           39,2;
%           40,1;
%           40,2;
%           ];
       ty_unknown=tx_unknown;
      tx_eq_no=[37,2,38,1;
          38,2,39,1;
          39,2,40,1;
          ];
      counter=1;
      smco=1;
      for i=1:30
          tx_eq_no(counter,smco)=tx_unknown(i+1,1);
          tx_eq_no(counter,smco+1)=tx_unknown(i+1,2);
          smco=smco+2;
          if smco==4
              smco=1;
              counter=counter+1;
          end
      end
      ty_eq_no=tx_eq_no;
%      %corner_node=[1,0.301044457;
%                   4,0.25;
%                   8,0.25;
%                   12,0.25;
%                   15,0.301044457;]
node=ele;
corner_node=zeros(node,2);
for i=1:ele
    corner_node(i,:)=[i,0.5];
end;
%               for i=1:ele
%           j=i-1;
%           k=i+1;
%           if i==1
%               j=ele;
%           elseif i==ele
%                   k=1;
%           end
%           P2=X_node(j,:);
%           P1=X_node(k,:);
%           P0=X_node(i,:);
%           ang = atan2((det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
%           angr=ang;
%           if angr<0
%                   angr=-1*angr;
%           end
%           corner_node(i,1)=i;
%           corner_node(i,2)=angr/(2*pi);
%               end
%           corner_node([16;17;18],2)=1-corner_node([16;17;18],2);
    p_type=1;
  %%Internal nodes for post processing
  
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
        
            postp=0;