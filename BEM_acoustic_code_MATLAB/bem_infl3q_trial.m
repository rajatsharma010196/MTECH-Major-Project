function [DISLK,DISMK,DISMKT,DISNK]=bem_infl3q_trial(coord,ex,ey,ez,ep,n,edof,coordno,nodecoord,~,gloop,gauss,weight,gloop2,gauss2,weight2)
rev=n;
k=ep(1)/ep(2);
%Node point
x=coord(1);
y=coord(2);
z=coord(3);
%Element centroid
qx=(ex(1)+ex(2)+ex(3))/3;
qy=(ey(1)+ey(2)+ey(3))/3;
qz=(ez(1)+ez(2)+ez(3))/3;

%Calculate Area
axy=[1 0 -1;0 1 -1]*[ex;ey]';
ayz=[1 0 -1;0 1 -1]*[ey;ez]';
azx=[1 0 -1;0 1 -1]*[ez;ex]';
Area=0.5*sqrt(det(axy).^2+det(ayz).^2+det(azx).^2);
%Initialize Normals
index=edof(coordno,2:4);
eix=nodecoord(index,1);
eiy=nodecoord(index,2);
eiz=nodecoord(index,3);
if coordno>size(edof,1)
    eix=-1*eix;
end
a=[ex(2)-ex(1) ey(2)-ey(1) ez(2)-ez(1)];
b=[ex(3)-ex(1) ey(3)-ey(1) ez(3)-ez(1)];
n=[a(2)*b(3)-a(3)*b(2); a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)]; 
NORMQ=rev*n/sqrt(n'*n);
c=[eix(2)-eix(1) eiy(2)-eiy(1) eiz(2)-eiz(1)];
d=[eix(3)-eix(1) eiy(3)-eiy(1) eiz(3)-eiz(1)];
ni=[c(2)*d(3)-c(3)*d(2) c(3)*d(1)-c(1)*d(3) c(1)*d(2)-c(2)*d(1)]; 

SK=k*k;
SKO2=SK/2;
%Initialize csuml to csumn
CSUML=0;
CSUMM=0;
CSUMMT=0;
CSUMN=0;
%Check if point lies on the element 
lpon=0;
if (x==qx)&&(y==qy)&&(z==qz)
    lpon=1;
end
%Initialise VECP
%if lpon==0
    VECP=(ni')./sqrt(ni*ni');
    %VECP=coord'./norm(coord);
% else
%     VECP=NORMQ;
% end
DNPNQ=NORMQ'*VECP;

QA=[ex(1) ey(1) ez(1)];
QB=[ex(2) ey(2) ez(2)];
QC=[ex(3) ey(3) ez(3)];
QBMQC=QB-QC;
QCMQA=QC-QA;
QAMQB=QA-QB;

%The analytical integration of singular components of DISLK and DISNK
if lpon==0
    CSUML=0;
    CSUMN=0;
else
    RQBQC=norm(QBMQC);
    RQCQA=norm(QCMQA);
    RQAQB=norm(QAMQB);
    RQAP=norm(coord-QA);
    RQBP=norm(coord-QB);
    RQCP=norm(coord-QC);
    ARO=[RQAP,RQBP,RQCP];
    ARA=[RQBP,RQCP,RQAP];
    AOPP=[RQAQB,RQBQC,RQCQA];
    for j=1:3
       RO=ARO(1,j);
       RA=ARA(1,j);
       OPP=AOPP(1,j);
       if(RO<RA)
           temp=RA;
           RA=RO;
           RO=temp;
       end
       SRO=RO^2;
       SRA=RA^2;
       SOPP=OPP^2;
       A=(acos((SRA+SRO-SOPP)/(2*RA*RO)));
       B=(atan2(RA*sin(A),(RO-RA*cos(A))));
       CSUML=CSUML+(RO*sin(B)*(log(tan((B+A)/2))-log(tan(B/2))))/Area-((0i*k))*(RO*sin(B))^2*(-cot(A+B)+cot(B))/Area;
       CSUMN=CSUMN+(cos(B+A)-cos(B))/(RO*sin(B)*Area)-((0i*k^3)/3)*(RO*sin(B))^2*(-cot(A+B)+cot(B))/Area;
    end
     CSUMN=(CSUMN+SKO2*CSUML);
end

%Set quadratures
for j=1:gloop
xi=gauss(j,1);
eta=gauss(j,2);
wa=2*weight(j,1);
xg=QA(1)-xi*QAMQB(1)+eta*QCMQA(1);
yg=QA(2)-xi*QAMQB(2)+eta*QCMQA(2);
zg=QA(3)-xi*QAMQB(3)+eta*QCMQA(3);
xdis=x-xg;
ydis=y-yg;
zdis=z-zg;
R=sqrt(xdis.^2+ydis.^2+zdis.^2);
SR=R*R;
RNQ=-[xdis ydis zdis]*NORMQ/R;
RNP=[xdis ydis zdis]*VECP/R;
FPG0=1/R;
FPG0R=-1/SR;
KR=k*R;
IKR=1i*k*R;
E=exp(-IKR);
FPG=E/R;
if lpon==0
    CSUML=CSUML+wa*FPG;
 else 
     CSUML=CSUML+wa*(FPG-FPG0);
end
EOSR=E/SR;
FPGR=-(1i*(EOSR)*KR+EOSR);
WFPGR=wa*FPGR;
if lpon==0
CSUMM=CSUMM+WFPGR*RNQ;
CSUMMT=CSUMMT+WFPGR*RNP;
 else
CSUMM=CSUMM+(WFPGR-wa*FPG0R)*RNQ;
CSUMMT=CSUMM+(WFPGR-wa*FPG0R)*RNP;
 end
end
%For N
for j=1:gloop2
xi=gauss2(j,1);
eta=gauss2(j,2);
wa=2*weight2(j,1);
xg=QA(1)-xi*QAMQB(1)+eta*QCMQA(1);
yg=QA(2)-xi*QAMQB(2)+eta*QCMQA(2);
zg=QA(3)-xi*QAMQB(3)+eta*QCMQA(3);
xdis=x-xg;
ydis=y-yg;
zdis=z-zg;
R=sqrt(xdis.^2+ydis.^2+zdis.^2);
SR=R*R;
CR=R*SR;
RNQ=-[xdis ydis zdis]*NORMQ/R;
RNP=[xdis ydis zdis]*VECP/R;
RNPRNQ=RNP*RNQ;
KR=k*R;
IKR=1i*k*R;
SKR=KR^2;
E=exp(-IKR);
ECR=E/CR;
if lpon==0
    CSUMN=CSUMN+wa*ECR*(RNPRNQ*(3+3*IKR-SKR)+DNPNQ*(1+IKR)); 
    %CSUMN=CSUMN+wa*(FPGR*RNPNQ+FPGRR*RNPRNQ);   
else
     %CSUMN=CSUMN+wa*((FPGR-FPG0R)*RNPNQ+(FPGRR-FPG0RR+SKO2*FPG0)*RNPRNQ);
     %CSUMN=CSUMN+wa*((FPGR-FPG0R)*(DNPNQ+RNPRNQ)/R+(FPGRR-FPG0RR)*RNPRNQ+SKO2*FPG0);
     CSUMN=CSUMN+wa*(ECR*(RNPRNQ*(3+3*IKR-SKR)+DNPNQ*(1+IKR))-1/R^3-SKO2/R+1i*k^3*(R^2)/6);
end
end
DISLK=CSUML*Area/(4*pi);
DISMK=CSUMM*Area/(4*pi);
DISMKT=CSUMMT*Area/(4*pi);
DISNK=CSUMN*Area/(4*pi);
%-----------------------------------end-----------------------