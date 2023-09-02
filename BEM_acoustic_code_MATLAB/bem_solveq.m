function [pr,nv]=bem_solveq(G,H,bcpr,bcnv,bcim,bcpi,bcvi,mu)
% [pr,nv]=bem_solveq(P,V,bcpr,bcnv,bcim)

n_size=size(G);
fpdof=[1:n_size]'; fvdof=[1:n_size]'; fidof=[1:n_size]';
[rowp]=size(bcpr,1);
[rowv]=size(bcnv,1);
[rowi]=size(bcim,1);
pr=zeros(size(fpdof)); nv=zeros(size(fvdof));
%======Segregating knowns and unknowns===============% 
%If pressure BC exists
if rowp~=0
ppdof=bcpr(:,1);
prp=bcpr(:,2);
fpdof(ppdof)=[];
%If velocity and pressure BC exist
if rowv~=0
pvdof=bcnv(:,1);
nvp=bcnv(:,2);
fvdof(pvdof)=[];
%If impedance BC also exists
if rowi~=0
pidof=bcim(:,1);
imp=bcim(:,2);
HG=G;
HG(:,pvdof)=0;
for s=1:rowi
HG(:,pidof(s))=HG(:,pidof(s))-H(:,pidof(s)).*imp(s);
end
HH=H;
HH(:,fvdof)=0;
x=(HG-HH)\(H(:,ppdof)*prp-G(:,pvdof)*nvp);
nv(pvdof)=nvp;
nv(pidof)=x(pidof);
nv(fvdof)=x(fvdof);
pr(ppdof)=prp;
pr(fpdof)=x(fpdof);
pr(pidof)=nv(pidof).*imp;
else
HG=G;
HG(:,pvdof)=0;
HH=H;
HH(:,ppdof)=0;
x=(HG-HH)\(H(:,ppdof)*prp-G(:,pvdof)*nvp);
pr(ppdof)=prp;
pr(fpdof)=x(fpdof);
nv(pvdof)=nvp;
nv(fvdof)=x(fvdof);
end
elseif rowi~=0
pidof=bcim(:,1);
imp=bcim(:,2);
HG=G;
for s=1:rowi
HG(:,pidof(s))=HG(:,pidof(s))-H(:,pidof(s)).*imp(s);
end
x=HG\(H(:,ppdof)*prp);
nv=x;
pr(ppdof)=prp;
pr(pidof)=nv(pidof).*imp;
else
x=G\H*prp;
pr=prp;
nv=x;
end
else
if rowv~=0    
pvdof=bcnv(:,1);
nvp=bcnv(:,2);
fvdof(pvdof)=[];
if rowi~=0
pidof=bcim(:,1);
imp=bcim(:,2);
HH=H;
for s=1:rowi
HH(:,pidof(s))=HH(:,pidof(s))-G(:,pidof(s))./imp(s);
end
x=HH\(G(:,pvdof)*nvp);
pr=x;
nv(pvdof)=nvp;
nv(pidof)=pr(pidof)./imp;
else
x=H\(G*(nvp));
nv=nvp;
pr=x;
end
else
x=H\(bcpi(:,2)-mu*bcvi(:,2));
nv=0*bcvi(:,2);
pr=x;
end
end
%-----------------------------------end--------------------------------