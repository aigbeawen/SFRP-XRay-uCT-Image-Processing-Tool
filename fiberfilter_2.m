function [I,E,T]=fiberfilter_2(MV,Fibers,fname)
Qtol=.975; Qvtol=.9995;  Xtol1=5.5; Xtol2=17.5; dXtol=5; Nf=5;
pxlim1=250; pxlim2=500; fname=['.\Data\',fname,'.mat'];
% adjust parameters for imajdustn and fibermetric as needed
nf=max(Fibers,[],'all'); sz=size(Fibers); 
if isfile(fname)
    load(fname,'I','E','T','jf'); kf=jf+1; ds=max(I(:)); de=max(E(:));
else
    I=zeros(sz); E=zeros(sz); T=zeros(sz); kf=1; ds=0; de=0; 
end
%
for jf=kf:nf
   ftf0= Fibers==jf;  ftf=imfill(ftf0,4,'holes');
   ctf=sum(ftf,'all'); MV(logical(ftf-ftf0))=1.; MVf=MV.*uint8(ftf);
   ftj=find(ftf); [xf,yf,zf]=ind2sub(sz,ftj) ; Xf=[xf,yf,zf];
   Rsf=regionprops3(ftf,'PrincipalAxisLength','Solidity'); Avf=Rsf.Solidity; 
   Lsf=Rsf.PrincipalAxisLength; Arf=Lsf(1)/sqrt(prod(Lsf(2:3)));
   if ctf>pxlim1
        Jfs=imsplit(MVf)  ;  fts=Jfs>0;  ftd=imfill(fts,4,"holes");
        cts=sum(fts,'all');   Ds = bwskel(ftd); Ds = bwmorph3(Ds,'remove');
        Ds = bwmorph3(Ds,'clean'); Ns=bwlabeln(Ds,26) ; ns=max(Ns,[],'all');
        if ns==1 || (ctf<pxlim2 && Arf>3.5) || cts/ctf<.15
            Ds=Dskel(ftf,Avf); Ns=bwlabeln(Ds,26); epts=bwmorph3(Ds,"endpoints"); 
        else
            [~,Ns,epts]=imstitch(ftf,Ns,[],.95,Qtol,Xtol1,Xtol2,0.);
            [~,Ns,epts]=imextend(Ns,epts,ftf,Xf,Qvtol,Qtol,Xtol1,Xtol2,0.25);
            %
            ns=max(Ns(:));freqs=histc(Ns(:),1:ns); bs=find(freqs>=Nf);
            Ns(~ismember(Ns,bs))=0; [~,ia,ic]=unique(Ns(:),'sorted');
            ns=length(ia);ls=0:ns-1;Ns=reshape(ls(ic),sz);epts=epts.*(Ns>0); 
        end
    else
        Ds=Dskel(ftf,Avf); Ns=bwlabeln(Ds,26); epts =bwmorph3(Ds,"endpoints"); 
    end
     nskel=find(Ns); ns=max(Ns(:)) ; I(nskel)=ds+Ns(nskel); ds=ds+ns; 
     [xd,yd,zd]=ind2sub(sz,nskel)  ; Xd=[xd,yd,zd];
    %
    if numel(nskel)==1
       npts=nskel; 
    else
        npts=find(epts);
        if isempty(npts)||numel(npts)<2
            Tx=round(regiontips(Xd,Ns(nskel))); 
            npts=sum((Tx-1).*sz(1).^(0:2),2)+1;
        end
    end
    [xe,ye,ze]=ind2sub(sz,npts) ; Xe=[xe,ye,ze]; 
    if ctf>pxlim1
        [fte,nte,ke]=ndxtip(ftj,Xf,Xe,dXtol); E(fte)=de+nte; de=de+ke;
    else
        Nprop=regionprops3(Ns,'PrincipalAxisLength','Solidity'); 
        Ls=Nprop.PrincipalAxisLength; Arf=Ls(1)/sqrt(prod(Ls(2:3)));
        Avf=Nprop.Solidity;
        if Arf>2.5 && Avf>.35
            [fte,nte,ke]=ndxtip(ftj,Xf,Xe,dXtol); E(fte)=de+nte; de=de+ke;
        else
            fte=ftj; E(fte)=de+1; de=de+1;
        end
    end 
    T(npts)=1;  % [xp,yp,zp]=ind2sub(sz,fte) ; Xp=[xp,yp,zp];
    save(fname,'I','E','T','jf');
    %
    % figure(1); clf; hold on;
    % plot3(xf,yf,zf,'k.'); plot3(xd,yd,zd,'r.');
    % plot3(xe,ye,ze,'bp'); plot3(xp,yp,zp,'c.');
    % view(60,30); axis equal; pause(.5);
    fprintf('Processing - %2d of %2d complete\n',jf,nf);
end
%
end
%
function [fte,nte,nk]=ndxtip(ftj,Xf,Xe,dXtol)
dXf=pdist2(Xf,Xe,'euclidean'); [dXj,jtj]=min(dXf,[],2);
itj=find(dXj<dXtol); fte=ftj(itj); 
ntj=size(Xe,1);  kfj=jtj(itj); ne=(1:ntj); ntj=ne(kfj);
[~,ja,jc]=unique(ntj,'stable'); nk=length(ja); ntk=1:nk; nte=ntk(jc);
end
%
function [Ds,Ns,epts]=imextend(Ns,epts,ftf,Xf,varargin)
[Qvtol,Qtol,Xtol1,Xtol2,extol]=varargin{:};
Ds=Ns>0; ns=max(Ns(:)); sz=size(ftf);
%
for jt=1:ns
    njt=find((Ns==jt).*epts); nkt=length(njt);
    %
    if nkt==2
        [xj,yj,zj]=ind2sub(sz,njt); Xj=[xj,yj,zj]; Xp=Xj;
        %
        iNs=(Ns~=jt).*Ns;   iNf=find(iNs>0|~ftf);
        [xk,yk,zk]=ind2sub(sz,iNf); iXk=[xk,yk,zk] ;
        %
        for jk=1:nkt
            nxj=(-1)^jk*diff(Xj); nxj=nxj/norm(nxj); dXkj=iXk-Xj(jk,:);
            nXkj=vecnorm(dXkj,2,2); nkj=dXkj./nXkj; cosQk=sum(nkj.*nxj,2);
            jtk=1; xbol=1;
            while xbol && jtk<30
                ckbol=cosQk>Qvtol^jtk;
                [dXj,kXj]=min(nXkj(ckbol));
                xbol=isempty(dXj);
                if ~xbol
                    jXk=iXk(ckbol,:); rXk=jXk(kXj,:);
                    % [Rx,Ry,Rz]=bresenham3(rXk,Xj(jk,:)); 
                    % jkt=sub2ind(sz,Rx,Ry,Rz);
                    Tx=[rXk;Xj(jk,:)]; jpts=sum((Tx-1).*sz(1).^(0:2),2)+1;
                    jkt=impatch(jpts,sz); jkt=jkt{:};
                    %
                    ktf=sqrt(sum(ftf(jkt))/numel(jkt));
                    xbol=ktf<.75;
                end
                jtk=jtk+1;
                %
            end
            if ~xbol
                dXfj=Xf-Xj(jk,:); nXfj=vecnorm(dXfj,2,2); nfj=dXfj./nXfj;
                cosQj=sum(nfj.*nxj,2); cjbol=cosQj>.98*Qvtol;
                %
                nXa=find((nXfj<dXj).*cjbol);
                if ~isempty(nXa)
                    [~,nXb]=max(nXfj(nXa));
                    if ~isempty(nXb)
                        Xp(jk,:)=Xf(nXa(nXb),:);
                    end
                end
            end
        end
        nptj=sum((Xp-1).*sz(1).^(0:2),2)+1;
        %
        for jk=1:nkt
            % [Px,Py,Pz]=bresenham3(Xj(jk,:),Xp(jk,:)); Zp=[Px',Py',Pz'];
            % fzp=sum((Zp-1).*sz.^(0:2),2)+1; 
            Tx=[Xj(jk,:);Xp(jk,:)]; jpts=sum((Tx-1).*sz(1).^(0:2),2)+1;
            fzp=impatch(jpts,sz); fzp=fzp{:}; 
            Zp=[]; [Zp(:,1),Zp(:,2),Zp(:,3)]=ind2sub(sz,fzp(:));
            %
            fds=find(Ds); [Cx,Cy,Cz]=ind2sub(sz,fds)  ; Zc=[Cx,Cy,Cz];
            dZ=pdist2(Zp,Zc,'euclidean');
            %
            [mzp,kds]=min(dZ,[],2);xds=find(mzp<3.); yds=kds(xds); 
            Ds(fzp(xds))=0; Ds(fds(yds))=0;
            Ns(fzp(xds))=0; Ns(fds(yds))=0;
            Zcd=Zp; Zcd(xds,:)=round(.5*(Zp(xds,:)+Zc(yds,:)));
            %
            %  Zcd=Zp(logical(1-any(dZ<2.,2)),:);
            Px=Zcd(:,1); Py=Zcd(:,2); Pz=Zcd(:,3);
            npj=sub2ind(sz,Px,Py,Pz); Ds(npj)=1; Ns(npj)=jt;
            epts(njt(jk))=0; Zpc=[Zp(xds,:);Zc(yds,:)];
            nptx=sum((Zpc-1).*sz(1).^(0:2),2)+1; epts(nptx)=0; 
            % ??
            kcd=knnsearch(Zcd,Xp(jk,:));
            nptk=sum((Zcd(kcd,:)-1).*sz(1).^(0:2),2)+1;  epts(nptk)=1; 
        end
    end
end
%
fs=histcounts(Ns(:),1:ns); epts(ismember(Ns,find(fs==1)))=1;npts=find(epts);
[Fx,Fy,Fz]=ind2sub(sz,npts); PX=[Fx,Fy,Fz];
pdx=pdist2(PX,PX,'euclidean'); ndx=size(pdx,1); pdx(logical(eye(ndx)))=inf;
[mpx,rdx]=min(pdx,[],2); npdx=[(1:ndx)',rdx];
bdx=logical((mpx>sqrt(3)).*(mpx<3.)); npdx=npdx(bdx,:); 
npdx=unique(sort(npdx,2),'rows')    ; kpts=npts(npdx); 
kdx=size(kpts,1) ; kdx=(numel(kpts)>2)*(kdx-1)+1;
if ~isempty(npdx)
    jpts=impatch(kpts,sz); Ds([jpts{:}])=1;  epts(kpts(:))=0;
    for jk=1:kdx
        Ns(jpts{jk})=Ns(kpts(jk,1)); % check
    end
end
[~,ia,ib]=unique(Ns(:),'sorted'); js=length(ia); 
ls=0:js-1; Ns=reshape(ls(ib),sz);
%
[Ds,Ns,epts]=imstitch(ftf,Ns,epts,1.,Qtol,Xtol1,Xtol2,extol);
% npts=endpnts(Ds); epts =zeros(sz); epts(npts)=1;
end
%
function [Ds,Ns,epts]=imstitch(ftf,Ns,epts,nfac,varargin)
[Qtol,Xtol1,Xtol2,extol]=varargin{:}; pmtj=nchoosek(1:4,2);
sz=size(Ns); Ds=Ns>0; npts=find(epts);
if isempty(epts)
    npts=endpnts(Ds); epts=zeros(sz); epts(npts)=1;
end
%
[xe,ye,ze]=ind2sub(sz,npts)       ; Xe=[xe,ye,ze]  ;
[Xk,kpts]=sort(Ns(npts),'ascend') ; Xe=Xe(kpts,:)  ;
%
Xs=unique(Xk,'sorted'); ns=length(Xs);
nk=.5*(ns-1)*ns; prmt=zeros(nk,2);
for i=1:ns
    for j=i+1:ns
        p=(j-i)+.5*(ns*(ns-1)-(ns-i)*(ns-i+1));
        prmt(p,:)=[i,j];
    end
end
prmt=Xs(prmt); if numel(prmt)==2, prmt=prmt';end
%
Js=regionprops3(Ns,"Centroid","EigenVectors");
Jc=Js.Centroid; Jc=Jc(:,[2,1,3]);
Jq=Js.EigenVectors; Jq=[Jq{:}]; Jq=Jq(:,1:3:end)';
% Co-linear fibers
P=Jq(prmt(:,1),:); Q=Jq(prmt(:,2),:);
R=diff(reshape(Jc(prmt,:)',3,nk,2),1,3)';
P=P./vecnorm(P,2,2); Q=Q./vecnorm(Q,2,2); R=R./vecnorm(R,2,2);
cPQ=dot(P,Q,2); cQR=dot(Q,R,2); cRP=dot(R,P,2);
% Qbol=logical((abs(cPQ)>Qtol).*(abs(cQR)>Qtol).*(abs(cRP)>Qtol));
Qbol=logical(mean(abs([cPQ,cQR,cRP]),2)>Qtol);
%
Ns2=Ns;
if sum(Qbol)>0
    creg=prmt(Qbol,:);
    cg=size(creg,1); lreg={creg(1,:),};lg=1;
    for j=2:cg
        jreg=creg(j,:);jbol=0;
        for k=1:lg
            jkreg=lreg{k}; kbol=any(ismember(jkreg,jreg));
            if kbol
                lreg{k}=unique(cat(2,jkreg,jreg));
            end
            jbol=jbol+kbol;
            if k==lg
                if ~jbol
                    lreg{k+1}=jreg; lg=lg+1;
                end
            end
        end
    end
    %
    for j1=1:lg
        jreg=lreg{j1}; nj=length(jreg); jregk=jreg;
        pregj=nchoosek(jreg,2); pregj=sortrows(sort(pregj,2),[1,2]);
        xregj=vecnorm(diff(reshape(Jc(pregj,:)',3,[],2),1,3))';
        for j2=1:nj-1
            sj2=.5*(nj*(nj-1)-(nj-j2-1)*(nj-j2)); sj1=sj2-nj+j2+1;
            jbol=0; sjs=sj1:sj2;
            while jbol==0 && ~isempty(sjs)
                [~,kj]=min(xregj(sjs)); sjk=sjs(kj);hj=pregj(sjk,:);
                zj1=Xe(Xk==hj(1),:); zj2=Xe(Xk==hj(2),:);
                [~,rj]=min(vecnorm(zj1-permute(zj2,[3,2,1]),2,2),[],'all');
                [rj1,rj2]=ind2sub([size(zj1,1),size(zj2,1)],rj);
                z1=zj1(rj1,:); z2=zj2(rj2,:); zn=z2-z1; nzn=norm(zn);
                Tx=[z1;z2]; jpts=sum((Tx-1).*sz(1).^(0:2),2)+1;
                zbj=[Jq(hj,:);diff(Jc(hj,:));zn]'; zbj=zbj./vecnorm(zbj);
                dzbj= dot(zbj(:,pmtj(:,1)),zbj(:,pmtj(:,2)),1);
                % Qjbol=prod(dzbj>Qtol,2);
                Qjbol=mean(abs(dzbj))>nfac*Qtol;
                Xjbol1=nzn<Xtol1; Xjbol2=nzn<Xtol2;
                jbol=Xjbol1|(Qjbol*Xjbol2); sjs=sjs(sjs~=sjk);
            end
            %
            if jbol>0
                % [Dx,Dy,Dz]=bresenham3(z1,z2);nt=sub2ind(sz,Dx,Dy,Dz);
                nt=impatch(jpts,sz); nt=nt{:};
                fftf=sqrt(sum(ftf(nt))/numel(nt));
                if fftf>=extol
                    jregk(jreg==hj(2))=jregk(jreg==hj(1)); 
                    Ns(nt)=hj(1); Ds(nt)=1;  epts(jpts)=0;
                end
            end
        end
        [idx1,idx2]=ismember(Ns,jreg); Ns2(idx1)=jregk(idx2(idx1));
    end
end
%
Ns=Ns2; [~,ia,ib]=unique(Ns(:),'sorted'); 
js=length(ia); ls=0:js-1; Ns=reshape(ls(ib),sz);
end
%
function npts=endpnts(Ds)
sz=size(Ds); Ns=bwlabeln(Ds,26); ns=max(Ns(:)); npts=[];
ept0=bwmorph3(Ns,"endpoints");
for js=1:ns
    nptj=find((Ns==js).*ept0);
    if isempty(nptj)|numel(nptj)<2
        Njs=find((Ns==js)); [xjs,yjs,zjs]=ind2sub(sz,Njs);
        Xjs=[xjs,yjs,zjs]; Tjs=regiontips(Xjs,nonzeros(Ns==js));
        nptj=sum((Tjs-1).*sz(1).^(0:2),2)+1;
    end
    npts=[npts;nptj];
end
end
%
function Ds=Dskel(ftf,Avf)
sz=size(ftf); Tx=regtips(ftf,sz); 
jpts=sum((Tx-1).*sz(1).^(0:2),2)+1;
nt=impatch(jpts,sz); nt=nt{:};
%
Ds=zeros(sz); Ds(nt)=1;
if Avf<.35
    Ds = bwskel(ftf); Ds = bwmorph3(Ds,'remove')   ;
end
end
%
function [Tx,Lc,npts]=regtips(ftf,sz)
Ls=regionprops3(ftf,"Centroid","EigenVectors","PrincipalAxisLength");
Lx=Ls.PrincipalAxisLength; Lc=Ls.Centroid; Lc=Lc(:,[2,1,3]);
Lq=Ls.EigenVectors; Lq=[Lq{:}];Lq=Lq(:,1:3:end)';
x1=Lc-.5*Lx(:,1).*Lq; x2=Lc+.5*Lx(:,1).*Lq;
[xs,ys,zs]=ind2sub(sz,find(ftf)); Xs=[xs,ys,zs];
Ys=[x1; x2]; idx=knnsearch(Xs,Ys); Tx=Xs(idx,:);
npts=sum((Tx-1).*sz(1).^(0:2),2)+1;
end
%
function f=regiontips(X,T)
Tj=unique(T,'sorted'); nj=length(Tj);
for j=1:nj
    Xj=X(T==Tj(j),:); Yj=fXtip(Xj); 
    idx=knnsearch(Xj,Yj); f(2*j-1:2*j,:)=Xj(idx,:);
end
end
%
function [f,X0]=fXtip(X)
n=size(X,1); X0=mean(X,1); dX=X-X0;
C=(dX'*dX)/(n-1);
[R,D]=svd(C,0); D=diag(D); R2=D(1)/sum(D);
Y=dX*R(:,1); Y1=min(Y);Y2=max(Y);
dY=Y2-Y1;
Xa=(Y1-0.05*dY)*R(:,1)' + X0;
Xb=(Y2+0.05*dY)*R(:,1)' + X0;
f=[Xa;Xb];
end
%
function jpts=impatch(kpts,sz)
kz=sz(1);
[xsz,ysz,zsz]=meshgrid([-1,0,+1],[-kz,0,+kz],[-kz^2,0,kz^2]);
psz=xsz(:)+ysz(:)+zsz(:); ndx=size(kpts,1) ; 
ndx=(numel(kpts)>2)*(ndx-1)+1; jpts=cell(ndx,1);
[Px,Py,Pz]=ind2sub(sz,kpts(:)); PX=[Px,Py,Pz]; 
PX=permute(reshape(PX',3,ndx,[]),[2,1,3]);
dPX=diff(PX,1,3); gPX=vecnorm(dPX,2,2); nPX=dPX./gPX;
%
for jp=1:ndx
    ip=1; jpt=kpts(jp,ip); jpts{jp}(ip)=jpt; jKX=0;
    while jKX<gPX(jp)
        ip=ip+1;
        kpt=jpt+psz; kpt=kpt(logical((kpt>=1).*(kpt<kz^3)));
        kpt=kpt(~ismember(kpt,jpts{jp}));
        [Kx,Ky,Kz]=ind2sub(sz,kpt); KX=[Kx,Ky,Kz];
        dKX=KX-PX(jp,:,1); gKX=vecnorm(dKX,2,2); nKX=dKX./gKX;
        cKX=sum(nKX.*nPX(jp,:),2); [~,nk]=max(cKX);
        jpt=kpt(nk);jpts{jp}(ip)=jpt; jKX=gKX(nk);
    end
end
end
%
function mtf=imsplit(MVf)
fs=fspecial3('prewitt'); mf1=mat2gray(MVf); mf1=exp(-6.*(1-mf1));
mf2=imfilter(MVf,fs,'symmetric'); mf3=mat2gray(mf2);
kj=.3; mfc=sqrt(kj*mf1.^2+(1-kj)*mf3.^2); zv=0.05; df=7;
mf4=fibermetric(mfc,1:df,'ObjectPolarity','bright','StructureSensitivity',zv);
mf5=medfilt3(mf4); mf6=imgaussfilt3(mf5,.35,'FilterSize',[3,3,3]);
mf7=mat2gray(mf4.*(mf6>.01)); 
mf8=fibermetric(mf7,1:df,'ObjectPolarity','bright','StructureSensitivity',zv);
mf9=exp(0.0558-173.815*exp(-8.043*mf8)); mtf=mf9>zv; 
end
%
%
% functions 
% [~,ks]=min(sqrt(sum((Y-permute(X,[3,2,1])).^2,2)),[],3) - knnsearch(X,Y)
% dXf2=permute(vecnorm(Xf-permute(Xe,[3,2,1]),2,2),[1,3,2]) -pdist2(Xf,Xe);
%

