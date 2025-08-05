clear, clc, % close all force
% check https://www.mathworks.com/help/images/3d-volumetric-image-processing.html
% -------------------------------------------------------------------------
fname='Gray_Value'; buffer=0; 
Tdata = tabularTextDatastore([fname,'.csv'],"TreatAsMissing","n/a","MissingValue",0);
% change "n/a" to "na" for 'cnt_calibrate.csv'
Tdata.SelectedVariableNames='x_MeanGrayValue'; sz=178; na=18;
Tdata=Tdata.readall; GV= Tdata.x_MeanGrayValue;
GV(GV<=0)=na; GV(GV>40)=na; MV=uint8(reshape(GV,sz,sz,sz));  
% % -----------------------------------------------------------------------
v1=min(GV); f2=max(GV);
bnd=[14.55,22.45  ; % ABS
        v1,14.50  ; % Void
     22.50,40.00 ]; % Fiber
for j=1:3, bol(:,j)=(GV>bnd(j,1)).*(GV<=bnd(j,2)); end
Q=reshape(sum((1:3).*bol,2),sz,sz,sz);
ABS_bw  =reshape(bol(:,1),sz,sz,sz);
Void_bw =reshape(bol(:,2),sz,sz,sz);
Fiber_bw=reshape(bol(:,3),sz,sz,sz);
%
bsz=buffer+1:sz-buffer ; Q=Q(bsz,bsz,bsz); ABS_bw=ABS_bw(bsz,bsz,bsz);
Void_bw=Void_bw(bsz,bsz,bsz) ; Fiber_bw=Fiber_bw(bsz,bsz,bsz) ;
sz=sz-2*buffer;
% % -----------------------------------------------------------------------
[xsz,ysz,zsz]=meshgrid([-1,0,+1],[-sz,0,+sz],[-sz^2,0,+sz^2]);
npsz=xsz(:)+ysz(:)+zsz(:); nqsz=[-1,+1,-sz,+sz,-sz^2,+sz^2];
[xsz,ysz,zsz]=meshgrid([-1,0,+1],[-sz-2,0,+sz+2],[-(sz+2)^2,0,+(sz+2)^2]);
psz=xsz(:)+ysz(:)+zsz(:); qsz=[-1,+1,-sz-2,+sz+2,-(sz+2)^2,+(sz+2)^2];
% -------------------------------------------------------------------------
% Q = imsegkmeans3(MV,3,NumAttempts=5,MaxIterations=200,Threshold=1.e-6);
% vs=volshow(Q);  title("Microstructural Features");
% ABS_bw=Q==1; Void_bw=Q==2; Fiber_bw=Q==3;
% for j=1:3
%     Pj=(Q==j); bnd(j,:)=[min(MV(Pj),[],'all'),max(MV(Pj),[],'all')];
% end
% -------------------------------------------------------------------------
vf=sum(Void_bw,'all')/sz^3*100; ff=sum(Fiber_bw,'all')/sz^3*100;
fprintf('Void  Volume Fraction = %.4f %%\n',vf);
fprintf('Fiber Volume Fraction = %.4f %%\n',ff);
% fnms={'ABS','Voids','Fibers'};
% for j=1:3
%     f1=uifigure('Name',fnms{j});  Qj=bwareaopen(Q==j,15,26);
%     Qj=medfilt3(Qj) ; vs=volshow(Qj,'Parent',f1);
%     sc=vs.Parent; sc.BackgroundColor='w'; sc.BackgroundGradient='off';
%     sc.Box='on'; vs.RenderingStyle='GradientOpacity'; vs.Colormap=gray;
%     vs.GradientOpacityValue=.1+.15*(j-1); vs.Alphamap=1.-.15*j;
% end
%
Voids  =bwlabeln(Void_bw ,26); Fibers=bwlabeln(Fiber_bw,26); 
pVoids =padarray(Voids ,[1,1,1],0,'both');
pFibers=padarray(Fibers,[1,1,1],0,'both');
nv=max(Voids,[],'all'); nf=max(Fibers,[],'all');
freqv=histc(Voids(:),1:nv); freqf=histc(Fibers(:),1:nf);
% --------------------- Isolated Void Regions -----------------------------
% % for i=1:nv
% %     voidi=(Voids==i);
% %     parfor j=1:nf
% %         fiberj=(Fibers==j);
% %         reg(i,j)=max(bwlabeln(voidi|fiberj,26),[],"all")==1;
% %     end 
% % end
% % [ri,cj]=find(reg); conn_reg=[ri,cj];
% -------------------------------------------------------------------------
% parfor i=1:nv
%     voidi=(Voids==i);
%     reg(i)=max(bwlabeln(voidi|Fiber_bw,26),[],"all")>nf;
% end
% iso_reg=find(reg); nvo=length(iso_reg);unconn_vd=ismember(Voids(:),iso_reg);
% unconn_v=sum(unconn_vd,'all')/sum(freqv)*100;
% ------------------------ Alternatively ----------------------------------
fbr=find(pFibers); sfbr=fbr + psz';efbr=sfbr(pFibers(sfbr)==0);
ivdf=pVoids(efbr); conn_idv=unique(ivdf(ivdf>0));
unconn_v=(1-sum(freqv(conn_idv))/sum(freqv))*100;
lstv=(1:nv)'; iso_reg=lstv(~ismember(lstv,conn_idv));
% %--------- Percentage Vol. & Average Dia. of Isolated Voids --------------
unconn_bol=ismember(1:nv,iso_reg);
unconn_frac=sum(freqv(unconn_bol))/sum(freqv)*100;
fprintf('Isolated Void = %.4f %%\n',unconn_frac);
vprops=regionprops3(Voids,"Volume","EquivDiameter","SurfaceArea");
v_vol=vprops.Volume; v_dia=vprops.EquivDiameter; v_surf=vprops.SurfaceArea;
nc_davg=mean(v_dia(unconn_bol)); cn_davg=mean(v_dia(~unconn_bol));
% uc_davg=(6/pi*mean(freqv( unconn_bol)))^(1/3); 
% cn_davg=(6/pi*mean(freqv(~unconn_bol)))^(1/3); 
fprintf(' Iso. void avg. pixel diameter = %.4f pix\n',nc_davg);
fprintf('Conn. void avg. pixel diameter = %.4f pix\n',cn_davg);
% %-------------------------------------------------------------------------
% parfor j=1:nv
%     vndx=find(pVoids==j); vndx=vndx+qsz;
%     vlb=vndx(logical((vndx(:)>=1).*(vndx(:)<=(sz+2)^3)));
%     vlb=vlb(pVoids(vlb)==0); vsrf{j}=vlb(:);
% end
% csrf=arrayfun(@(x) numel(x{:}),vsrf);
% %
% v_sph=(36*pi*v_vol.^2).^(1/3)./csrf';
% nc_savg=mean(v_sph(unconn_bol)); cn_savg=mean(v_sph(~unconn_bol));
% fprintf(' Iso. void avg. Wadell sphericity = %.4f \n',nc_savg);
% fprintf('Conn. void avg. Wadell sphericity = %.4f \n',cn_savg);
% ----------------------- Tip Voids & interacting Tips --------------------
[Fskel,Epts,Tips]=fiberfilter_2(MV,Fibers,fname);  
%load(['.\Data\',fname,'.mat'],'I','E','T'); Fskel=I; Epts=E; Tips=T;
Fskel=Fskel(bsz,bsz,bsz); Epts=Epts(bsz,bsz,bsz); Tips=Tips(bsz,bsz,bsz);
% zero Epts Boundary Buffer zone 
Etol=5; gsz=reshape(1:sz^3,sz,sz,sz); Erng=Etol+1:sz-Etol;
Egsz=gsz(Erng,Erng,Erng); Epts(~ismember(gsz,Egsz))=0;
%
pEpts=padarray(Epts,[1,1,1],0,'both'); pfte=find(pEpts); pfjz=pfte + psz';
ifjz=pfjz(pFibers(pfjz)==0); mdfv=ifjz(pVoids(ifjz)>0);
ndfv=unique(pVoids(mdfv)); conn_vtip=sum(freqv(ndfv))/sum(freqv)*100;
fprintf('Tip Voids = %.4f %%\n',conn_vtip);
tipfrac=sum(Epts(:)>0)/sum(Fiber_bw(:))*100;
fprintf('Frac of Fiber Vol. considered Tip Regions = %.4f %%\n',tipfrac);
%
% vnms={'Tip Voids','Isolated void'}; vval={ndfv,iso_reg};
% for j=1:2
%     % ival{j}=lstv(~ismember(lstv,vval{j}));
%     f1=uifigure('Name',vnms{j});  Rj=ismember(Voids,vval{j});
%     Rj=medfilt3(Rj) ; vs=volshow(Rj,'Parent',f1);
%     sc=vs.Parent; sc.BackgroundColor='w'; sc.BackgroundGradient='off';
%     sc.Box='on'; vs.RenderingStyle='GradientOpacity'; vs.Colormap=gray;
%     vs.GradientOpacityValue=.15; vs.Alphamap=.5 ;
% end
%
vfte=pfte(any(ismember(pfjz,mdfv),2)); 
jpte=unique(pEpts(vfte)); kEpts=ismember(pEpts,jpte);
nEpts=pEpts.*kEpts;cEpts=nEpts(2:sz+1,2:sz+1,2:sz+1);
fcnts=(cEpts>0).*(Tips>0).*Fskel; 
fpte=unique(fcnts(:)); fpte=fpte(fpte>0); wEpts=ismember(Fskel,fpte); 
% vfrac1=sum(kEpts,'all')/sum(pEpts>0,'all');
% vfrac2=length(jpte)/max(pEpts,[],'all')   ; 
% vfrac3=sum(wEpts,'all')/sum(Fskel>0,'all') ;
% vfrac4=length(fpte)/max(Fskel,[],'all') ;
% vfrac5=sum(fcnts>0,'all')/sum(Tips>0,'all');
% fprintf('Fiber w/ Tip Voids = %.4f %%\n',vfrac1);
%
% kj=355; Qf=Fibers==kj; Rf=zeros(sz,sz,sz); Rf(Qf&(E>0))=1;
% pFj=pFibers==kj; cfte=find(pEpts.*pFj); cfjz=cfte + psz';
% sfjz=cfjz(pFibers(cfjz)==0); sdfv=sfjz(pVoids(sfjz)>0);
% udfv=unique(pVoids(sdfv)); Rf(ismember(Voids,udfv))=2;
% [vs,sc]=VolPlot(Rf,Qf,[],.40,.30,.70,[],{'gray','hsv'});
% % ------------------------- Visulaization ---------------------------------
R=zeros(sz,sz,sz); R(Epts>0)=0; R(ismember(Voids,ndfv))=2; 
R(Void_bw & (R~=2))=3;
[vs,sc]=VolPlot(R,Q-(Q==1),[],.40,.30,.70,[],{'gray','hsv'});
sc.CameraPosition=[375, 375, 300];
%
lf=50;Lf=ismember(Fibers,find(freqf<lf));
Q2=Q; Q2(Lf)=1; Q2(Q2==1)=0; % Q2(R==2)=0;
R2=R; R2(Lf)=0; % R2(R==2)=0;
ndx=zoomView(3,[2,2,2],sz,Q2,R2);
%
P=zeros(sz,sz,sz);P(ndx(1,:),ndx(2,:),ndx(3,:))=1; np=find(P); 
nps= np+npsz'; rp=find(sum(1-P(nps),2)>9); ep=np(rp); 
%
Rk=R; Rk((Q==3)&(~Epts))=4; 
kz=round(sz/3); Qk=Q; Qk(kz+1:sz,kz+1:sz,kz+1:sz)=0; Qk(ep)=3;
vsk=VolPlot(Rk,Qk,P,.65,.45,.15,[2,.25],{'bone','hsv'});
% % ---------------------------- Slice Planes -------------------------------
% slc={'YZ', 'XZ', 'XY'}; Ax={'X', 'Y', 'Z'};
% xk={125,125,125}; Rs=R; Rs(Q==3)=4;
% for j=1:3
%         figure('Name',...
%         sprintf('%s - Slice @ %d voxel unit on %s Axis',slc{j},xk{j},Ax{j}))
%         MV_slice(MV,Rs, xk{j},j); 
% end
% ----------------------------------------------------------------------
% volf=regionprops3(Fibers,'Solidity','PrincipalAxisLength');
% sdty=volf.Solidity; [fracf,idf]=sort(sdty,'descend');
% vdf=idf(fracf>.45); Lf=volf.PrincipalAxisLength; lenf=Lf(vdf); 
% ivdf=vdf(lenf>50); vFbrs=ismember(Fibers,ivdf); vEpts=(E>0).*vFbrs;
% %
% vsf=VolPlot(vEpts,vFbrs,[],.35,.85,.35,[],{'gray','hsv'});
% vsf.OverlayColormap=[34,139,34]/255;  
% %
% Z=zeros(sz,sz,sz); vEpts=(Epts>0).*vFbrs; Z(vEpts>0)=1;
% sEpts=padarray(vEpts,[1,1,1],0,'both'); sfte=find(sEpts); sfjz=sfte + psz';
% jfjz=sfjz(pFibers(sfjz)==0);hdfv=jfjz(pVoids(jfjz)>0);udfv=unique(pVoids(hdfv)); 
% svol=6000; wdfv=udfv(freqv(udfv)<svol); Z(ismember(Voids,wdfv))=2; 
% %
% adfv=lstv(~ismember(lstv,ndfv)); xsfv=randi([1,length(adfv)],700,1); 
% xdfv=adfv(xsfv); pdfv=xdfv(freqv(xdfv)>7); Z(ismember(Voids,pdfv))=3;
% VolPlot(Z,vFbrs,[],.35,.35,.25,[],{'gray','hsv'});
% ----------------------------------------------------------------------
dLx=10; [Mskel,cnt]=bwskelseg(Fskel,dLx); freqm=histc(Mskel(:),1:cnt);
fprops=regionprops3(Mskel,"EigenVectors","PrincipalAxisLength");
fx=fprops.PrincipalAxisLength; Lx=fx(:,1);
fQ=fprops.EigenVectors; fQ=[fQ{:}]; fQ=fQ(:,1:3:end)';
avg_Qf=mean(fQ(freqm>=dLx,:).^2,1);
wavg_Qf=sum(freqm(freqm>=dLx).*fQ(freqm>=dLx,:).^2,1)/sum(freqm(freqm>=dLx));
fprintf(['Average Fiber Orientation \n'...
         'Axx = %4.2f, Ayy=%4.2f, Azz=%4.2f \n'],avg_Qf);
% -------------------------------------------------------------------------
function C=MV_slice(MV,R,k,flg)
Q=R==2; [n1,n2,n3]=size(MV); m={1:n1,1:n2,1:n3};
ord=1:3; ord([flg,3])=[3,flg]; m{flg}=k;
C=labeloverlay(permute(MV(m{:}),ord),permute(R(m{:}),ord));
s1=subplot(1,3,1); image(s1,permute(MV(m{:}),ord),'AlphaData',0.85);
colormap(s1,"hsv"); axis equal off
s2=subplot(1,3,2); image(s2,C,'AlphaData',0.85);
colormap(s2,jet) ; axis equal off
s3=subplot(1,3,3); imshow(permute(Q(m{:}),ord));
set(gcf,'Color','w');
end
%
function [vs,sc]=VolPlot(R,Q,P,aQ,aR,gv,gs,cmap)
fs=uifigure('Name','All Regions'); vs=volshow(Q,'Parent',fs);
sc=vs.Parent; sc.Lighting='off'; sc.BackgroundColor='w'; sc.Box='on';
sc.BackgroundGradient='off'; vs.Colormap=colormap(fs,cmap{1});
vs.RenderingStyle='GradientOpacity'; vs.GradientOpacityValue=gv;
if isempty(P)
    vs.Alphamap=aQ; P=1; style='Label';
else
    vs.AlphaData=(aQ*(P-1.)+1.); style='Volume';
end
vs.OverlayData=R.*P; vs.OverlayColormap=colormap(fs,cmap{2});  
vs.OverlayRenderingStyle=[style,'Overlay']; vs.OverlayAlphamap=aR;
vs.OverlayThreshold=1.e-5;
if ~isempty(gs)
    sc.Denoising='on'; sc.DenoisingSigma =gs(1);
    sc.DenoisingDegreeOfSmoothing=gs(2);
end
end
% -------------------------------------------------------------------------
function ndx=zoomView(gsz,gdx,sz,Q,R)
% gsz=5; gdx=[3,4,2];
jsz=round(sz/gsz); jdx=jsz*((1:gsz)'-1)+(1:jsz); ndx=jdx(gdx,:);
jx=ndx(1,:); jy=ndx(2,:); jz=ndx(3,:); Qj=Q(jx,jy,jz); Rj=R(jx,jy,jz);
fs=uifigure('Name','Sub Region'); vs=volshow(Qj,'Parent',fs);
sc=vs.Parent; sc.Lighting='off';sc.BackgroundColor='w';
sc.Box='on'; sc.BackgroundGradient='off';  vs.OverlayData=Rj; 
vs.RenderingStyle='GradientOpacity';  vs.Colormap=gray;
vs.OverlayColormap=hsv; vs.OverlayAlphamap=.35;
vs.OverlayThreshold=1.e-5;vs.Alphamap=.45; vs.GradientOpacityValue=.75;
sc.Denoising='on';sc.DenoisingSigma =2;sc.DenoisingDegreeOfSmoothing=.25;
end
% -------------------------------------------------------------------------
function [Mskel,cnt]=bwskelseg(Nskel,Lx)
ns=max(Nskel(:)); sz=size(Nskel); Mskel=zeros(sz);
cnt=0;
for js=1:ns
    fjs=Nskel==js; cjs=find(fjs); Vx=sum(fjs,'all');
    if Vx<Lx
        Mskel(cjs)=cnt+1; cnt=cnt+1;
    else
        kjs=round(Vx/Lx); [xs,ys,zs]=ind2sub(sz,cjs); Xs=[xs,ys,zs];
        %
        Ts=spectralcluster(Xs,kjs,'Distance','euclidean','Radius',inf,...
            'SimilarityGraph','epsilon','KNNGraphType','complete',...
            'LaplacianNormalization','symmetric','ClusterMethod','kmeans');
        %
        Mskel(cjs)=cnt+Ts; cnt=cnt+max(Ts(:));
    end
    fprintf('%d of %d\n',js,ns);
end
end