% author: Giuseppe Giacopelli
% pre-print: A full-scale agent-based model of Lombardy COVID-19 dynamics 
% to explore social networks connectivity and vaccine impact on epidemic
% license: GPL-3.0

clear all
close all

%Set current mode
mode='fl'; header='Descent fitting'; %falling
%mode='vc'; header='Vaccine scenario'; %vaccine

vpd=10^6; %vaccines per day

spath=['Sims/sim_f_',mode];
if ~exist(spath, 'dir')
       mkdir(spath)
end

Map=(rgb2gray(imread('lombardia.png'))<100);
Ms=size(Map);
Map(1,:)=0; Map(end,:)=0; Map(:,1)=0; Map(:,end)=0;
Den=imread('density.png');
se = offsetstrel('ball',3,3);
Dc=zeros(size(Den,1),size(Den,2),size(Den,3));
Dc(:,:,1)=imdilate(Den(:,:,1),se);
Dc(:,:,2)=imdilate(Den(:,:,2),se);
Dc(:,:,3)=imdilate(Den(:,:,3),se);
Dmap=imresize(rgb2gray(uint8(Dc)),Ms);
Dvals=2000*((255-(double(Dmap)))/255)+50*Map;
N=10.06*10^6;
Dvals=flipud(round(N*Dvals/sum(Dvals(:))));

Map=flipud(Map);
Area_Lomb=23844; %km^2
pxdim=sqrt(Area_Lomb/sum(Map(:)))

rng(1)
dims=Ms*pxdim;
Area=dims(1)*dims(2);

N=sum(Dvals(:));
[ri,ci]=find(Dvals>0);
Pop=zeros(N,2);

li=length(ri(:));
cc=0;
WaitMessage = parfor_wait(li,'Waitbar', true);
for j=1:li
    cit=ci(j); rit=ri(j);
    Npt=Dvals(rit,cit);
    Pop((cc+1):(cc+Npt),:)=rand(Npt,2)*pxdim+[cit,rit]*pxdim;
    cc=cc+Npt;
    WaitMessage.Send;
end
WaitMessage.Destroy;


density=N/Area; 
FPD=6; %6 frames per day
Nit=14*FPD; %3 days

%data from 31/05/2020 to 14/06/2020
%source: http://opendatadpc.maps.arcgis.com/apps/opsdashboard/index.html#/b0c68bce2cce478eaac82fe38d4138b1
Idata=[20996,20861,20255,20224,20224,19853,19499,19420,19319,18297,17857,17340,17024,16875,15989,15976];
Rdata=[51860,52026,52807,53046,53101,53853,54322,54505,54768,55967,56474,57218,57775,58201,59220,59484];
Ddata=[16112,16131,16143,16172,16201,16222,16249,16270,16302,16317,16349,16374,16405,16428,16449,16457];
DNC=[449,478,513,533,540,552,610,633,653,682,701,752,783,868,910,928,957,1035,1093,1172,1126,1194,1278,1394,1486,1552,1645,1789,1864];
days=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
tdays=[0:Nit]/FPD;
edays=5*7;
eframes=edays*FPD;
vpframe=round(vpd/FPD);
%phist=fliplr(kron(DNC(1:edays)/sum(DNC(1:edays)),ones(1,FPD)/FPD));
phist=ones(1,eframes)/eframes;
dels=fliplr([0:(eframes-1)]/FPD);

Iint=interp1(days,Idata,tdays,'cubic');
Rint=interp1(days,Rdata,tdays,'cubic');
Dint=interp1(days,Ddata,tdays,'cubic');

figure(11); hold on;
sgtitle('Original data');

subplot(2,2,1); hold on;
plot(tdays,Iint,'--r');
xlabel('days')
ylabel('Infected')

subplot(2,2,2); hold on;
plot(tdays,Rint,'--g');
xlabel('days')
ylabel('Recovered')

subplot(2,2,3); hold on;
plot(tdays,Dint,'--k');
xlabel('days')
ylabel('Deaths')

subplot(2,2,4); hold on;
RRatio=Rdata./Idata;
DRatio=Ddata./Idata;
KRatio=Ddata./Rdata;
pdeath=0.05; %source: https://www.worldometers.info/coronavirus/country/italy/
plot(days,KRatio,'-xr');
plot([days(1),days(end)],pdeath*[1 1],'--r');

S=ones(N,1,'logical');
I=zeros(N,1,'logical');
R=zeros(N,1,'logical');
D=zeros(N,1,'logical');
dT=zeros(N,1,'double');
Check=zeros(Nit,1);

if or(or(mode=='ld',mode=='fl'),mode=='vc')
    cdist=0.3; %lock down
    kmd=15; %lock down
else
    cdist=1.0; %canonical case
    kmd=43; %Ossertvatorio UnipolSai Lombardia: 43 Km per day
end

if or(or(mode=='2m',mode=='fl'),mode=='vc') 
    mult=0.5; %source: Physical distancing, face masks, and eye protection to
    %prevent person-to-person transmission of SARS-CoV-2 and
    %COVID-19: a systematic review and meta-analysis
else
    mult=1; %canonical case
end

pi=mult*1/40500; %source: Estimation of Individual Probabilities
% of COVID-19 Infection, Hospitalization, and
% Death From A County-level Contact of Unknown infection Status
Sn=zeros(1,Nit);
In=zeros(1,Nit);
ncell=5;
Cdim=ncell*pxdim;
CArea=Cdim^2;
Xedges=[0:Cdim:dims(1)];
Yedges=[0:Cdim:dims(2)];
Xp=(Xedges(1:end-1)+Xedges(2:end))/2;
Yp=(Yedges(1:end-1)+Yedges(2:end))/2;
lxp=length(Xp(:));
lyp=length(Yp(:));
CMap=imresize(Map',[lxp,lyp]);
Np=histcounts2(Pop(:,1),Pop(:,2),Xedges,Yedges);
xvi=min(max(floor(Pop(:,1)/(4*Cdim)),1),lxp/4);
yvi=min(max(floor(Pop(:,2)/(4*Cdim)),1),lyp/4);
ivi=sub2ind([lyp/4,lxp/4],yvi,xvi);
%Red_Zone=[sub2ind([lyp/4,lxp/4],6,4)];
Red_Zone=unique(ivi);
inf=find(ismember(ivi,Red_Zone)); li=length(inf);
I0=Idata(1);
R0=Rdata(1);
D0=Ddata(1);
lif=randperm(length(inf),I0+R0+D0);
Igid=inf(lif(1:I0)); lI=length(Igid);
dtid=randperm(lI,lI);
I(Igid)=1;
rhist=round(phist*lI);
dT(Igid(dtid(1:(lI-sum(rhist(2:end))))))=dels(1);
cc=lI-sum(rhist(2:end));
for i=2:eframes
    dT(Igid(dtid(cc+1:cc+rhist(i))))=dels(i);
    cc=cc+rhist(i);
end

figure(31); hold on;
ht=histogram(dT(dT>0),25);

R(inf(lif((I0+1):(R0+I0))))=1;
D(inf(lif((R0+I0+1):(D0+R0+I0))))=1;
S(inf(lif))=0;
vd=(kmd/40)*wblrnd(6,1.5,N,1);
Init=Pop;
initI=sum(I)

if mode=='vc'
    idS=find(S==1); lS=length(idS);
    vperm=randperm(lS,round(lS*0.7));
    idvc=idS(vperm);
    S(idvc)=0;
    R(idvc)=1;
end
    

figure(10); hold on;
hw=histogram(vd(:),50);
xlabel('drift kmd')
ylabel('frequency')
title('drift kmds histogram')
mvd=mean(vd(:)) %target 5.2

i=0;
figure(1); hold on;
title('initial condition')
imagesc(Xp,Yp,Np'/CArea)
colorbar
axis equal
dist=zeros(N,1);
nS=zeros(Nit+1,1);
nI=zeros(Nit+1,1);
nR=zeros(Nit+1,1);
nD=zeros(Nit+1,1);

nS(1)=sum(S); nI(1)=sum(I);
nR(1)=sum(R); nD(1)=sum(D);

Np=histcounts2(Pop(:,1),Pop(:,2),Xedges,Yedges);
NI=CMap.*(100*histcounts2(Pop(I==1,1),Pop(I==1,2),Xedges,Yedges)./max(Np,1));
NS=CMap.*(100*histcounts2(Pop(S==1,1),Pop(S==1,2),Xedges,Yedges)./max(Np,1));
figure(2); hold off;

subplot(2,3,1); hold on;
imagesc(Xp,Yp,Np'/CArea)
hold on; sgtitle(['frame ',num2str(i)]);
title('density')
axis([0 dims(1) 0 dims(2)]);
caxis([0 1200])
colorbar;
colormap(gca,'parula')

subplot(2,3,2); hold on;
title('susceptible fraction')
imagesc(Xp,Yp,log10(NS'))
axis([0 dims(1) 0 dims(2)]);
caxis([-3 2])
colorbar;
colormap(gca,'winter')

subplot(2,3,3); hold on;
title('infected fraction')
imagesc(Xp,Yp,log10(NI'))
axis([0 dims(1) 0 dims(2)]);
caxis([-3 2])
colorbar;
colormap(gca,'hot')

subplot(2,3,4); hold on;
title('Infected'); hold on;
plot([0]/FPD,nI(1),'r-x');

subplot(2,3,5); hold on;
title('Recovered'); hold on;
plot([0]/FPD,nR(1),'g-x');

subplot(2,3,6); hold on;
title('Deaths'); hold on;
plot([0]/FPD,nD(1),'k-x');


path='Sims/';
v=VideoWriter([path,'sim_f_sp_',num2str(round(10*kmd)),'_',mode,'.avi']);
v.FrameRate=FPD;
open(v);
disp('writing movie...')
h = waitbar(0,'Please wait...');
hI=nI(1);
hg=[];
for i=1:Nit
    [X,Y]=meshgrid(Xp,Yp);
    if mod(Nit,FPD)==0
        Pop=Init;
    end
    
    Ppast=Pop;
    lXY=length(X(:));
    Px=[X(:),Y(:)];
    Rd=squareform(pdist(Px));
    Rd(Rd(:)==0)=1;
    Gx=0.5*mean(Np(:).*(Px(:,1)-Px(:,1)')./Rd.^3,1);
    Gy=0.5*mean(Np(:).*(Px(:,2)-Px(:,2)')./Rd.^3,1);
    xv=min(max(floor(Pop(:,1)/Cdim),1),lxp);
    yv=min(max(floor(Pop(:,2)/Cdim),1),lyp);
    iv=sub2ind([lyp,lxp],yv,xv);
    Pop(:,1)=Pop(:,1)+(kmd/40)*Gx(iv)';
    Pop(:,2)=Pop(:,2)+(kmd/40)*Gy(iv)';
    Pop=move(Pop,vd.*randn(N,2),Map,Ms,pxdim);
    lxy=lxp*lyp;
    dcells=4;
    samples_per_area=dcells^2*N/lxy;
    xv2=min(max(floor(Pop(:,1)/(dcells*Cdim)),1),lxp/dcells);
    yv2=min(max(floor(Pop(:,2)/(dcells*Cdim)),1),lyp/dcells);
    iv2=sub2ind([lyp/dcells,lxp/dcells],yv2,xv2);
    lxy2=lxp*lyp/dcells^2;
    
    idzs={};
    WaitMessage = parfor_wait(lxy2,'Waitbar', true);
    parfor idv=1:lxy2
        idI=find(and(iv2==idv,I==1));
        idS=find(and(iv2==idv,S==1));
        PI=Pop(idI,:); PS=Pop(idS,:);
        temp_idzs = rangesearch(PS,PI,cdist);
        if length(temp_idzs)>0
            idT=find(~cellfun(@isempty,temp_idzs));
            nidzs=temp_idzs(idT);
            temp=idS(cell2mat(nidzs'));
            idzs(idv)={temp(:)};
        end
        WaitMessage.Send;
    end
    WaitMessage.Destroy;
    
    idT=find(~cellfun(@isempty,idzs));
    idz=cell2mat(idzs(idT)');
    
    lz=length(idz);
    rpz=unique(randperm(lz,round(lz*pi)));
    lrz=length(rpz);
    S(idz(rpz))=0;
    I(idz(rpz))=1;
    
    idI=find(I==1); lI=length(idI);
    dT(idI)=dT(idI)+1/FPD;
    idT=find(and(dT>=edays,I==1)); lT=length(idT);
    rT=randperm(lT,lT); npD=round(pdeath*lT); npR=lT-npD;
    I(idT)=0;
    D(idT(rT(1:npD)))=1;
    R(idT(rT(npD+1:end)))=1;
    
%     if mode=='vc'
%         idS=find(S==1); lS=length(idS);
%         vperm=randperm(lS,vpframe);
%         idvc=idS(vperm);
%         S(idvc)=0;
%         R(idvc)=1;
%     end
    
    nS(i+1)=sum(S); nI(i+1)=sum(I); nD(i+1)=sum(D); nR(i+1)=sum(R);
    Check(i)= nS(i+1)+nI(i+1)+nD(i+1)+nR(i+1);
    
    Np=histcounts2(Pop(:,1),Pop(:,2),Xedges,Yedges);
    NI=CMap.*(100*histcounts2(Pop(I==1,1),Pop(I==1,2),Xedges,Yedges)./max(Np,1));
    NS=CMap.*(100*histcounts2(Pop(S==1,1),Pop(S==1,2),Xedges,Yedges)./max(Np,1));
    
    td=sqrt((Ppast(:,1)-Pop(:,1)).^2+(Ppast(:,2)-Pop(:,2)).^2);
    dist=dist+td;
    disp(['path length: ',num2str(mean(dist))])
    
    if length(hg)>0
        close(hg);
    end
    
    hg=figure('units','normalized','outerposition',[0 0 1 1]);
    clf('reset')
    subplot(2,4,[1 2]); hold on;
    imagesc(Xp,Yp,Np'/CArea)
    
    if tdays(i+1)>=1
        dstr=strcat(num2str(floor(tdays(i+1))),' June 2020');
    else
        dstr='31 May 2020';
    end
    
    hstr=strcat(num2str((24/FPD)*floor(mod(i,FPD))),':00 CEST');
    
    tl=[header,': ',hstr,' ',dstr,' (Population ',num2str(N),')'];
    hold on; sgtitle(tl);
    title('Density')
    axis([0 dims(1) 0 dims(2)]);
    axis equal
    caxis([0 1200])
    colorbar;
    colormap(gca,'parula')
    xlabel('Km')
    ylabel('Km')
    
    subplot(2,4,[3 4]); hold on;
    title('Infected tracking')
    imagesc(Xp,Yp,log10(NI'))
    caxis([-3 2])
    axis([0 dims(1) 0 dims(2)]);
    axis equal
    colorbar;
    colormap(gca,'hot')
    xlabel('Km')
    ylabel('Km')
    
    sf=1/1000; 
    slabel='thousands';
    
    subplot(2,4,5); hold on;
    title('Infected'); hold on;
    plot([0:i]/FPD,sf*nI(1:(i+1)),'r-');
    plot(tdays(1:i+1),sf*Iint(1:i+1),'--r');
    ylim_curr = get(gca,'ylim'); ylim([0 ylim_curr(2)+5]);
    ylabel(slabel)
    xlabel('days')
    
    subplot(2,4,6); hold on;
    title('Recovered'); hold on;
    plot([0:i]/FPD,sf*nR(1:(i+1)),'g-');
    plot(tdays(1:i+1),sf*Rint(1:i+1),'--g');
    ylim_curr = get(gca,'ylim'); ylim([0 ylim_curr(2)+5]);
    ylabel(slabel)
    xlabel('days')
    
    subplot(2,4,7); hold on;
    title('Deaths'); hold on;
    plot([0:i]/FPD,sf*nD(1:(i+1)),'k-');
    plot(tdays(1:i+1),sf*Dint(1:i+1),'--k');
    ylim_curr = get(gca,'ylim'); ylim([0 ylim_curr(2)+5]);
    ylabel(slabel)
    xlabel('days')
    
    subplot(2,4,8); hold on;
    title('Recover Ratio'); hold on;
    plot([0:i]/FPD,nR(1:(i+1))./(nR(1:(i+1))+nD(1:(i+1))),'-b');
    plot(tdays(1:i+1),Rint(1:(i+1))./(Rint(1:(i+1))+Dint(1:(i+1))),'--b');
    ylim([0 1])
    
    if mod(i,FPD)==0
        saveas(gcf,[spath,'/day_',num2str(i/FPD),'.png'])
    end
    
    FrameOutput=getframe(gcf);
    writeVideo(v,FrameOutput);
    waitbar(i/Nit,h)
end
close(h)
close(v);
clear v
disp('movie saved')

disp(['path length per day: ',num2str(mean(dist)*FPD/Nit)])
disp(['Integrity: ',num2str(round(100*(sum(Check==N)/Nit))),'%']);

function Pend=move(Pin,step,Map,Ms,pxdim)
Pend=Pin+step;
xv=min(max(floor(Pend(:,2)/pxdim),1),Ms(1));
yv=min(max(floor(Pend(:,1)/pxdim),1),Ms(2));
iv=sub2ind(Ms,xv,yv); idxs=(Map(iv)==0);
Pend(idxs,:)=Pend(idxs,:)-2*step(idxs,:);
end
