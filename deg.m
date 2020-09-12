% author: Giuseppe Giacopelli
% pre-print: A full-scale agent-based model of Lombardy COVID-19 dynamics 
% to explore social networks connectivity and vaccine impact on epidemic
% license: GPL-3.0

clear all
close all

%Set current mode
%mode='st'; header='Standard connectivity'; %standard
mode='ld'; header='Lock down connectivity'; %lock down

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

spath=['Sims/degree_',mode];
if ~exist(spath, 'dir')
       mkdir(spath)
end

density=N/Area; 
FPD=6; %6 frames per day
Nit=14*FPD; %3 days

%data from 04/03/2020 to 11/03/2020
%source: http://opendatadpc.maps.arcgis.com/apps/opsdashboard/index.html#/b0c68bce2cce478eaac82fe38d4138b1
Idata=[552,887,1077,1326,1497,1777,2008,2742,3327,4490,4427,5763,6896,7732,9059,10043];
Rdata=[40,73,139,139,250,376,469,524,550,646,896,900,1085,1198,1660,2011];
Ddata=[23,24,38,55,73,98,135,154,267,333,468,617,744,890,966,1218];
days=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
tdays=[0:Nit]/FPD;
Iint=interp1(days,Idata,tdays,'cubic');
Rint=interp1(days,Rdata,tdays,'cubic');
Dint=interp1(days,Ddata,tdays,'cubic');

S=ones(N,1,'logical');
I=zeros(N,1,'logical');
Check=zeros(Nit,1);

if mode=='st'
    cdist=1.0;
    kmd=43; %Ossertvatorio UnipolSai Lombardia: 43 Km per day
end

if mode=='ld'
    cdist=0.1; %lockdown
    kmd=5; %lockdown
end

pi=1/40500; %source: Estimation of Individual Probabilities 
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
Red_Zone=[sub2ind([lyp/4,lxp/4],6,4)];
inf=find(ismember(ivi,Red_Zone));
I0=10^3;
R0=0;
D0=0; 
Ll=I0+R0+D0;
lif=randperm(length(inf),Ll);
tif=inf(lif);
I(tif)=1;
S(tif)=0;
vd=(kmd/40)*wblrnd(6,1.5,N,1);
Init=Pop;
initI=sum(I)

LUT={};
for j=1:Ll
    LUT{tif(j)}=j;
end

figure(10); hold on;
hw=histogram(vd(:),50);
xlabel('drift kmd')
ylabel('frequency')
title('drift kmds histogram')
mvd=mean(vd(:)) %target 5.2


figure(1); hold on;
title('initial condition')
imagesc(Xp,Yp,Np'/CArea)
colorbar
axis equal
dist=zeros(N,1);

Np=histcounts2(Pop(:,1),Pop(:,2),Xedges,Yedges);

path='Sims/';
v=VideoWriter([path,'deg_sp_',num2str(round(10*kmd)),'_',mode,'.avi']);
v.FrameRate=FPD;
open(v);
disp('writing movie...')
h = waitbar(0,'Please wait...');
Conn={};

for j=1:Ll
    Conn{j}=[];
end

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
    TConn={};
    parfor idv=1:lxy2
        idI=find(and(iv2==idv,I==1));
        idS=find(iv2==idv);
        PI=Pop(idI,:); PS=Pop(idS,:);
        temp_idzs = rangesearch(PS,PI,cdist);
        FConn={};
        for j=1:Ll
            FConn{j}=[];
        end
        if length(temp_idzs)>0
            lt=length(temp_idzs);
            for j=1:lt
                Iids=LUT{idI(j)};
                temp_a=idS(temp_idzs{j}');
                FConn{Iids}=unique([FConn{Iids};temp_a]);
            end
        end
        TConn{idv}=FConn;
        WaitMessage.Send;
    end
    WaitMessage.Destroy;
    
    for j=1:Ll
        for idv=1:lxy2
            temp_b=TConn{idv}{j};
            if length(temp_b)>0
                Conn{j}=unique([temp_b;Conn{j}]);
            end
        end
    end
    
    degree=cellfun(@length,Conn);
    
    
    Np=histcounts2(Pop(:,1),Pop(:,2),Xedges,Yedges);
    NI=CMap.*(100*histcounts2(Pop(I==1,1),Pop(I==1,2),Xedges,Yedges)./max(Np,1));
    td=sqrt((Ppast(:,1)-Pop(:,1)).^2+(Ppast(:,2)-Pop(:,2)).^2);
    dist=dist+td;
    disp(['path length: ',num2str(mean(dist))])
    
    if length(hg)>0
        close(hg);
    end
    
    hg=figure('units','normalized','outerposition',[0 0 1 1]);
    clf('reset')
    
    subplot(2,2,1); hold on;
    imagesc(Xp,Yp,Np'/CArea)
    
    if tdays(i+1)>=1
        dstr=strcat(num2str(floor(tdays(i+1))),' March 2020');
    else
        dstr='29 February 2020';
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
    
    subplot(2,2,2); hold on;
    title('Group tracking')
    imagesc(Xp,Yp,log10(NI'))
    axis([0 dims(1) 0 dims(2)]);
    axis equal
    caxis([-3 2])
    colorbar;
    colormap(gca,'pink')
    xlabel('Km')
    ylabel('Km')
    
    sd=25;
    subplot(2,2,3); hold on;
    title('Total degree')
    hd=histogram(degree,sd,'normalization','pdf');
    xlabel('degree')
    ylabel('probability density')
    
    nd=ceil(tdays(i+1));
    subplot(2,2,4); hold on;
    title('Daily degree')
    hd=histogram(degree/nd,sd,'normalization','pdf');
    xlabel('degree')
    ylabel('probability density')
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

function Pend=move(Pin,step,Map,Ms,pxdim)
Pend=Pin+step;
xv=min(max(floor(Pend(:,2)/pxdim),1),Ms(1));
yv=min(max(floor(Pend(:,1)/pxdim),1),Ms(2));
iv=sub2ind(Ms,xv,yv); idxs=(Map(iv)==0);
Pend(idxs,:)=Pend(idxs,:)-2*step(idxs,:);
end
