% author: Giuseppe Giacopelli
% pre-print: A full-scale agent-based model of Lombardy COVID-19 dynamics 
% to explore social networks connectivity and vaccine impact on epidemic
% license: GPL-3.0

close all
clear all

N=3*10^5;
cdist=1.0;
Pop=rand(N,2)*40;
temp_idzs = rangesearch(Pop,Pop,cdist);
Nconn=length(cell2mat(temp_idzs'));

PopA=Pop(and(Pop(:,1)<20,Pop(:,2)<20),:);
temp_idzs = rangesearch(PopA,PopA,cdist);
NconnA=length(cell2mat(temp_idzs'));

PopB=Pop(and(Pop(:,1)<20,Pop(:,2)>=20),:);
temp_idzs = rangesearch(PopB,PopB,cdist);
NconnB=length(cell2mat(temp_idzs'));

PopC=Pop(and(Pop(:,1)>=20,Pop(:,2)<20),:);
temp_idzs = rangesearch(PopC,PopC,cdist);
NconnC=length(cell2mat(temp_idzs'));

PopD=Pop(and(Pop(:,1)>=20,Pop(:,2)>=20),:);
temp_idzs = rangesearch(PopD,PopD,cdist);
NconnD=length(cell2mat(temp_idzs'));

NconnT=NconnA+NconnB+NconnC+NconnD;
Acc=NconnT/Nconn;