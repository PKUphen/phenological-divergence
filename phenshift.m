clear all
clc
load 'I:\workspace\phendata.mat'
phendata(phendata(:,6)<1980,:)=[];

ll=phendata(:,2:3);%站点经纬度信息
stlist(:,2:3)=unique(ll,'rows');%站点经纬度列表
stlist(:,1)=1:length(stlist);
st=nan(length(phendata),1);
for i=1:length(stlist)
    st(ll(:,1)==stlist(i,2)&ll(:,2)==stlist(i,3))=stlist(i);
end
splist=unique(phendata(:,4));%物种列表
spcode=[1:length(splist)]';%物种编号
phlist=unique(phendata(:,5));%物候期列表
phcode=[1:length(phlist)]';%物候期编号
spc=phendata{1,4};
ppha=phendata{1,5};
spid=nan(length(phendata),1);
phid=nan(length(phendata),1);
for i=1:length(splist)
    spid(strcmp(spc,splist{i,1})==1,1)=spcode(i);
end
for i=1:length(phlist)
    
    phid(strcmp(ppha,phlist{i,1})==1,1)=phcode(i);
end
datas=phendata(:,1);
dataset=unique(phendata(:,1));
dscode=[1:length(dataset)]';%数据集编号
for i=1:length(dataset)
    dsid(strcmp(datas,dataset{i,1})==1,1)=dscode(i);
end
phenstan=[dsid st spid phid phendata(:,6:7)];
timlist=unique(phenstan(:,1:4),'rows');
cd 'I:\workspace\';
for i=1:length(splist)
    [x0,~]=find(phenstan(:,53)==spcode(i,1));
    tim=unique(phenstan(x0,1:4),'rows');
    tran=phenstan(x0,:);
    for j=1:size(tim,1)
        disp(['Processing the',num2str(i),'species and ',num2str(j),'time series;'])
        [xt,~]=find(tran(:,1)==tim(j,1)&tran(:,2)==tim(j,2)&tran(:,3)==tim(j,3)&tran(:,4)==tim(j,4));
        ptran=tran(xt,5:6);
        [xi,~]=find(timlist(:,1)==tim(j,1)&timlist(:,2)==tim(j,2)&timlist(:,3)==tim(j,3)&timlist(:,4)==tim(j,4));
        sta(xi,1:11)=[i timlist(xi,1:4) length(ptran(:,1)) min(ptran(:,1)) max(ptran(:,1)) nanmean(ptran(:,2)) nanstd(ptran(:,2))];
%         [sl,z,lc,uc]=mktrend(ptran(:,2),ptran(:,1));
%     sta(xi,12:13)=[sl z];
    [b,bint,r,rint,stats]=regress(ptran(:,2),[ptran(:,1) ones(length(ptran(:,1)),1)]);
    sta(xi,14:17)=[b(1) b(2) stats(1) stats(3)];
    end
end
     save('I:\workspace\Result.mat','sta');
%     save('I:\workspace\plantsite.mat','stlist');
%     save('I:\workspace\plantsplist.mat','splist');
%     save('I:\workspace\plantphlist.mat','phlist');
