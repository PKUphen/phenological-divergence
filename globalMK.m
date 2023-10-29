clear all
clc
%导入animal数据
load 'I:\workspace\animalphennumdata.mat';
[phc,phli]=xlsread('I:\workspace\动物物候期分类.xlsx','Sheet1','J2:O448');
[spc,spli]=xlsread('I:\workspace\动物物候期分类.xlsx','Sheet1','F2:I951');

%导入plant数据
% load 'I:\workspace\plantphennumdata.mat';
% phenstan=unique(phenstan,'rows');
% [spc,spli]=xlsread('I:\workspace\植物物种_物候期分类.xlsx','植物','A2:H1651');%物种列表和编码
% [phc,phli]=xlsread('I:\workspace\植物物种_物候期分类.xlsx','植物','I2:O254');%物候期列表和编码



for i=1:length(spc)
    phenstan(phenstan(:,3)==spc(i,1),7)=spc(i,4);
end

for i=1:length(phc)
    phenstan(phenstan(:,4)==phc(i,1),8)=phc(i,6);
end
% varS=0;
phcl=[1 1 1 1 1 2 2 2 2 3 3 3 3 4 5 5 6 6 6 6;2 4 6 7 10 2 6 7 10 2 6 7 10 2 2 6 2 3 6 10]';%animal物候期编码
% phcl=[2 3 4 5 6 7]';%plant物候期编码
for i=1:length(phcl)
    datatran=phenstan(phenstan(:,7)==phcl(i,1)&phenstan(:,8)==phcl(i,2),2:8);%animal
%     datatran=phenstan(phenstan(:,8)==phcl(i,1),2:8);%plant
    tims=unique(datatran(:,1:3),'rows');
    s=nan(size(tims,1),2500);
    nn=0;
    for j=1:size(tims,1)
        disp(['Processing the ',num2str(j),'/',num2str(length(tims)),'time series;'])
        tran=datatran(datatran(:,1)==tims(j,1)&datatran(:,2)==tims(j,2)&datatran(:,3)==tims(j,3),4:5);
        
        ti=1;
        n=length(tran);
        if n>=5
nn=nn+n;
        for k=1:n-1
            for ki=k+1:n
                s(j,ti) = ( tran(ki,2) - tran(k,2) ) / ( tran(ki,1) - tran(k,1) );
                ti = ti + 1;
            end
        end
        end
    end

nc=length(s(~isnan(s)));
nn=floor((1+sqrt(1+8*nc))/2);

ss=sum(sign(s(~isnan(s))));


v = (( nn * ( nn - 1 ) * ( 2 * nn + 5 ) )) / 18;
if ss == 0
z = 0;
elseif ss > 0
z = ( ss - 1 ) / sqrt( v );
else
z = ( ss + 1 ) / sqrt( v );
end
nor = 1.96;%z的绝对值在大于1.64、 1.96、 2.58时，分别表示通过了置信度90%、95%和99%的显著性检验
m1 = fix( ( nc - nor * sqrt( v ) ) / 2 );
m2 = fix( ( nc  + nor * sqrt( v ) ) / 2 );
s1 = sort( s (~isnan(s)));
lc = s1( m1 );
uc = s1( m2 + 1 );

sss=s1(~isnan(s1));
sta(i,1:4)=[median(sss) lc uc z];
end
sta(:,1:3)=round(sta(:,1:3)*10,1);
sta(:,4)=round(sta(:,4),3);



