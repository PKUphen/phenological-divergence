clear all
clc
Ta=ncread('I:\ERA monthly\ERA5 monthly averaged data on single levels from 1979 to present .nc','t2m')-273.15;%1979年以来的气温
yrtmp=permute(Ta(:,:,1,:),[2,1,4,3]);
Ta=ncread('F:\ERA monthly\ERA5 monthly averaged data on single levels from 1979 to present .nc','tp');%降水
yrtp=permute(Ta(:,:,1,:),[2,1,4,3]);
clear Ta

phenstan=importdata('I:\workspace\animalphennumdata.mat');%1-7列分别为数据集序号，站号，物种号，物候期号，年份，物候期，物候期大类编号
[phcode,phlist]=xlsread('I:\workspace\动物物候期分类.xlsx','Sheet1','P2:R11'); %phcode=1:1:10; phlist={'Spring departure','First occurrence','Peak occurrence','Nest','Foraging','Breeding','Young individual','Molt','Autumn arrive','Last occurrence'}';
stlist=importdata('I:\workspace\animalsite.mat');%物候站点序号，纬度，经度
stlist(stlist(:,3)<0,3)=360+stlist(stlist(:,3)<0,3);

%植物
% phenstan=importdata('I:\workspace\plantphennumdata.mat');
% [phcode,phlist]=xlsread('I:\workspace\植物物种_物候期分类.xlsx','植物','Q2:T27');
% [phc,phli]=xlsread('I:\workspace\植物物种_物候期分类.xlsx','植物','I2:O254');
% stlist=importdata('I:\workspace\plantsite.mat');
% stlist(stlist(:,3)<0,3)=360+stlist(stlist(:,3)<0,3);
% for i=1:length(phc)
%     disp(['Processing the',num2str(i),'phenophase;']);
%     phc(i,4)=phcode(phcode(:,1)==phc(i,5),3);
%     phenstan(phenstan(:,4)==phc(i,1),7)= phc(i,4);
% end

sta=[];
nrc=length(phlist);
tims=unique(phenstan(:,[2 3 4 7]),'rows');
for ph1=1:nrc-1
    for ph2=ph1+1:nrc
        tim_tran_e=tims(tims(:,4)==ph1,:);%较早物候期的站-种-物候期序列
        tim_tran_l=tims(tims(:,4)==ph2,:);%较晚物候期的站-种-物候期序列
        phen_tran_e=phenstan(phenstan(:,7)==ph1,2:7);%较早物候期的站-种-物候期-年份数据
        phen_tran_l=phenstan(phenstan(:,7)==ph2,2:7);%较晚物候期的站-种-物候期-年份数据
        stsp_tim=unique(tim_tran_e(:,[1 2 4]),'rows');%较早物候期的站-种_物候分类序列
        for stp=1:size(stsp_tim,1)
            disp(['Processing the',num2str(ph1),'phenclass',num2str(stp),'tim series/',num2str(size(stsp_tim,1)),';']);
            pp_e=unique(tim_tran_e(tim_tran_e(:,1)==stsp_tim(stp,1)&tim_tran_e(:,2)==stsp_tim(stp,2),3));%较早物候的站点-物种对应的物候期数量
            pp_l=unique(tim_tran_l(tim_tran_l(:,1)==stsp_tim(stp,1)&tim_tran_l(:,2)==stsp_tim(stp,2),3));%较晚物候的站点-物种对应的物候期数量

            yr=[1980:1:2021]';
            ph_olp=nan(length(yr),2);%设定物候序列为空值变量

            if isempty(pp_l)==0
                rtran=nan(length(pp_e),length(pp_l));%物候期之间简单相关系数
                p=nan(length(pp_e),length(pp_l));%物候期之间简单相关显著性p-value


                prtran=nan(length(pp_e),length(pp_l));%物候期与较早物候期的偏相关系数，去除气温和降水影响
                pp=nan(length(pp_e),length(pp_l));%物候期与较早物候期的偏相关显著性p-value，去除气温和降水影响
                tprtran=nan(length(pp_e),length(pp_l));%物候期与气温的偏相关系数，去除较早物候期和降水影响
                tpp=nan(length(pp_e),length(pp_l));%物候期与气温的偏相关显著性p-value，去除较早物候期和降水影响
                pprtran=nan(length(pp_e),length(pp_l));%物候期与降水的偏相关系数，去除较早物候期和气温影响
                ppp=nan(length(pp_e),length(pp_l));%物候期与降水的偏相关显著性p-value，去除较早物候期和气温影响

                EPhtran=nan(length(pp_e),length(pp_l));%物候期与较早物候期的偏相关系数，去除环境的影响（温度和降水取两者最大偏相关系数）
                ENtran=nan(length(pp_e),length(pp_l));%物候期与环境因子的偏相关系数，去除较早物候期的影响


                slotran=nan(length(pp_e),length(pp_l));%物候期间隔的变化趋势
                slop=nan(length(pp_e),length(pp_l));%物候期间隔趋势的显著性

                for i=1:length(pp_e)
                    for j=1:length(pp_l)
                        phe=phen_tran_e(phen_tran_e(:,1)==stsp_tim(stp,1)&phen_tran_e(:,2)==stsp_tim(stp,2)&phen_tran_e(:,3)==pp_e(i),4:5);
                        ph_olp(phe(:,1)-1979,2)=phe(:,2);%较早物候期数据
                        phl=phen_tran_l(phen_tran_l(:,1)==stsp_tim(stp,1)&phen_tran_l(:,2)==stsp_tim(stp,2)&phen_tran_l(:,3)==pp_l(j),4:5);
                        ph_olp(phl(:,1)-1979,1)=phl(:,2);%较晚物候期数据
                        rl=stlist(stlist(:,1)==stsp_tim(stp,1),2);%站点纬度
                        cl=stlist(stlist(:,1)==stsp_tim(stp,1),3);%站点经度
                        rl=round((90.125-rl)/0.25);%站点所在像元行号
                        cl=round((0.125+cl)/0.25);%站点所在像元列号

                        nyr=mean(ph_olp,2);
                        ph_olp(isnan(nyr),1)=nan;%数据不重叠的年份物候日期设置为nan
                        ph_olp(isnan(nyr),2)=nan;%数据不重叠的年份物候日期设置为nan

                        phml=nanmean(ph_olp(:,1));%较晚物候期均值
                        phme=nanmean(ph_olp(:,2));%较早物候期均值

                        if day(phml)>15 %较晚日期平均值晚于15日的采用当月，早于15日的采用前一个月
                            lmon=month(phml)+12;
                        else
                            lmon=month(phml)-1+12;
                        end
                        if phml<0 %物候发生日期为负数，则为前一年的月份
                            lmon=lmon-12;
                        else if phml>=367%物候发生日期大于367天，则为后一年的月份
                                lmon=lmon+12;
                        end
                        end

                        if day(phme)>15 %较早日期平均值晚于15日的采用次月，早于15日的采用当月
                            emon=month(phme)+1+12;
                        else
                            emon=month(phme)+12;
                        end

                        if phme<0%物候发生日期为负数，则为前一年的月份
                            emon=emon-12;
                        else if phme>=367%物候发生日期大于367天，则为后一年的月份
                                emon=emon+12;
                        end
                        end

                        if lmon-emon>=0&&lmon-emon<=12
                            mon=emon:1:lmon;
                            for mm=1:length(mon)
                                tran=yrtmp(rl,cl,mon(mm):12:mon(mm)+12*41);
                                tyr(:,mm)=tran(:);
                                tran=yrtp(rl,cl,mon(mm):12:mon(mm)+12*41);
                                pyr(:,mm)=tran(:);
                            end
                            tyr=nanmean(tyr,2);
                            pyr=nansum(pyr,2);
                        end


                        if length(nyr(~isnan(nyr)))>=5&&phml>phme
                            [rtran(i,j),p(i,j)]=corr(ph_olp(~isnan(nyr),1),ph_olp(~isnan(nyr),2));
                            [prtran(i,j),pp(i,j)]=partialcorr(ph_olp(~isnan(nyr),1),ph_olp(~isnan(nyr),2),[tyr(~isnan(nyr)) pyr(~isnan(nyr))]);
                            [tprtran(i,j),tpp(i,j)]=partialcorr(ph_olp(~isnan(nyr),1),tyr(~isnan(nyr)),[ph_olp(~isnan(nyr),2) pyr(~isnan(nyr))]);
                            [pprtran(i,j),ppp(i,j)]=partialcorr(ph_olp(~isnan(nyr),1),pyr(~isnan(nyr)),[ph_olp(~isnan(nyr),2) tyr(~isnan(nyr))]);

                            [Ct,~]=partialcorr(ph_olp(~isnan(nyr),1),tyr(~isnan(nyr)),ph_olp(~isnan(nyr),2));
                            [Cp,~]=partialcorr(ph_olp(~isnan(nyr),1),pyr(~isnan(nyr)),ph_olp(~isnan(nyr),2));
                            [ENtran(i,j),loc]=max([abs(Ct) abs(Cp)]);
                            if loc==1
                                EPhtran(i,j)=abs(partialcorr(ph_olp(~isnan(nyr),1),ph_olp(~isnan(nyr),2),tyr(~isnan(nyr))));
                            else
                                EPhtran(i,j)=abs(partialcorr(ph_olp(~isnan(nyr),1),ph_olp(~isnan(nyr),2),pyr(~isnan(nyr))));
                            end

                            [b,~,~,~,stats]=regress(abs(ph_olp(:,1)-ph_olp(:,2)),[yr ones(42,1)]);
                            slotran(i,j)=b(1);
                            slop(i,j)=stats(3);
                        end
                    end
                end

                if ~isnan(nanmean(p(:)))
                    de=abs(prtran)+abs(tprtran)+abs(pprtran);
                    Prtran=abs(prtran)./de;
                    Etprtran=abs(tprtran)+abs(pprtran)./de;
                    sta=[sta; [stsp_tim(stp,:) ph2 nanmedian(rtran(:)) nanmean(rtran(:)) length(p(p<0.05))./length(p(~isnan(p))) nanmedian(prtran(:)) nanmean(prtran(:)) length(pp(pp<0.05))./length(pp(~isnan(pp))) nanmedian(tprtran(:)) nanmean(tprtran(:)) length(tpp(tpp<0.05))./length(tpp(~isnan(tpp))) nanmedian(pprtran(:)) nanmean(pprtran(:)) length(ppp(ppp<0.05))./length(ppp(~isnan(ppp))) nanmedian(Prtran(:)) nanmean(Prtran(:)) nanmedian(Etprtran(:)) nanmean(Etprtran(:)) nanmedian(EPhtran(:)) nanmean(EPhtran(:)) nanmedian(ENtran(:)) nanmean(ENtran(:)) nanmedian(slotran(:)) nanmean(slotran(:)) length(slop(slop<0.05))./length(slop(~isnan(slop)))]];
                end

            end
        end
    end
end

sta(~isnan(sta(:,5))&isnan(sta(:,10)),10)=0;
sta(isnan(sta(:,5)),:)=[];
stace= [cell({'stID','spID','phclase','phclasl','rmedian','rmean','persignr', 'prmedian','prmean','ppersignr','tprmedian','tprmean','tppersignr','pprmedian','pprmean','pppersignr','Phrela3median','Phrela3mean','ENrela3median','ENrela3mean','Phmax2median','Phmax2mean','ENmax2median','ENmax2mean','slomedian','slomean','persignslo'}); mat2cell(sta,ones(size(sta,1),1),ones(1,27))];

save('I:\workspace\animal_phase_pcorrtp.mat','stace');
%save('I:\workspace\plant_phase_pcorrtp.mat','stace');