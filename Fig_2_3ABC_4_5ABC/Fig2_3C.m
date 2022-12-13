% @ 2022-12-02
clear;clc;close all
load('raw_data.mat')
load('psignifit_options.mat')
N=5;% five animals
M3T=[255,131,104,150]; M94W=[35,189,255,150]; M71V=[178,138,0,150];
M76X=[142,17,137,150]; M63W=[0,153,54,150]; C=[M3T;M94W;M71V;M76X;M63W]/255;

x_deg=[7.5;15;22.5;45;90];
fitX=nan(N, 1200);
fitY=nan(N, 2);
fitValues=nan(N, 1200);
fitPoints=nan(N, length(x_deg));
fitParas=nan(N,4);
cond=3; % azi_Gau_front/RSS_front/Gau_rear***ele_Gau_2to32/4to26/4to12kHz
result_save=cell(N,1);
options.confP=0.68;
for n = 1 : N
nCorrect=T{n+(cond-1)*N, 5:11}';nCorrect([4,5])=[]; 
total=T{n+(cond-1)*N, 13:19}'; total([4,5])=[]; 
if sum(nCorrect)>0 % values==-1 for empty subject/animal
%miss_deg=find(nCorrect==-1); 
% if ~isempty(miss_deg)
% nCorrect(miss_deg)=[];total(miss_deg)=[];
% end    
fitPoints(n,:)=nCorrect./total;
data=[x_deg,nCorrect,total]; % columns==3(x | nCorrect | total)
result = psignifit(data,options);
[fitY(n,:),fitX(n,:),fitValues(n,:)]=plotPsych(result,plotOptions);
fitParas(n,:)=result.Fit(1:4);%threshold+width+lapse+guess
result_save{n}=result;
end
end
x_deg_2d=repmat(x_deg,1,5)';
%% dependent on specific figures
x_deg_2d(4,2)=16;
x_deg_2d(4,4)=47;
%% 50% method
figure
% CR1=T{1+(cond-1)*5, 3};CR2=T{2+(cond-1)*5, 3};CR3=T{3+(cond-1)*5, 3};
% CR4=T{4+(cond-1)*5, 3};CR5=T{5+(cond-1)*5, 3};
CR=nan(N,1);
for n = 1 : N
CR(n)=T{n+(cond-1)*N, 3};  
plot(CR(n), 0.5,'LineStyle','none','Marker','square','MarkerSize',10,...
    'MarkerFaceColor','k','MarkerEdgeColor',C(n,1:3));hold on    
HR=fitPoints(n,:);
FAR=T{n+(cond-1)*N, 21}/T{n+(cond-1)*N, 22};
HRc=(HR-FAR)/(1-FAR);
plot(x_deg_2d(n,:),HRc,"LineStyle","-.",'LineWidth',3,'Marker','.','MarkerSize',25,...
    'Color',C(n,:));hold on
end   
plot([0 100],[0.5 0.5],"LineStyle","-.",'LineWidth',1,'Color','k')
lgd= legend({'',['M3T MAA=',num2str(round(CR(1),2)),'\circ'],...
    '',['M94W MAA=',num2str(round(CR(2),2)),'\circ'],...
    '',['M71V MAA=',num2str(round(CR(3),2)),'\circ'],...
    '',['M76X MAA=',num2str(round(CR(4),2)),'\circ'],...
    '',['M63W MAA=',num2str(round(CR(5),2)),'\circ']}); 
lgd.Box='off'; lgd.FontSize=12; lgd.Location='southeast';
xlim([0 90]); xlabel('Angle (\circ)'); ylim([0 1]); ylabel('Corrected hit rates');
xticks(0:x_deg(1):x_deg(end))
%%
figure
for n = 1 : N
% %plot([fitParas(n,1), fitParas(n,1)],fitY(n,:),'LineWidth',1,'Color','k');hold on
plot(fitParas(n,1), fitY(n,2),'LineStyle','none','Marker','square','MarkerSize',10,...
    'MarkerFaceColor','k','MarkerEdgeColor',C(n,1:3));hold on
plot(x_deg_2d(n,:),fitPoints(n,:),"LineStyle","none",'Marker','.','MarkerSize',25,...
    'Color',C(n,:));hold on
plot(fitX(n,:),fitValues(n,:),'LineStyle','-.','LineWidth',3,'Color',C(n,:));hold on
end  
lgd= legend({'',['M3T MAA=',num2str(round(fitParas(1,1),2)),...
    '\circ lapse=',num2str(round(fitParas(1,3),2))],...
    '','',['M94W MAA=',num2str(round(fitParas(2,1),2)),...
    '\circ lapse=',num2str(round(fitParas(2,3),2))],'','',...
    ['M71V MAA=',num2str(round(fitParas(3,1),2)),...
    '\circ lapse=',num2str(round(fitParas(3,3),2))],'','',...
    ['M76X MAA=',num2str(round(fitParas(4,1),2)),...
    '\circ lapse=',num2str(round(fitParas(4,3),2))],'','',...
    ['M63W MAA=',num2str(round(fitParas(5,1),2)),...
    '\circ lapse=',num2str(round(fitParas(5,3),2))],''}); 
lgd.Box='off'; lgd.FontSize=12; lgd.Location='southeast';
xlim([0 90]); xlabel('Angle (\circ)'); ylabel('Percent Correct');
xticks(0:x_deg(1):x_deg(end))
%%
W = 17.4 ; %centimeter
H = 5 ; %centimeter
F_posi = [10, 10, W, H] ; %X-Units to right of monitor, Y-Units above bottom of monitor
fig = figure; 
fig.Units = "centimeters";
fig.Color = "White"; %evenif using 'None', there is a still white-filled square
fig.InnerPosition = F_posi ; %X-Units to right, Y-Units above bottom; figure wide & tall
fig.PaperSize = fig.Position(3:4) ;
fig.PaperUnits = "centimeters" ;

subject={'M3T';'M94W';'M71V';'M76X';'M63W'};
plotOptions.priorColor     = [.2,.2,.2];          % color of the prior distibution
plotOptions.xLabel         = 'MAA (\circ)';    % X-Axis label
plotOptions.yLabel         = [];
plotOptions.lineWidth      = 1;
tiledlayout(1,5,TileSpacing="tight",Padding="tight");
for n = 1 : N
nexttile;title(subject{n})  
plotOptions.lineColor = C(n,1:3); % color of the density
if ~isempty(result_save{n})
plotMarginal(result_save{n}, 1, plotOptions); 
xlim([-25 50])
end
end  