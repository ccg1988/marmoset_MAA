% @ 2022-12-01
clear;clc;close all
load('raw_data.mat')
load('psignifit_options.mat')
% 'norm'        a cummulative Gaussian distribution
% 'logistic'    a logistic function
% 'gumbel'      a cummulative gumbel distribution
% 'rgumbel'     a reversed gumbel distribution*****USED*****
% 'tdist'       a t-distribution with df=1 as a heavytail distribution
%******* for positive stimulus levels which make sense on a log-scale:
% 'logn'        a cumulative lognormal distribution
% 'Weibull'     a Weibull function
N=5;% five animals
M3T=[255,131,104,150]; M94W=[35,189,255,150]; M71V=[178,138,0,150];
M76X=[142,17,137,150]; M63W=[0,153,54,150]; C=[M3T;M94W;M71V;M76X;M63W]/255;

x_deg=([7.5;15;22.5;30;37.5;45;90]);
fitX=nan(N, 1200);
fitY=nan(N, 2);
fitValues=nan(N, 1200);
fitPoints=nan(N, 7);
fitParas=nan(N,4);
cond=1; % azi_Gau_front/RSS_front/Gau_rear***ele_Gau_2to32/4to26/4to12kHz
result_save=cell(N,1);
for n = 1 : N
nCorrect=T{n+(cond-1)*5, 5:11}';   
total=T{n+(cond-1)*5, 13:19}';
fitPoints(n,:)=nCorrect./total;
data=[x_deg,nCorrect,total]; % columns==3(x | nCorrect | total)
result = psignifit(data, options);
[fitY(n,:),fitX(n,:),fitValues(n,:)]=plotPsych(result,plotOptions);
fitParas(n,:)=result.Fit(1:4);%threshold+width+lapse+guess
result_save{n}=result;
end
%%
x_deg_2d=repmat(x_deg,1,5)';
x_deg_2d([1 2],5)=38.5;
x_deg_2d(3,6)=43;
x_deg_2d([3 5],7)=88;
%%
figure; CR=nan(N,1);
for n = 1 : N
CR(n)=T{n+(cond-1)*N, 3};  
plot(CR(n), 0.5,'LineStyle','none','Marker','square','MarkerSize',10,...
    'MarkerFaceColor','k','MarkerEdgeColor',C(n,1:3));hold on       
HR=fitPoints(n,:);
FAR=T{n+(cond-1)*5, 21}/T{n+(cond-1)*5, 22};
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
    'Color',C(n,1:3));hold on
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
%% marginal posterior density for a single parameter
% A sharp curve with narrow width means better fitting
% The gray patch shadow shows the 95% chosen confidence interval 
% The black line shows the point estimate of paramater(used in previous figure)

% The gray dashed line shows the prior. Psignifit chooses a prior which assumes that 
% you sampled the whole psychometric function.Specifically, it assumes that the threshold 
% is within the range of the data and with decreasing probability up to half the range 
% above or below the measured data. 
% Idealy, we should record a trial well above and a trial well below threshold. 
% and posterior concentrates on an area (aks gray patch) for which the prior was constant.
% For good fitting, there is no area where the prior dominates the posterior. 
% For bad fitting, the prior goes down where there is still posterior probability. 
% Thus, the prior has an influence on the outcome.
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
% 1 = threshold, 2 = width, 3 = lambda, 4 = gamma, 5 = eta
plotMarginal(result_save{n}, 1, plotOptions); 
xlim([-25 50])
end
end  
%% 2 dimensional posterior marginals
% plotOptions.label1    = 'Threshold (\circ)';   % label for the first parameter
% plotOptions.label2    = 'Lapse rate';   % label for the second parameter
% pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
% figure('position',[X_size*0.02 Y_size*0.06 X_size*0.95 Y_size*0.55]);
% tiledlayout(1,5);
% for n = 1 : N
% nexttile    
% plot2D(result_save{n},1,3,plotOptions);title(subject{n})
% % xlim([0 0.3]);ylim([0 90])
% end    