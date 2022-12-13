% @ 2022-12-01
clear;clc;close all
load('raw_data.mat')
load('psignifit_options.mat')
N=1;% ONE animals
M3T=[255,131,104,150]; M94W=[35,189,255,150]; M71V=[178,138,0,150];
M76X=[142,17,137,150]; M63W=[0,153,54,150]; C=[M3T;M94W;M71V;M76X;M63W]/255;

x_deg=[7.5;15;22.5;45;90];
n=1:4; %only first 3 animals
cond=3; % azi_Gau_front/RSS_front/Gau_rear***ele_Gau_2to32/4to26/4to12kHz
nCorrect=sum(T{n+(cond-1)*5, [5:7, 10,11]})';   
total=sum(T{n+(cond-1)*5, [13:15, 18,19]})';
fitPoints=nCorrect./total;
data=[x_deg,nCorrect,total]; % columns==3(x | nCorrect | total)
result = psignifit(data, options);
[fitY,fitX,fitValues]=plotPsych(result,plotOptions);
fitParas=result.Fit(1:4);%threshold+width+lapse+guess
result_save=result;
%%
figure;
CR=mean(T{n+(cond-1)*N, 3});  
plot(CR, 0.5,'LineStyle','none','Marker','square','MarkerSize',10,...
    'MarkerFaceColor','k','MarkerEdgeColor',[0 0 0]+0.5);hold on       
HR=fitPoints;
FAR=sum(T{n+(cond-1)*5, 21})/sum(T{n+(cond-1)*5, 22});
HRc=(HR-FAR)/(1-FAR);
plot(x_deg,HRc,"LineStyle","-.",'LineWidth',3,'Marker','.','MarkerSize',25,...
    'Color',rgb('Grey'));hold on   
plot([0 100],[0.5 0.5],"LineStyle","-.",'LineWidth',1,'Color','k')
lgd= legend({'',['All=',num2str(round(CR,2)),'\circ']});  
lgd.Box='off'; lgd.FontSize=12; lgd.Location='southeast';
xlim([0 90]); xlabel('Angle (\circ)'); ylim([0 1]); ylabel('Corrected hit rates');
xticks(0:x_deg(1):x_deg(end))
%%
figure
plot(fitParas(1), fitY(2),'LineStyle','none','Marker','square','MarkerSize',10,...
    'MarkerFaceColor','k','MarkerEdgeColor','k');hold on
plot(x_deg,fitPoints,"LineStyle","none",'Marker','.','MarkerSize',25,...
    'Color', 'k');hold on
plot(fitX,fitValues,'LineStyle','-.','LineWidth',3,'Color',rgb('Grey'));hold on 
lgd= legend({'',['Averaged',num2str(round(fitParas(1),2)),...
    '\circ lapse=',num2str(round(fitParas(3),2))]}); 
lgd.Box='off'; lgd.FontSize=12; lgd.Location='southeast';
xlim([0 90]); xlabel('Angle (\circ)'); ylabel('Percent Correct');
xticks(0:x_deg(1):x_deg(end))
%% marginal posterior density for a single parameter
W = 28/5 ; %centimeter
H = 5 ; %centimeter
F_posi = [10, 10, W, H] ; %X-Units to right of monitor, Y-Units above bottom of monitor
fig = figure; 
fig.Units = "centimeters";
fig.Color = "White"; %evenif using 'None', there is a still white-filled square
fig.InnerPosition = F_posi ; %X-Units to right, Y-Units above bottom; figure wide & tall
fig.PaperSize = fig.Position(3:4) ;
fig.PaperUnits = "centimeters" ;

subject='All';
plotOptions.priorColor     = [.2,.2,.2];          % color of the prior distibution
plotOptions.xLabel         = 'MAA (\circ)';    % X-Axis label
plotOptions.yLabel         = [];
plotOptions.lineWidth      = 1;
tiledlayout(1,1,TileSpacing="tight",Padding="tight");
nexttile;title(subject)  
plotOptions.lineColor = [0 0 0]; % color of the density
% 1 = threshold, 2 = width, 3 = lambda, 4 = gamma, 5 = eta
plotMarginal(result_save, 1, plotOptions); 
xlim([0 25])
 