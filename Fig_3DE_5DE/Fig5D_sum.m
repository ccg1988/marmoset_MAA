%modified the code from
%PAL_PFHB_SingleSubjectsMultipleConditions_Demo  
%This demo code fits multiple psychometric functions to data derived 
%from a single subject testing in multiple conditions using a Bayesian 
%criterion. The model that is fitted is shown here:
%www.palamedestoolbox.org/hierarchicalbayesian.html)
%JAGS (http://mcmc-jags.sourceforge.net/) must first be
%installed before this will work. JAGS will perform the MCMC sampling of the posterior. 

clc; clear; 
close all;
load('raw_data.mat')
engine = 'jags';
subject1=[1,2,3];%subject 1-3
subject2=[1,2,3,4,5];%subject 1-5
Ncond=2;
M = PAL_Contrasts(Ncond,'polynomial');
T{17,5}=T{16,5}/2+T{18,5}/2; T{17,13}=T{16,13}/2+T{18,13}/2; %M94W missed 7.5 degree
T{25,5}=mean(T{21:24,5}); T{25,13}=mean(T{21:24,13}); %M63W missed 7.5 degree
% Gaussion 2-32kHz vs 4-26kHz
data.y = round([sum(T{subject1+5*3,5:12}), sum(T{subject2+5*4,5:12})]);
data.n = round([sum(T{subject1+5*3,13:20}), sum(T{subject2+5*4,13:20})]);
data.x = ([7.5 15 22.5 30 37.5 45 90 180 7.5 15 22.5 30 37.5 45 90 180]);
data.c = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2] ;        
pfhb = PAL_PFHB_fitModel(data,'PF','gumbel',...
    'lapse','unconstrained','guess','fixed',0,'engine',engine,...
    'nsamples',10000);
PAL_PFHB_inspectFit(pfhb,'condition',1);
PAL_PFHB_inspectFit(pfhb,'condition',2);
PAL_PFHB_inspectParam(pfhb,'a','condition',1,'a','condition',2);
%%
%The following shows diagnostics and derives summary statistics for the 
%difference between the location parameters (aka 'threshold') in conditions 
%1 and 2. Also shows scatter plot between these two parameters.

% PAL_PFHB_inspectParam(pfhb,'a','condition',1,'l','condition',1);

%The following should lead to (essentially) identical fits, but samples
%from posterior of reparameterized location parameters.
pfhb_reparam=PAL_PFHB_fitModel(data,'PF','gumbel',...
    'lapse','unconstrained','guess','fixed',0,'engine',engine,...
    'a',M,'nsamples',10000);
% Dot corresponds to central tendency ({mode}, median, or mean) in 
% marginal posterior distribution, line covers 68% high-density interval 
% (hdi, or credible interval), curves show posterior across its 95% hdi.
PAL_PFHB_drawViolins(pfhb_reparam,'centralTendency','mean')