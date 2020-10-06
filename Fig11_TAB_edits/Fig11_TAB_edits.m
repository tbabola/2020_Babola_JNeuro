%% Figure 11: alpha9 GOF and KO mice
%% load MATLAB path
addpath(genpath('..\MATLAB Functions'));
%% example images for Panel A
% a9 KO
img = loadTif('.\Data\a9 KO example\e0293_sub2.tif',16);
%% normalize image
[~,Fo,img] = normalizeImg(img,5,1);
%%
blsubt = img - Fo;
%% 
meanImg = mean(blsubt(:,:,760:780),3);
figure; imagesc(meanImg); colormap gfb;
caxis([-50 800]);
writeTif(single(meanImg), '.\Data\a9 KO example\KO_example_760_780.tif',32)

%% WT image
img = loadTif('.\Data\WT example\Experiment-205_CJK_20190214 exp205 mouse 1 wt littermate-1_sub.tif',16);
%% normalize image
[~,Fo,img] = normalizeImg(img,5,1);
blsubt = img - Fo;
%%
meanImg = mean(blsubt(:,:,720:740),3);
figure; imagesc(meanImg); colormap gfb;
caxis([-50 800]);
writeTif(single(meanImg),'.\Data\WT example\WT_example_720_740.tif',32)
%% GOF image
img = loadTif('.\Data\a9 GOF Example\Experiment-1214_sub.tif',16);
%% normalize image
[~,Fo,img] = normalizeImg(img,5,1);
blsubt = img - Fo;
%%
meanImg = mean(blsubt(:,:,1232:1252),3);
figure; imagesc(meanImg); colormap gfb;
caxis([-50 800]);
writeTif(single(meanImg),'.\Data\a9 GOF Example\GOF_example_720_740.tif',32)

%% Actual Panels for A
%KO
img1 = loadTif('.\Data\a9 KO example\KO_example_760_780.tif',32);
figure; imagesc(imresize(imgaussfilt(img1),2)); truesize; colormap gfb; axis off;
caxis([-50 800]); figQuality(gcf,gca,[3,2]); truesize;
export_fig('.\EPS_Panels\A_KO_exemplar.eps')

%WT
img2 = loadTif('.\Data\WT example\WT_example_720_740.tif',32);
figure; imagesc(imresize(imgaussfilt(img2),2)); truesize; colormap gfb; axis off;
caxis([-50 800]); figQuality(gcf,gca,[3,2]); truesize;
export_fig('.\EPS_Panels\A_WT_exemplar.eps')

%GOF
img3 = loadTif('.\Data\a9 GOF example\GOF_example_720_740.tif',32);
figure; imagesc(imresize(imgaussfilt(img3),2));  colormap gfb; axis off;
caxis([-50 800]); figQuality(gcf,gca,[3,2]); truesize;
export_fig('.\EPS_Panels\A_GOF_exemplar.eps')

%% alpha 9 GOF
h = load('.\Data\a9 GOF Example\Experiment-1214\ICinfo16_dFoF.mat');
[LIC_bladj RIC_bladj] = getBladj(h.ICsignal)

color1 = [57 107 43]/255;
color2 = color1 + [0.4 0.4 0.4];
figure; plot(LIC_bladj, 'color',color1,'LineWidth',1.25);
hold on; plot(RIC_bladj,'color',color2,'LineWidth',1.25);
xlim([0 2400]);
ylim([-.05 .3])
figQuality(gcf,gca,[6 0.75]);
export_fig('.\EPS_Panels\a9GOF_example.eps')


lollipopPlot(h.pkData)
%% WT
h = load('.\Data\WT example\Experiment-205\ICinfo16_dFoF.mat');
[LIC_bladj RIC_bladj] = getBladj(h.ICsignal)

color1 = [0 0 0];
color2 = [0.6 0.6 0.6];
figure; plot(LIC_bladj, 'color',color1,'LineWidth',1.25);
hold on; plot(RIC_bladj,'color',color2,'LineWidth',1.25);
xlim([0 2400]);
ylim([-.05 .3])
figQuality(gcf,gca,[6 0.75]);
export_fig('.\EPS_Panels\WT_example.eps')

lollipopPlot(h.pkData)
%% a9 KO
h = load('.\Data\a9 KO example\Experiment-0293\ICinfo16_dFoF.mat');
[LIC_bladj RIC_bladj] = getBladj(h.ICsignal)

color1 = [1 0 0];
color2 = [1 0.6 0.6];
figure; plot(LIC_bladj, 'color',color1,'LineWidth',1.25);
hold on; plot(RIC_bladj,'color',color2,'LineWidth',1.25);
xlim([2400 4800]);
ylim([-.05 .3])
figQuality(gcf,gca,[6 0.75]);
export_fig('.\EPS_Panels\a9KO.eps')

lollipopPlot(h.pkData)
%%

function [LIC_bladj RIC_bladj] = getBladj(ICsignal)
    LIC = ICsignal(:,1);
    RIC = ICsignal(:,2);
    time = [1:1:size(LIC,1)]';
    LIC_bladj = msbackadj(time, LIC);
    RIC_bladj = msbackadj(time,RIC);
end