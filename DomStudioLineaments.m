function DomStudioLineaments(~)
% @ 2023 by Andrea Bistacchi & Stefano Casiraghi, distributed under the GNU AGPL v3.0 license.
%
% Function used for the analysis of lineaments imported as SHP files
%
% Last update 2023/04/19

% _____________
% 0- initialize
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')
rad = pi/180;
deg = 180/pi;

% ______________
% 1 - input data

% load SHP with lineament data as a Geographic Data Structure
% (https://www.mathworks.com/help/map/geographic-data-structures.html)
[file, path] = uigetfile('*.shp');
filename = [path file];
[path,file,~] = fileparts(filename);
inDataStuct = shaperead(filename);
disp(['File ' file '.shp read.'])
disp(' ')

% input dip/direction of mean fracture plane
disp('Input mean attitude of fracture set.')
meanDip = -1;
while meanDip<=0
    meanDip = input('Mean DIP of fracture set [90]: ');
    if isempty(meanDip), meanDip = 90; end
end

meanDir = -1;
while meanDir<=0
    meanDir = input('Mean DIRECTION of fracture set [90]: ');
    if isempty(meanDir), meanDir = 90; end
end
disp(' ');

% outcrop size (just used for visualization)
outcropHeight = -1;
while outcropHeight<=0
    outcropHeight = input('Outcrop height along dip [10]: ');
    if isempty(outcropHeight), outcropHeight = 10; end
end

outcropLength = -1;
while outcropLength<=0
    outcropLength = input('Outcrop length along strike [10]: ');
    if isempty(outcropLength), outcropLength = 19; end
end
disp(' ');

% ________________________________________________________________________
% 1 - calculate mean plunge/trend and xyz components of normal unit vector
meanP = 90 - meanDip;
meanT = meanDir+180 -360*(meanDir+180>360);
meanNormal = [sin(meanT*rad)*cos(meanP*rad)...
    cos(meanT*rad)*cos(meanP*rad)...
    sin(meanP*rad)];

% ___________________________________
% 2 - extract and clean polyline data

% Ndata = number of data points
Nlines = size(inDataStuct);

% create cell array for polylines, removing NaNs and other invalid values
lines = {};
for i = 1:Nlines
    lines{i} = [inDataStuct(i).X(isfinite(inDataStuct(i).X))' inDataStuct(i).Y(isfinite(inDataStuct(i).X))' zeros(size(find(isfinite(inDataStuct(i).X))))'];
end

% _______________________________________________________________
% 3 - calculate centroids, outcrop center, and shift to re-center

% calculate centroids
centroids = [];
for i = 1:Nlines
    centroids(i,:) = mean(lines{i});
end

% center of dataset from mean of facet centers - all following processing
% will assume this point as origin of the reference frame
outcropCenter = mean(centroids);

% translate lines in order to have outcropCenter = [0 0 0]
for i = 1:Nlines
    lines{i} = lines{i} - outcropCenter;
end

% re-calculate centroids
centroids = [];
for i = 1:Nlines
    centroids(i,:) = mean(lines{i});
end

% re-calculate outcrop center to check
outcropCenter = mean(centroids);

% _____________________________________
% 4 - interpolare average outcrop plane

% interpolate average outcrop plane with pincipal components
[coeff,~,~] = pca(centroids);
outcropNormal = coeff(:,3)';

% plunge/trend of outcrop plane
outcropP = asin(-outcropNormal(3))*deg;
outcropT = atan2(outcropNormal(1),outcropNormal(2))*deg;
outcropT = outcropT+(outcropT<0)*360;
outcropP = outcropP*(outcropP>=0)-outcropP*(outcropP<0);
outcropT = outcropT-(outcropP<0)*180;
outcropT = outcropT+(outcropT<0)*360;
outcropT = outcropT-(outcropT>360)*360;
disp(['outcropP = ' num2str(outcropP)]);
disp(['outcropT = ' num2str(outcropT)]);

% dip/dip azimuth of outcrop plane
outcropDip = 90-outcropP;
outcropDipAzimuth = outcropT+180;
outcropDipAzimuth = outcropDipAzimuth+(outcropDipAzimuth<0)*360;
outcropDipAzimuth = outcropDipAzimuth-(outcropDipAzimuth>360)*360;
disp(['outcropDip = ' num2str(outcropDip)]);
disp(['outcropDipAzimuth = ' num2str(outcropDipAzimuth)]);

% outcrop strike
outcropStrike = (outcropDipAzimuth<90).*(outcropDipAzimuth+270) + (outcropDipAzimuth>=90).*(outcropDipAzimuth-90);
disp(['outcropStrike = ' num2str(outcropStrike)]);

% ______________________________
% 5 - calculate projection plane

% intersection (unit vector) of mean fracture plane and outcrop plane
facetOutcropIntersection = cross(meanNormal,outcropNormal);
facetOutcropIntersectionMod = sqrt(facetOutcropIntersection(1)^2+facetOutcropIntersection(2)^2+facetOutcropIntersection(3)^2);
facetOutcropIntersection = facetOutcropIntersection/facetOutcropIntersectionMod;

% plunge/trend of intersection
intersectionPlaneP = asin(-facetOutcropIntersection(3))*deg;
intersectionPlaneT = atan2(facetOutcropIntersection(1),facetOutcropIntersection(2))*deg;
intersectionPlaneT = intersectionPlaneT+(intersectionPlaneT<0)*360;
intersectionPlaneP = intersectionPlaneP*(intersectionPlaneP>=0)-intersectionPlaneP*(intersectionPlaneP<0);
intersectionPlaneT = intersectionPlaneT-(intersectionPlaneP<0)*180;
intersectionPlaneT = intersectionPlaneT+(intersectionPlaneT<0)*360;
intersectionPlaneT = intersectionPlaneT-(intersectionPlaneT>360)*360;
disp(['intersectionPlaneP = ' num2str(intersectionPlaneP)]);
disp(['intersectionPlaneT = ' num2str(intersectionPlaneT)]);

% normal to mean normal and intersection
projectionPlaneNormal = cross(facetOutcropIntersection,meanNormal);
projectionPlaneNormalMod = sqrt(projectionPlaneNormal(1)^2+projectionPlaneNormal(2)^2+projectionPlaneNormal(3)^2);
projectionPlaneNormal = projectionPlaneNormal/projectionPlaneNormalMod;

% plunge/trend of projection plane
projectionPlaneP = asin(-projectionPlaneNormal(3))*deg;
projectionPlaneT = atan2(projectionPlaneNormal(1),projectionPlaneNormal(2))*deg;
projectionPlaneT = projectionPlaneT+(projectionPlaneT<0)*360;
projectionPlaneP = projectionPlaneP*(projectionPlaneP>=0)-projectionPlaneP*(projectionPlaneP<0);
projectionPlaneT = projectionPlaneT-(projectionPlaneP<0)*180;
projectionPlaneT = projectionPlaneT+(projectionPlaneT<0)*360;
projectionPlaneT = projectionPlaneT-(projectionPlaneT>360)*360;
disp(['projectionPlaneP = ' num2str(projectionPlaneP)]);
disp(['projectionPlaneT = ' num2str(projectionPlaneT)]);

% dip/dip azimuth of projection plane
projectionPlaneDip = 90-projectionPlaneP;
projectionPlaneDipAzimuth = projectionPlaneT+180;
projectionPlaneDipAzimuth = projectionPlaneDipAzimuth+(projectionPlaneDipAzimuth<0)*360;
projectionPlaneDipAzimuth = projectionPlaneDipAzimuth-(projectionPlaneDipAzimuth>360)*360;
disp(['projectionPlaneDip = ' num2str(projectionPlaneDip)]);
disp(['projectionPlaneDipAzimuth = ' num2str(projectionPlaneDipAzimuth)]);

% projection plane strike
projectionPlaneStrike = (projectionPlaneDipAzimuth<90).*(projectionPlaneDipAzimuth+270) + (projectionPlaneDipAzimuth>=90).*(projectionPlaneDipAzimuth-90);
disp(['projectionPlaneStrike = ' num2str(projectionPlaneStrike)]);

% _______________________________________________
% 6 - rectangular frame for average outcrop plane

% downdip vector with length = 1/2 Height
halfOutcropPlaneDip(1) = cos(outcropDip*rad)*sin(outcropDipAzimuth*rad)*outcropHeight/2;
halfOutcropPlaneDip(2) = cos(outcropDip*rad)*cos(outcropDipAzimuth*rad)*outcropHeight/2;
halfOutcropPlaneDip(3) = sin(outcropDip*rad)*outcropHeight/2;  % in the facets version there was a "-" in front of this

% alongstrike vector with length = 1/2 Length
halfOutcropPlaneStrike(1) = sin(outcropStrike*rad)*outcropLength/2;
halfOutcropPlaneStrike(2) = cos(outcropStrike*rad)*outcropLength/2;
halfOutcropPlaneStrike(3) = 0;

% four corners of outcrop best fit plane, which is centerd in [0 0 0] = translated outcropCenter
outcropCorner1 = - halfOutcropPlaneStrike - halfOutcropPlaneDip;
outcropCorner2 = - halfOutcropPlaneStrike + halfOutcropPlaneDip;
outcropCorner3 = + halfOutcropPlaneStrike + halfOutcropPlaneDip;
outcropCorner4 = + halfOutcropPlaneStrike - halfOutcropPlaneDip;

% __________________________________________
% 7 - rectangular frame for projection plane

% downdip vector with length = 1/2 Height
halfProjectionPlaneDip(1) = cos(projectionPlaneDip*rad)*sin(projectionPlaneDipAzimuth*rad)*outcropHeight/2;  % uses outcrop height
halfProjectionPlaneDip(2) = cos(projectionPlaneDip*rad)*cos(projectionPlaneDipAzimuth*rad)*outcropHeight/2;
halfProjectionPlaneDip(3) = sin(projectionPlaneDip*rad)*outcropHeight/2;  % in the facets version there was a "-" in front of this

% alongstrike vector with length = 1/2 Length
halfProjectionPlaneStrike(1) = sin(projectionPlaneStrike*rad)*outcropLength/2;  % uses outcrop length
halfProjectionPlaneStrike(2) = cos(projectionPlaneStrike*rad)*outcropLength/2;
halfProjectionPlaneStrike(3) = 0;

% four corners of outcrop best fit plane, which is centerd in outcropCenter
projectionPlaneCorner1 = - halfProjectionPlaneStrike - halfProjectionPlaneDip;
projectionPlaneCorner2 = - halfProjectionPlaneStrike + halfProjectionPlaneDip;
projectionPlaneCorner3 = + halfProjectionPlaneStrike + halfProjectionPlaneDip;
projectionPlaneCorner4 = + halfProjectionPlaneStrike - halfProjectionPlaneDip;

% _____________________________________________
% 8 - project lines onto projection plane in 3D

lines_prj_3D = {};
for i = 1:Nlines
    Nnodes = size(lines{i},1);
    lines_prj_3D{i} = zeros(size(lines{i}));

    for j = 1:Nnodes
        lines_prj_3D{i}(j,:) = lines{i}(j,:) + (-lines{i}(j,:) * projectionPlaneNormal') * projectionPlaneNormal;
    end
end

% ___________
% 9 - 3D plot

% initialize figure 1
figure(1); hold on; axis equal; grid on
title({'3D view of projected fracture facets - mean normal in cyan, outcrop in magenta,';...
    'projection plane in green, intersection of outcrop and mean facet plane in red'})
xlabel('East'); ylabel('North'); zlabel('Z');

% plot input lines and centroids
for i = 1:Nlines
    plot3(lines{i}(:,1), lines{i}(:,2), lines{i}(:,3));
    plot3(centroids(:,1), centroids(:,2), centroids(:,3), 'd', MarkerSize=3);
end

% plot projected lines
for i = 1:Nlines
    plot3(lines_prj_3D{i}(:,1), lines_prj_3D{i}(:,2), lines_prj_3D{i}(:,3));
end

% plot normal unit vectors scaled by outcrop length
quiver3(0,0,0,meanNormal(1)*outcropLength/3,meanNormal(2)*outcropLength/3,meanNormal(3)*outcropLength/3,'LineWidth',2,'Color','cyan');
quiver3(0,0,0,outcropNormal(1)*outcropLength/3,outcropNormal(2)*outcropLength/3,outcropNormal(3)*outcropLength/3,'LineWidth',2,'Color','magenta');
quiver3(0,0,0,projectionPlaneNormal(1)*outcropLength/3,projectionPlaneNormal(2)*outcropLength/3,projectionPlaneNormal(3)*outcropLength/3,'LineWidth',2,'Color','green');
quiver3(0,0,0,facetOutcropIntersection(1)*outcropLength/3,facetOutcropIntersection(2)*outcropLength/3,facetOutcropIntersection(3)*outcropLength/3,'LineWidth',2,'Color','red');

% plot frame for average outcrop plane
plot3([outcropCorner1(1) outcropCorner2(1) outcropCorner3(1) outcropCorner4(1) outcropCorner1(1)]',...
    [outcropCorner1(2) outcropCorner2(2) outcropCorner3(2) outcropCorner4(2) outcropCorner1(2)]',...
    [outcropCorner1(3) outcropCorner2(3) outcropCorner3(3) outcropCorner4(3) outcropCorner1(3)]',...
    'LineWidth',2,'Color','magenta');

% plot frame for projection plane
plot3([projectionPlaneCorner1(1) projectionPlaneCorner2(1) projectionPlaneCorner3(1) projectionPlaneCorner4(1) projectionPlaneCorner1(1)]',...
    [projectionPlaneCorner1(2) projectionPlaneCorner2(2) projectionPlaneCorner3(2) projectionPlaneCorner4(2) projectionPlaneCorner1(2)]',...
    [projectionPlaneCorner1(3) projectionPlaneCorner2(3) projectionPlaneCorner3(3) projectionPlaneCorner4(3) projectionPlaneCorner1(3)]',...
    'LineWidth',2,'Color','green');

% _____________________________
% 10 - projected lines 3D -> 2D

lines_prj_2D = {};

% U coorddinate along scanline unit vector
% V coordinate along intersection unit vector
for i = 1:Nlines
    Nnodes = size(lines{i},1);
    lines_prj_2D{i} = zeros(size(lines{i},1), 2);

    for j = 1:Nnodes
        lines_prj_2D{i}(j,1) = dot(lines_prj_3D{i}(j,:),meanNormal);
        lines_prj_2D{i}(j,2) = dot(lines_prj_3D{i}(j,:),facetOutcropIntersection);
    end
end

% _______________________________
% 10 - plot projected lines in 2D
figure(2); hold on; axis equal; grid on
%set(gca,'XDir','reverse'); set(gca,'YDir','reverse')
title('2D view of projected fracture facets')
xlabel('mean normal axis (U)'); ylabel('outcrop - mean facet intersection (V)');

for i = 1:Nlines
    plot(lines_prj_2D{i}(:,1), lines_prj_2D{i}(:,2),'LineWidth',2);
end

% _________________________________________________________________________
% 11 - input number of scanlines and find intersections for each facet trace

% input number of scanlines
disp(' ');
Nscan = -1;
while Nscan<=0
    Nscan = input('input number of scanlines [100]: ');
    if isempty(Nscan), Nscan = 100; end
    Nscan = round(Nscan);
end
disp(' ');

% __________________
% 12 - scanline area
% scanline max and min U are given by dataset max and min U plus some tolerance
maxUs = zeros(Nlines);
minUs = zeros(Nlines);
maxVs = zeros(Nlines);
minVs = zeros(Nlines);

for i = 1:Nlines
    maxUs(i) = max(lines_prj_2D{i}(:,1));
    minUs(i) = min(lines_prj_2D{i}(:,1));
    maxVs(i) = max(lines_prj_2D{i}(:,2));
    minVs(i) = min(lines_prj_2D{i}(:,2));
end

maxU = max(maxUs);
minU = min(minUs);
maxV = max(maxVs);
minV = min(minVs);

scanMaxU = maxU + outcropLength*0.005;
scanMinU = minU - outcropLength*0.005;
scanMaxV = maxV - outcropHeight*0.005;
scanMinV = minV + outcropHeight*0.005;

% scanline step along V is given dividing maxV - minV minus some tolerance by number of scanlines - 1
scanStep = (scanMaxV - scanMinV)/(Nscan-1);

% _____________________________
% 14 - intersection and spacing
intersectionCount = zeros(1,Nscan);
spacing = [];

for i = 1:Nscan
    scanV = scanMinV + scanStep*(i-1);
    plot([scanMaxU scanMinU],[scanV scanV],'Color',[0.5 0.5 0.5]);
    recordedIntersection = [];

    for j = 1:Nlines
        [intersectionU,intersectionV] = polyxpoly(lines_prj_2D{j}(:,1), lines_prj_2D{j}(:,2), [scanMaxU scanMinU], [scanV scanV]);
        if isfinite(intersectionU)
            plot(intersectionU,intersectionV,'k.');
            intersectionU = intersectionU';
            % intersectionV = intersectionV';  % this is used just for plotting
            intersectionCount(i) = intersectionCount(i)+length(intersectionU);            
            recordedIntersection = [recordedIntersection intersectionU];
        end
    end

    if intersectionCount(i)>1
        disp(['Scanline ' num2str(i) ' - ' num2str(intersectionCount(i)) ' intersections']);
        recordedIntersection = sort(recordedIntersection);
        thisScanSpacing = recordedIntersection(2:end) - recordedIntersection(1:end-1);
        %disp(num2str(thisScanSpacing));
        spacing = [spacing thisScanSpacing];
    end
end

spacingPrctile = prctile(spacing,0:5:100);
spacingMean = mean(spacing);
spacingMode = mode(spacing);
spacingStDev = std(spacing);

% ___________________________________
% 15 - fit distributions and K-S test
% fit distributions
LognormDist = fitdist(spacing','Lognormal');
NormDist = fitdist(spacing','Normal');
ExpDist = fitdist(spacing','Exponential');
BurrDist = fitdist(spacing','Burr');
GammaDist = fitdist(spacing','Gamma');
LogisticDist = fitdist(spacing','Logistic');

% K-S tests

[LognormKSHo,LognormKSPval,LognormKSKstat] = kstest(spacing', LognormDist);
[NormKSHo,NormKSPval,NormKSKstat] = kstest(spacing', NormDist);
[ExpKSHo,ExpKSPval,ExpKSKstat] = kstest(spacing', ExpDist);
[BurrKSHo,BurrKSPval,BurrKSKstat] = kstest(spacing', BurrDist);
[GammaKSHo,GammaKSPval,GammaKSKstat] = kstest(spacing', GammaDist);
[LogisticKSHo,LogisticKSPval,LogisticKSKstat] = kstest(spacing', LogisticDist);
% K-S test text
%Lognorm

Lognorm_txt1 =['K-S test for LogNormal spacing distribution with μ = ' num2str(LognormDist.mu) ' - σ = ' num2str(LognormDist.sigma)];
Lognorm_txt2 = ['K statistics = ' num2str(LognormKSKstat)];
Lognorm_txt3 = ['p-value = ' num2str(LognormKSPval)];
Lognorm_txt4 = ['Ho = ' num2str(LognormKSHo)];
if LognormKSHo == 0
    Lognorm_txt5 = ['No evidence against Ho of no difference with LogNormal distr. -> Ho retained -> a LogNormal length distribution is detected at 5% significance'];
else
    Lognorm_txt5 = ['Strong evidence against Ho of no difference with LogNormal distr. -> Ho rejected -> no LogNormal length distribution is detected at 5% significance'];
end

%Norm

Norm_txt1 =['K-S test for Normal length distribution with Mean = ' num2str(NormDist.mu) ' - St. dev = ' num2str(NormDist.sigma)];
Norm_txt2 = ['K statistics = ' num2str(NormKSKstat)];
Norm_txt3 = ['p-value = ' num2str(NormKSPval)];
Norm_txt4 = ['Ho = ' num2str(NormKSHo)];
if NormKSHo == 0
    Norm_txt5 = ['No evidence against Ho of no difference with Normal distr. -> Ho retained -> a Normal length distribution is detected at 5% significance'];
else
    Norm_txt5 = ['Strong evidence against Ho of no difference with Normal distr. -> Ho rejected -> no Normal length distribution is detected at 5% significance'];
end

%Exponential

Exp_txt1 =['K-S test for Normal length distribution with Mean = ' num2str(ExpDist.mu)];
Exp_txt2 = ['K statistics = ' num2str(ExpKSKstat)];
Exp_txt3 = ['p-value = ' num2str(ExpKSPval)];
Exp_txt4 = ['Ho = ' num2str(ExpKSHo)];
if ExpKSHo == 0
    Exp_txt5 = ['No evidence against Ho of no difference with Exponential distr. -> Ho retained -> a Exponential length distribution is detected at 5% significance'];
else
    Exp_txt5 = ['Strong evidence against Ho of no difference with Exponential distr. -> Ho rejected -> no Exponential length distribution is detected at 5% significance'];
end

%Burr

Burr_txt1 =['K-S test for Burr length distribution with c = ' num2str(BurrDist.c) ' - k = ' num2str(BurrDist.k) '-  alpha =' num2str(BurrDist.alpha)];
Burr_txt2 = ['K statistics = ' num2str(BurrKSKstat)];
Burr_txt3 = ['p-value = ' num2str(BurrKSPval)];
Burr_txt4 = ['Ho = ' num2str(BurrKSHo)];
if BurrKSHo == 0
    Burr_txt5 = ['No evidence against Ho of no difference with Burr distr. -> Ho retained -> a Burr length distribution is detected at 5% significance'];
else
    Burr_txt5 = ['Strong evidence against Ho of no difference with Burr distr. -> Ho rejected -> no Burr length distribution is detected at 5% significance'];
end

%Gamma

Gamma_txt1 =['K-S test for Gamma length distribution with a = ' num2str(GammaDist.a) ' - b = ' num2str(GammaDist.b)];
Gamma_txt2 = ['K statistics = ' num2str(GammaKSKstat)];
Gamma_txt3 = ['p-value = ' num2str(GammaKSPval)];
Gamma_txt4 = ['Ho = ' num2str(GammaKSHo)];
if GammaKSHo == 0
    Gamma_txt5 = ['No evidence against Ho of no difference with Gamma distr. -> Ho retained -> a Gamma length distribution is detected at 5% significance'];
else
    Gamma_txt5 = ['Strong evidence against Ho of no difference with Gamma distr. -> Ho rejected -> no Gamma length distribution is detected at 5% significance'];
end

%Logistic

Logistic_txt1 =['K-S test for Logistic length distribution with μ = ' num2str(LogisticDist.mu) ' - σ = ' num2str(LogisticDist.sigma)];
Logistic_txt2 = ['K statistics = ' num2str(LogisticKSKstat)];
Logistic_txt3 = ['p-value = ' num2str(LogisticKSPval)];
Logistic_txt4 = ['Ho = ' num2str(LogisticKSHo)];
if LogisticKSHo == 0
    Logistic_txt5 = ['No evidence against Ho of no difference with Logistic distr. -> Ho retained -> a Logistic length distribution is detected at 5% significance'];
else
    Logistic_txt5 = ['Strong evidence against Ho of no difference with Logistic distr. -> Ho rejected -> no Logistic length distribution is detected at 5% significance'];
end

% show K-S test text

disp(char(Lognorm_txt1,Lognorm_txt2,Lognorm_txt3,Lognorm_txt4,Lognorm_txt5))
disp(char(Norm_txt1,Norm_txt2,Norm_txt3,Norm_txt4,Norm_txt5))
disp(char(Exp_txt1,Exp_txt2,Exp_txt3,Exp_txt4,Exp_txt5))
disp(char(Burr_txt1,Burr_txt2,Burr_txt3,Burr_txt4,Burr_txt5))
disp(char(Gamma_txt1,Gamma_txt2,Gamma_txt3,Gamma_txt4,Gamma_txt5))
disp(char(Logistic_txt1,Logistic_txt2,Logistic_txt3,Logistic_txt4,Logistic_txt5))


% _______________
% 15b - plot stats
figure(3); hold on; grid on
%histogram(spacing,'BinWidth',1,'Normalization','pdf');
histogram(spacing,'Normalization','pdf');
plot(linspace(min(spacing),max(spacing),20), pdf(LognormDist,linspace(min(spacing),max(spacing),20)),'r',LineWidth=2)
plot(linspace(min(spacing),max(spacing),20), pdf(NormDist,linspace(min(spacing),max(spacing),20)),'g',LineWidth=2)
plot(linspace(min(spacing),max(spacing),20), pdf(ExpDist,linspace(min(spacing),max(spacing),20)),'b',LineWidth=2)
plot(linspace(min(spacing),max(spacing),20), pdf(BurrDist,linspace(min(spacing),max(spacing),20)),'c',LineWidth=2)
plot(linspace(min(spacing),max(spacing),20), pdf(GammaDist,linspace(min(spacing),max(spacing),20)),'m',LineWidth=2)
plot(linspace(min(spacing),max(spacing),20), pdf(LogisticDist,linspace(min(spacing),max(spacing),20)),'y',LineWidth=2)
% xlim([0 14])
box on
legend('Spacing data','Lognormal','Normal','exponential','Burr','Gamma','Logistic','Location','NE');
xlabel('Spacing'); ylabel('Density');
title({['Spacing histogram from ' num2str(Nscan) ' scanlines'],...
    ['Mean = ' num2str(spacingMean) '  Mode = ' num2str(spacingMode) '  St Dev = ' num2str(spacingStDev)],...
    ['Percentiles:'],...
    ['  0% = ' num2str(spacingPrctile(1),3) '   5% = ' num2str(spacingPrctile(2),3) '  10% = ' num2str(spacingPrctile(3),3) '  15% = ' num2str(spacingPrctile(4),3) '  20% = ' num2str(spacingPrctile(5),3)],...
    [' 25% = ' num2str(spacingPrctile(6),3) '  30% = ' num2str(spacingPrctile(7),3) '  35% = ' num2str(spacingPrctile(8),3) '  40% = ' num2str(spacingPrctile(9),3) '  45% = ' num2str(spacingPrctile(10),3)],...
    [' 50% = ' num2str(spacingPrctile(11),3) '  55% = ' num2str(spacingPrctile(12),3) '  60% = ' num2str(spacingPrctile(13),3) '  65% = ' num2str(spacingPrctile(14),3) '  70% = ' num2str(spacingPrctile(15),3)],...
    [' 75% = ' num2str(spacingPrctile(16),3) '  80% = ' num2str(spacingPrctile(17),3) '  85% = ' num2str(spacingPrctile(18),3) '  90% = ' num2str(spacingPrctile(19),3) '  95% = ' num2str(spacingPrctile(20),3)],...
    ['100% = ' num2str(spacingPrctile(21),3) ]});

figure(4); hold on; grid on
%histogram(spacing,'BinWidth',1,'Normalization','cdf');
histogram(spacing,'Normalization','cdf');
plot(linspace(min(spacing),max(spacing),20), cdf(LognormDist,linspace(min(spacing),max(spacing),20)),'r',LineWidth=2)
xlabel('Spacing'); ylabel('Cumulative Frequency');
title({['Cumulative frequency from ' num2str(Nscan) ' scanlines'],...
    ['Mean = ' num2str(spacingMean) '  Mode = ' num2str(spacingMode) '  St Dev = ' num2str(spacingStDev)],...
    ['Percentiles:'],...
    ['  0% = ' num2str(spacingPrctile(1),3) '   5% = ' num2str(spacingPrctile(2),3) '  10% = ' num2str(spacingPrctile(3),3) '  15% = ' num2str(spacingPrctile(4),3) '  20% = ' num2str(spacingPrctile(5),3)],...
    [' 25% = ' num2str(spacingPrctile(6),3) '  30% = ' num2str(spacingPrctile(7),3) '  35% = ' num2str(spacingPrctile(8),3) '  40% = ' num2str(spacingPrctile(9),3) '  45% = ' num2str(spacingPrctile(10),3)],...
    [' 50% = ' num2str(spacingPrctile(11),3) '  55% = ' num2str(spacingPrctile(12),3) '  60% = ' num2str(spacingPrctile(13),3) '  65% = ' num2str(spacingPrctile(14),3) '  70% = ' num2str(spacingPrctile(15),3)],...
    [' 75% = ' num2str(spacingPrctile(16),3) '  80% = ' num2str(spacingPrctile(17),3) '  85% = ' num2str(spacingPrctile(18),3) '  90% = ' num2str(spacingPrctile(19),3) '  95% = ' num2str(spacingPrctile(20),3)],...
    ['100% = ' num2str(spacingPrctile(21),3) ]});

% _______________________________
% 16 - save figures for reporting
for i=1:4
    figure(i);
    %savefig([path '\' file '_fig_' num2str(i) '.fig']);
    saveas(gcf,[path '\' file '_fig_' num2str(i) '.jpg'])
end

% save spacing
spacing = spacing';
save([path '\' file '_spacing.txt'],'spacing','-ascii');

disp('Files saved');

