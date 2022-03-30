%% script_simstat
% script_simstat loads results of dimensionality reduction in simulation
% and computes R^2 between estimated latent variables and ground truth
% latent variables to produce Fig. 3d.
%
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : script_simstat.m
%% LAST CHECKED: 220330 (YYMMDD)

allRedmeth = {'cilds'};
dataFolder = 'sim_datasample';
resultFolder = "sim_results";
statFolder = "sim_stat";
if ~exist(statFolder)
    mkdir(statFolder);
end

numFold = 2; % Number of cross-validation folds
simIdx = 6; gammaIdx = 1; tauIdx = 1:6;

for iSim = simIdx
    for jGamma = gammaIdx
        statFile = sprintf("./%s/sim%.03d_gamma%d_latentrsquared.mat",...
            statFolder,simIdx,gammaIdx);
        if ~exist(statFile)
            FileArr=dir(strcat('./',dataFolder,sprintf('/sim%.03d*.mat',iSim)));
            load(strcat(FileArr(1).folder,'\',FileArr(1).name));
            for kMeth = 1:size(allRedmeth,2)
                redmeth = allRedmeth{kMeth};
                for lTau =tauIdx
                    % load ground truth latent variables
                    dataFile = sprintf("./%s/sim%.03d_tau%.02d_sampledlatent.mat",...
                        dataFolder,iSim,lTau);
                    load(dataFile);
                    % load dim red permutation indices
                    permIndFile = sprintf("./%s/sim%.03d_tau%.02d_gamma%d_%s_permutedindices.mat",...
                        resultFolder,iSim,lTau,jGamma,redmeth);
                    load(permIndFile);
                    TrueLatent = data(ind); %rearrange data in the same way as dimred cv
                    % load and combine estimated latent variables
                    for mFold = 1:numFold
                        resultFile = sprintf("./%s/sim%.03d_tau%.02d_gamma%d_010_%d_%s_result_test.mat",...
                            resultFolder,iSim,lTau,jGamma,mFold,redmeth);
                        load(resultFile);
                        EstLatent = data;
                        
                        if strcmp(redmeth,'cilds')
                            trueTime = 2:size(TrueLatent(1).z,2);
                            estTime = 1:size(EstLatent(1).z,2); %cilds starts at z_2
                        else
                            trueTime = 1:size(TrueLatent(1).z,2)-1;
                            estTime = 1:size(EstLatent(1).z,2)-1;
                        end
                        
                        % Align estimated latent variables to ground truth
                        % latent variables using linear regression
                        TransLatent((mFold-1)*100+1:mFold*100) = transformlatent(TrueLatent, trueTime, EstLatent,estTime, mFold);
                    end
                    
                    % compute r-square between true latent variable and
                    % estimated latent variable
                    stat(kMeth,lTau,:) = computersquared(TransLatent,'ztr',estTime,TrueLatent,'z',trueTime);
                    
                end
            end
            save(statFile,'stat','statType','allRedmeth');
        else
            load(statFile);
        end
    end
end

% plot results in similar format to fig 3 panel d
plotData(stat,allRedmeth);

%% === Begin functions ====

% Function calculatetrans finds best transformation using linear regression
function T = calculatetrans(Y,yTime,X,xTime)
for i = 1:size(Y,2)
    Y(i).z = Y(i).z(:,yTime);
    X(i).z = X(i).z(:,xTime);
end
y = [Y(:).z];
x = [X(:).z];
T = (y*x')*inv(x*x');
end

% Function transformlatent aligns the estimated latent variables to the
% true latent variables
function TransLatent = transformlatent(TrueLatent,trueTime, EstLatent,estTime,mFold)

CurrTrue = TrueLatent((mFold-1)*100+1:mFold*100);
foldDiv = floor(linspace(1, size((mFold-1)*100+1:mFold*100,2)+1, 2+1));
for nFold = 1:2 % get transformation T from one half and test on other half
    testMask = false(1, size((nFold-1)*100+1:nFold*100,2));
    testMask(foldDiv(nFold):foldDiv(nFold+1)-1) = true;
    trainMask = ~testMask;
    TrainLatent = EstLatent(trainMask);
    TestLatent = EstLatent(testMask);
    TrainTrue = CurrTrue(trainMask);
    Trans(nFold).T = calculatetrans(TrainTrue,trueTime,TrainLatent,...
        estTime);
    for iTrial = 1:size(TestLatent,2)
        CurrTransLatent(iTrial).ztr = Trans(nFold).T*TestLatent(iTrial).z(1:size(Trans(nFold).T,1),:);
    end
    TransLatent(testMask)=CurrTransLatent;
end
end

% Function computersquared computes the R^2 between the estimated latent
% variable and true latent variable
function rSquared = computersquared(EstData,estField,estTime,TrueData,trueField,trueTime)
for i = 1:size(EstData,2)
    EstData(i).(estField) = EstData(i).(estField)(:,estTime);
    TrueData(i).(trueField) = TrueData(i).(trueField)(:,trueTime);
    estLat = EstData(i).(estField);
    trueLat = TrueData(i).(trueField);
    SSR(i,:) = sum((trueLat-estLat).^2,2);
    SST(i,:) = sum((trueLat-mean(trueLat,2)).^2,2);
    SSE(i,:) = sum((estLat-mean(trueLat,2)).^2,2);
    rSquared = 1-SSR./SST;
end
% rSquared = SSE./SST;
% rSquared = mean(rSquared,1);
rSquared = reshape(rSquared,1,size(rSquared,1)*size(rSquared,2));

end

% Function plotdata takes R^2 and plots it
function plotData(stat,allRedmeth)
%% plot settings
set(0,'DefaultAxesTitleFontWeight','bold');
set(0,'defaultfigurecolor','w');
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultAxesFontSize',8);
traces = figure('Renderer', 'painters', 'units','centimeters','Position', [0 0 25 8],'PaperUnits','inches');

ldsColor = [0,206,209]./(255);
lds2stageColor = [138,43,226]./(255);
cildsColor = [255, 165, 0]./(255);
cifaColor = cildsColor./2;
allColors = [cildsColor;lds2stageColor;ldsColor];
yLimVal = [0 1];
methPosition = [-0.1, 0, 0.1];

subplot(1,5,1:2);
hold on;
for iMeth = 1:size(allRedmeth,2)
    xTau = 1:6;
    currPos = methPosition(iMeth);
    currStat = squeeze(stat(iMeth,:,:));
    meanStat = mean(currStat,2);
    stdStat = std(currStat,0,2);
    semStat = std(currStat,0,2)./sqrt(size(currStat,2));
    a=plot(xTau+currPos,meanStat,'o-','MarkerFaceColor',allColors(iMeth,:),'MarkerSize',6,'LineWidth',1.5);
    a.Color = allColors(iMeth,:);
    a =errorbar(xTau+currPos,meanStat,stdStat,'LineStyle','none','Color',allColors(iMeth,:),'CapSize',0,'LineWidth',1.5);
    b =errorbar(xTau+currPos,meanStat,semStat,'LineStyle','none','Color','k','CapSize',0,'LineWidth',1.5);
    %     plot(xTau(tLatent)+currPos,meanStat(tLatent),'o','MarkerEdgeColor','k','MarkerSize',8,'LineWidth',2);
end
xlim([0,7])
ylim(yLimVal);
ylabel('R^2','FontWeight','bold');
xticklabels({'','50','100','200','1000','2000','5000'});
set(gca,'TickDir','out');
xlabel('Latent timescale (ms)','FontWeight','bold');
end