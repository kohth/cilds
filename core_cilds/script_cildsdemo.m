%% script_cildsdemo
% Script cilds demo takes generated data from cilds model, splits
% it into training and testing, and performs several sanity checks
%
% Tests:
% 1. Check that given true parameters, an estimate close to ground truth
%    latents is constructed
% 2. Check that using the posteriors given by the an e-step estimated with
%    ground truth parameters, the m-step parameters are close to ground
%    truth
% 3. Check that the loglikelihood is non-decreasing and approaching the
%    loglikelihood the ground truth parameters give
% 4. Check that given the true parameters, the posteriors and estimated
%    parameters don't shift significantly from ground truth after N
%    iterations
% 5. Check that without the true parameters, the posteriors approach the
%    ground truth latent variables with enough EM iterations (given a
%    reasonable SNR)
% 6. Check that without the true parameters, the estimated parameters
%    approach the ground truth parameters with enough EM iterations (given a
%    reasonable SNR)
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : script_cildsdemo.m
%% LAST CHECKED: 220331 (YYMMDD)
tic
test =5;
rng(1);
% Load data generated from CILDS model
addpath(genpath('../core_fa'));
addpath(genpath('../core_lds'));
addpath(genpath('../oasis_matlab'));
addpath(genpath('../utilities'));
if isfile('cildsgenerateddata.mat')
    load('cildsgenerateddata.mat');
else
    script_generatecildsdata;
    load('cildsgenerateddata.mat');
end
N_NEURON = size(Observation(1).y,1); N_LATENT = size(Latent(1).z,1);
N_TRIAL = size(Observation,2);


trialNo = 1; latentNo = [1,2]; % Specify trial and latent index to plot


% User specified parameter initializations
InitParam.B = diag(ones(N_NEURON,1));
FixParam = struct;

% Split training and testing
RunParam.N_LATENT = N_LATENT; RunParam.N_NEURON = N_NEURON;
RunParam.TRAININD = 1:N_TRIAL/2;
RunParam.TESTIND = ceil(N_TRIAL)/2+1:N_TRIAL;

switch test
    case 1
        fileHeader = '1_cildssanitycheck';
        maxIter = 1;
        % Initialize with true parameters and run 1 expectation step
        [~, TrueResult] = cilds(Observation, RunParam,...
            'InitParam',TrueParam,'maxIter',maxIter,'fileHeader',...
            fileHeader,'initType','fixedInit','splitTrainTest',true);
        TestLatent = Latent(RunParam.TESTIND);
        
        %% Plot results
        figure();
        for iPlot = 1:size(latentNo,2)
            subplot(size(latentNo,2),1,iPlot)
            hold on;
            plot(TestLatent(trialNo).z(latentNo(iPlot),2:end),'k','LineWidth',1.5);
            plot(TrueResult(trialNo).z(latentNo(iPlot),1:end-1),'r','LineWidth',1.5);
            ylabel('latent variable');
        end
        xlabel('time steps');
        subplot(size(latentNo,2),1,1)
        title('latent variables estimated with ground truth parameters');
        disp('test 1 complete');
    case 2
        fileHeader = '2_cildssanitycheck';
        maxIter = 2;
        % Initialize with true parameters and run 1 EM step
        [TrueEstParam] = cilds(Observation, RunParam,...
            'InitParam',TrueParam,'maxIter',maxIter,'FixParam',...
            FixParam,'fileHeader',...
            fileHeader,'initType','fixedInit');
        %% Plot results
        parameter = 'G';
        chosenEstParam = TrueEstParam(end).(parameter);
        chosenTrueParam = TrueParam.(parameter);
        if strcmp(parameter,'P') || strcmp(parameter,'Q') || strcmp(parameter,'R')...
                || strcmp(parameter,'G') || strcmp(parameter,'B')
            chosenEstParam = diag(chosenEstParam);
            chosenTrueParam = diag(chosenTrueParam);
        else
            chosenEstParam = reshape(chosenEstParam,size(chosenEstParam,1)*...
                size(chosenEstParam,2),1);
            chosenTrueParam = reshape(chosenTrueParam,size(chosenTrueParam,1)*...
                size(chosenTrueParam,2),1);
        end
        figure();
        hold on;
        plot(chosenTrueParam,chosenEstParam,'r.','MarkerSize',12);
        h = refline(1,0);
        h.Color = 'k';
        xlabel(sprintf('true %s',parameter));
        ylabel(sprintf('est %s', parameter));
        disp('test 2 complete');
    case 3
        fileHeader = '3_cildssanitycheck_gt';
        % Get ll with ground truth parameters
        maxIter = 1;
        [~, ~,~,truell] = cilds(Observation(RunParam.TRAININD), RunParam,...
            'InitParam',TrueParam,'maxIter',maxIter,'fileHeader',...
            fileHeader,'initType','fixedInit');
        % Get ll with non ground truth parameters
        fileHeader = '3_cildssanitycheck_ngt';
        maxIter = 100;
        [~,~,~,estll] = cilds(Observation, RunParam,...
            'InitParam',InitParam,'initType','randInit','maxIter',maxIter,'FixParam',...
            FixParam,'fileHeader',...
            fileHeader,'splittraintest',true);
        %% Plot result
        estll = estll(estll<0);
        figure();
        hold on;
        plot(estll,'r-');
        h = refline(0,truell(1));
        h.Color = 'k';
        xlabel('iterations'); ylabel('loglikelihood');
        xlim([1 size(estll,1)]);
        disp('test 3 complete');
    case 4
        fileHeader = '4_cildssanitycheck';
        maxIter = 100;
        [TrueEstParam, TrueResult] = cilds(Observation, RunParam,...
            'InitParam',TrueParam,'maxIter',maxIter,'fileHeader',...
            fileHeader,'FixParam',FixParam,...
            'initType','fixedInit','splittraintest',true);
        TestLatent = Latent(RunParam.TESTIND);
        
        %% Plot latent results
        figure();
        for iPlot = 1:size(latentNo,2)
            subplot(size(latentNo,2),1,iPlot)
            hold on;
            plot(TestLatent(trialNo).z(latentNo(iPlot),2:end),'k','LineWidth',1.5);
            plot(TrueResult(trialNo).z(latentNo(iPlot),1:end-1),'r','LineWidth',1.5);
            ylabel('latent variable');
        end
        xlabel('time steps');
        subplot(size(latentNo,2),1,1)
        title('latent variables estimated with ground truth parameters');
        %% Plot parameters
        parameter = 'A';
        chosenEstParam = TrueEstParam(end).(parameter);
        chosenTrueParam = TrueParam.(parameter);
        if strcmp(parameter,'P') || strcmp(parameter,'Q') || strcmp(parameter,'R')...
                || strcmp(parameter,'G') || strcmp(parameter,'B') || strcmp(parameter,'D');
            chosenEstParam = diag(chosenEstParam);
            chosenTrueParam = diag(chosenTrueParam);
        else
            chosenEstParam = reshape(chosenEstParam,size(chosenEstParam,1)*...
                size(chosenEstParam,2),1);
            chosenTrueParam = reshape(chosenTrueParam,size(chosenTrueParam,1)*...
                size(chosenTrueParam,2),1);
        end
        figure();
        hold on;
        plot(chosenTrueParam,chosenEstParam,'r.','MarkerSize',12);
        h = refline(1,0);
        h.Color = 'k';
        xlabel(sprintf('true %s',parameter));
        ylabel(sprintf('est %s', parameter));
        
        disp('test 4 complete');
    case 5
        fileHeader = '5_cildssanitycheck_first';
        maxIter = 1;
        [~, FirstResult,~,~,~,FirstTrainResult] = cilds(Observation, RunParam,...
            'InitParam',InitParam,'maxIter',maxIter,'fileHeader',...
            fileHeader,'FixParam',FixParam,...
            'splittraintest',true,'returnTrain',true,'initType','randInit');
        fileHeader = '5_cildssanitycheck_end';
        maxIter = 100;
        [~, Result,~,~,~,TrainResult] = cilds(Observation, RunParam,...
            'InitParam',InitParam,'maxIter',maxIter,'fileHeader',...
            fileHeader,'FixParam',FixParam,...
            'splittraintest',true,'returnTrain',true,'initType','randInit');
        
        
        % Align estimated latent variables to ground truth latent variables
        TestLatent = Latent(RunParam.TESTIND);
        TrainLatent = Latent(RunParam.TRAININD);
        [Cfirst,Tfirst] = gettransform(TrainLatent,FirstTrainResult);
        [C,T] = gettransform(TrainLatent,TrainResult);
        
        FirstResult = transformlatent(FirstResult,Cfirst,Tfirst);
        Result = transformlatent(Result,C,T);
        
        %% Plot results
        figure();
        for iPlot = 1:size(latentNo,2)
            subplot(size(latentNo,2),1,iPlot)
            hold on;
            plot(TestLatent(trialNo).z(latentNo(iPlot),2:end),'k','LineWidth',1.5);
            plot(FirstResult(trialNo).z(latentNo(iPlot),1:end-1),'b','LineWidth',1.5);
            plot(Result(trialNo).z(latentNo(iPlot),1:end-1),'r','LineWidth',1.5);
            ylabel('latent variable');
        end
        xlabel('time steps');
        subplot(size(latentNo,2),1,1)
        disp('test 5 complete');
    case 6
        fileHeader = '5_cildssanitycheck_first';
        maxIter = 1;
        [FirstEstParam,FirstResult,~,~,~,FirstTrainResult] = cilds(Observation, RunParam,...
            'InitParam',InitParam,'maxIter',maxIter,'fileHeader',...
            fileHeader,'FixParam',FixParam,...
            'splittraintest',true,'returnTrain',true,'initType','randInit');
        fileHeader = '5_cildssanitycheck_end';
        maxIter = 100;
        [EstParam, Result,~,~,~,TrainResult] = cilds(Observation, RunParam,...
            'InitParam',InitParam,'maxIter',maxIter,'fileHeader',...
            fileHeader,'FixParam',FixParam,...
            'splittraintest',true,'returnTrain',true,'initType','randInit');
        
        TestLatent = Latent(RunParam.TESTIND);
        TrainLatent = Latent(RunParam.TRAININD);
        [Cfirst,Tfirst] = gettransform(TrainLatent,FirstTrainResult);
        [C,T] = gettransform(TrainLatent,TrainResult);
        
        FirstEstParam = transformparameter(FirstEstParam,Cfirst,Tfirst);
        EstParam = transformparameter(EstParam,C,T);
        
        %% Plot parameters
        parameter = 'A';
        chosenFirstParam = FirstEstParam(end).(parameter);
        chosenTrueParam = TrueParam.(parameter);
        chosenEstParam = EstParam(end).(parameter);
        if strcmp(parameter,'P') || strcmp(parameter,'Q') || strcmp(parameter,'R')...
                || strcmp(parameter,'G') || strcmp(parameter,'B') || strcmp(parameter,'D');
            chosenEstParam = diag(chosenEstParam);
            chosenTrueParam = diag(chosenTrueParam);
            chosenFirstParam = diag(chosenFirstParam);
        else
            chosenEstParam = reshape(chosenEstParam,size(chosenEstParam,1)*...
                size(chosenEstParam,2),1);
            chosenTrueParam = reshape(chosenTrueParam,size(chosenTrueParam,1)*...
                size(chosenTrueParam,2),1);
            chosenFirstParam = reshape(chosenFirstParam,size(chosenFirstParam,1)*...
                size(chosenFirstParam,2),1);
        end
        figure();
        hold on;
        plot(chosenTrueParam,chosenFirstParam,'b.','MarkerSize',12);
        plot(chosenTrueParam,chosenEstParam,'r.','MarkerSize',12);
        
        h = refline(1,0);
        h.Color = 'k';
        xlabel(sprintf('true %s',parameter));
        ylabel(sprintf('est %s', parameter));
        disp('test 6 complete');
end
toc
%% Functions
% Function gettransform finds the linear transformation that matches the
% ground truth latent using linear regression
function [C,T] = gettransform(TrueLatent,EstLatent)

for iTrial = 1:size(TrueLatent,2)
    TrueLatent(iTrial).z = TrueLatent(iTrial).z(:,2:end);
    EstLatent(iTrial).z = EstLatent(iTrial).z(:,1:end-1);
end
c = concatenate(TrueLatent,'c');
estc = concatenate(EstLatent,'c');
z = concatenate(TrueLatent,'z');
estz = concatenate(EstLatent,'z');

C = diag(diag(c*estc'*inv(estc*estc')));
T = z*estz'*inv(estz*estz');
end

% Function transformlatent uses linear transformations to best align
% estimated latent variables to ground truth latent variables
function TransLatent = transformlatent(EstLatent,C,T)
TransLatent = struct;
for iTrial = 1:size(EstLatent,2)
    TransLatent(iTrial).z = T*EstLatent(iTrial).z;
    TransLatent(iTrial).c = C*EstLatent(iTrial).c;
end
end


% Function transformparameter uses linear transformations to best align
% estimated parameters to ground truth parameters
function TransParam = transformparameter(EstParam,C,T)
TransParam.B = EstParam.B*inv(C);
TransParam.G = C*EstParam.G*inv(C);
TransParam.b = C*EstParam.b;
TransParam.R = EstParam.R;
TransParam.Q = C*EstParam.Q*C';
TransParam.A = C*EstParam.A*inv(T);
TransParam.mu_1 = C*EstParam.mu_1;
TransParam.P = T*EstParam.P*T';
TransParam.D = T*EstParam.D*inv(T);
end