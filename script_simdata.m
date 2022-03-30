%% script_simdata
% generate fluorescence from gaussian process latent variables
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : script_simdata.m
%% LAST CHECKED: 220331 (YYMMDD)

tic
rng('default');
script_simparam; % get data generation parameters
paramFolder = "sim_parameter";

%% === Start parallel pool ===
% StartParPool(10);

for simIdx = 4
    %% === load run parameters ===
    load(sprintf('./%s/sim%.03d_runparam.mat',paramFolder,simIdx));
    isSaveFull = false; isSaveSampled = true;
    NTau = 4; MTau = 7; %which latent timescales to save full
    
    tauCombination = repmat(RunParam.TAU,RunParam.N_LATENT,1)'; %create timescales for latents
    
    NLatent = RunParam.N_LATENT; NProj = RunParam.N_PROJ;
    T = RunParam.T_GEN; NTrial = RunParam.N_TRIAL;
    initType = RunParam.INITTYPE;
    
    %% === Create folders ===
    if ~exist(RunParam.PLOTFOLDER)
        mkdir(RunParam.PLOTFOLDER);
    end
    
    if ~exist(RunParam.SAMPLEDDATAFOLDER)
        mkdir(RunParam.SAMPLEDDATAFOLDER);
    end
    
    if ~exist(RunParam.FULLDATAFOLDER)
        mkdir(RunParam.FULLDATAFOLDER);
    end
    
    %% === Begin simulation run ===
    % === Generate model parameters ===
    InitParam = struct;
    InitParam.B = diag(ones(NProj,1));
    InitParam.k = zeros(NLatent,1);
    InitParam.d = zeros(NProj,1);
    InitParam.A = faParams.L;
    InitParam.b = faParams.d;
    GenParam = generate_param(NLatent,NProj,'InitParam',InitParam,'coeffA',...
        RunParam.AVAR,'coeff_b',RunParam.BVAR,'shift_b',RunParam.BSHIFT);
    
    % === Save parameters in file genparam ===
    saveFolder = RunParam.SAMPLEDDATAFOLDER;
    paramHeader = sprintf('sim%.03d_',simIdx);
    paramFile = strcat('./',saveFolder,'/',paramHeader,'genparam.mat');
    GenParam(1:size(tauCombination,1)) = GenParam;
    
%     for iTau =1:size(tauCombination,1) %potentially change to parfor
    for iTau =1 % Just simulate tau = 50ms for example
        % === Selected timescale of latents ===
        tau = tauCombination(iTau,:);
        
        % === Generate data ===
        % === File naming ===
        FileName = createfilename(RunParam,simIdx, iTau);
        
        % === Generate latents and projections ===
        Model = simmodel(RunParam,tau,GenParam(iTau));
        
        % === Generate spike trains ===
        Spktrain = simspktrain(RunParam,Model,'units','s');
        
        % === Bin latents and observations ===
        Binned_spk = struct; Sampled_lat = struct; Binned_fr = struct;
        Binned_lat = struct;
        for jSplit = 1:RunParam.N_SPLIT
            Binned_spk(jSplit).y = bindata(RunParam,Spktrain(jSplit).y,1);
            Binned_lat(jSplit).z = bindata(RunParam,Model(jSplit).z);
            Sampled_lat(jSplit).z = Model(jSplit).z(:,RunParam.BIN_VEC);
            Binned_fr(jSplit).y = bindata(RunParam,Model(jSplit).y);
        end
        
        
        for n = 1:length(RunParam.GAMMA)
            % === Generate fluorescence ===
            [Ftrace] = simfluorescence(RunParam,Spktrain,n);
            
            % === Sample fluorescence ===
            Sampled_ftrace = struct;
            for jSplit = 1:RunParam.N_SPLIT
                Sampled_ftrace(jSplit).y = Ftrace(jSplit).y(:,RunParam.BIN_VEC);
                Sampled_ftrace(jSplit).c = Ftrace(jSplit).c(:,RunParam.BIN_VEC);
            end
            
            % === Split data into trials ===
            Ftrace = splittrial(Ftrace,RunParam.N_TRIAL,'y','c');
            Sampled_ftrace = splittrial(Sampled_ftrace,RunParam.N_TRIAL,'y','c');
            
            % === Save data ===
            if isSaveFull && (iTau == NTau || iTau == MTau)
                parsave(Ftrace(RunParam.N_REMOVE+1:RunParam.N_TRIAL),FileName.fluo{n});
            end
            if isSaveSampled
                parsave(Sampled_ftrace(RunParam.N_REMOVE+1:RunParam.N_TRIAL),FileName.sampledfluo{n});
            end
            
        end
        
        % === Split data into trials ===
        Binned_lat = splittrial(Binned_lat,RunParam.N_TRIAL,'z');
        Sampled_lat = splittrial(Sampled_lat,RunParam.N_TRIAL,'z');
        Binned_spk = splittrial(Binned_spk,RunParam.N_TRIAL,'y');
        Binned_fr = splittrial(Binned_fr,RunParam.N_TRIAL,'y');
        
        
        % === Save data ===
        if isSaveFull && (iTau == NTau || iTau == MTau)
            parsave(Spktrain,FileName.spk);
            parsave(Model,FileName.latent);
        end
        
        if isSaveSampled
            parsave(Binned_lat(RunParam.N_REMOVE+1:RunParam.N_TRIAL),FileName.binnedlatent);
            parsave(Sampled_lat(RunParam.N_REMOVE+1:RunParam.N_TRIAL),FileName.sampledlatent);
            parsave(Binned_spk(RunParam.N_REMOVE+1:RunParam.N_TRIAL),FileName.binnedspk)
            parsave(Binned_fr(RunParam.N_REMOVE+1:RunParam.N_TRIAL), FileName.binnedfr);
        end
        disp(iTau);
    end
    
    RunParam.N_TRIAL = RunParam.N_TRIAL-RunParam.N_REMOVE;
    save(paramFile,'RunParam','tauCombination','GenParam');
    
end
toc


