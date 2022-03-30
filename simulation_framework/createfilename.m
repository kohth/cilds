function FileName = createfilename(RunParam,runIdx,iTau)

tauHeader = sprintf('sim%.03d_tau%.02d_',runIdx,iTau);
FileName.spk = strcat('./',RunParam.FULLDATAFOLDER,'/',tauHeader,'spiketrain.mat');
FileName.sampledlatent = strcat('./',RunParam.SAMPLEDDATAFOLDER,'/',tauHeader,'sampledlatent.mat');
FileName.binnedspk = strcat('./',RunParam.SAMPLEDDATAFOLDER,'/',tauHeader,'binned_spiketrain.mat');
FileName.binnedfr = strcat('./',RunParam.FULLDATAFOLDER,'/',tauHeader,'binned_firingrate.mat');
FileName.latent = strcat('./',RunParam.FULLDATAFOLDER,'/',tauHeader,'latent.mat');
FileName.binnedlatent = strcat('./',RunParam.SAMPLEDDATAFOLDER,'/',tauHeader,'binnedlatent.mat');
FileName.fluo = cell(1,size(RunParam.GAMMA,2)); FileName.sampledfluo = FileName.fluo;
FileName.decon = cell(1,size(RunParam.GAMMA,2));

for gammaIdx = 1:size(RunParam.GAMMA,2)
    gammaHeader = sprintf('sim%.03d_tau%.02d_gamma%d_',runIdx,iTau,gammaIdx);
    FileName.decon{gammaIdx} = strcat('./',RunParam.SAMPLEDDATAFOLDER,'/',gammaHeader,'deconvolved_spiketrain.mat');
    FileName.fluo{gammaIdx} = strcat('./',RunParam.FULLDATAFOLDER,'/',gammaHeader,'fluorescence.mat');
    FileName.sampledfluo{gammaIdx} = strcat('./',RunParam.SAMPLEDDATAFOLDER,'/',gammaHeader,'sampled_fluorescence.mat');
end
end
