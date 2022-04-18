CILDS (Calcium Imaging Linear Dynamical System)
=================
Code accompanying the paper "[Dimensionality reduction of calcium-imaged neuronal population activity](https://www.biorxiv.org/content/10.1101/2022.03.11.480682v1)". Runs on MATLAB2018a.

CILDS generative model, see methods for variable definitions and dimensions.
<p align="left">
<img src="figures/cilds_equations.png" width="400" />
</p>

CILDS allows estimates of shared activity among neurons (i.e., the latent variables) to influence the estimates of deconvolved spiking activity, and vice versa. In other words CILDS performs deconvolution for all neurons and dimensionality reduction jointly, in a unified framework. This is in contrast to deconv-LDS, which deconvolves the activity of each neuron independently. With low-dimensional latent variables that are jointly estimated with the model of calcium decay, CILDS is better able to peer through the calcium decay to more clearly identify the shared activity among neurons, as compared to deconv-LDS and LDS applied directly on fluorescence. The differences between models are illustrated in the figure below.

<p align="center">
<img src="figures/fig_2a.PNG" width="600" />
</p>

Getting Started
-----------
### Installation 
add CILDS and necessary folders to the search path of MATLAB. <br />
File location: main

  `>> cilds_setup`

If using deconvolution or 'ldsInit' for CILDS, add OASIS[[1]](#1) function found in oasis_matlab folder to the search path of MATLAB. For more details, look at https://github.com/zhoupc/OASIS_matlab. The demo does not need this. <br />
File location: oasis_matlab

  `>> oasis_setup`

Demo
-----------
### Apply CILDS on CILDS model generated data in core_cilds folder
File location: core_cilds

  `>> script_cilds_demo`

An example of CILDS being used is 

```matlab 
[EstParam, Result] = cilds(data, RunParam);
```
#### *Input*:
* data - Structure containing recorded data (fluorescence traces)
  - Dimensions: 1 x N_TRIAL  
  - Fields: 
    + y (N_NEURON x T) -- neural data

* RunParam - Structure containing dimension of observations (no. neurons) and dimension of latent variables, training indices and testing indices
  - Dimensions: 1 x 1
  - Fields:
    + N_LATENT -- desired latent dimension
    + TRAININD (only if splitting training and testing) -- trial indices of training data
    + TESTIND  (only if splitting training and testing) -- trial indices of testing data

Note that to use 'ldsInit' for initializing CILDS, OASIS needs to be installed.

In the given demo, a few tests are given for running CILDS. For instance, running test 5 checks that without the true parameters, the posteriors approach the ground truth latent variables with enough EM iterations (given a reasonable SNR)

#### *Output (Takes about 30s)*: 
* EstParam - Structure containing estimated model parameters from maximization step
  - Dimensions: 1 x 1
  - Fields: 
    + A (N_NEURON x N_LATENT),
    + B (N_NEURON x N_NEURON),
    + G (N_NEURON x N_NEURON) - Gamma,
    + D (N_LATENT x N_LATENT),
    + Q (N_NEURON x N_NEURON),
    + R (N_NEURON x N_NEURON),
    + P (N_NEURON x N_LATENT),
    + b (N_NEURON x 1),
    + mu_1 (N_NEURON x 1),
    + cov_1 (N_NEURON x N_NEURON)
    + h_2 (N_LATENT x 1),
    + G_2 (N_LATENT x N_LATENT)

* Result - Structure containing estimated posteriors from expectation step
  - Dimensions - 1 x N_TRIAL
  - Fields (If leaveoneout toggled, flProj and frProj saved (predicted fluorescence and firing rate): 
    + z (N_LATENT x T)
    + c (N_NEURON x T)

Sample figure produced, in a manner similar to Fig. 2 of paper (also shown in this document), right panel with estimated latent variables plotted over ground truth latent variables over time: 
<p align="left">
<img src="figures/cilds_check5.PNG" width="400" />
</p>

Additional expected outputs are included in the figures folder under the header cilds_check*. 

Scripts used in paper
-----------
### Simulate data using simulation framework
File location: main

  `>> script_simdata`
  
  
### Run CILDS/CIFA/LDS/deconv-LDS on simulated data (note that this won't run if you haven't run script_simdata)
File location: main

  `>> script_simdimred`
  
This performs a 2-fold cross-validation using the chosen dimensionality reduction method. For example, if using CILDS, the function called would be
```matlab 
cilds_crossvalidate(data,RunParam(iData),'zDimList',RunParam(iData).N_LATENT,...
                    'numFold',2,'maxIter',maxIter,'fileHeader',resultFile,...
                    'initParam',InitParam);
```
The example performs 2-fold cross-validation using N_LATENT number of latent dimensions, a maximum number of EM iterations as specified, saving the results in resultFile and initializing certain parameters as specified in InitParam. 

The two required parameters are data and RunParam, and the remaining parameters specified by '' are optional. data is a structure array of trials containing N_NEURONxTIMEPOINTS of fluorescence data. RunParam is a structure containing information about the data, necessarily specifiying N_LATENT (number of latent dimensions), TRAININD (trials used for training set), and TESTIND (trials used for testing set). 


### Compute R<sup>2</sup> between estimated latent variables and ground truth latent variables in given sample simulated data
Note that the process of data generation and dimensionality reduction takes awhile so the results from the run are provided in sim_stat just for plotting purposes. This step produces something like the figure below. 
File location: main

  `>> script_simstat`  
  
<p align="left">
<img src="figures/fig_3d.PNG" width="400" />
</p>
  

References
-----------
<a id="1">[1]</a> 
Vogelstein, J. T. et al. Fast nonnegative deconvolution for spike train inference from population calcium imaging.
J. Neurophysiol. 104, 3691–3704, 10.1152/jn.01073.2009 (2010). PMID: 20554834, https://doi.org/10.1152/jn.01073.2009.



