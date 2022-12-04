%% cilds_setup
% script for adding the necessary folders required for running cilds
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : cilds_setup.m
%% LAST CHECKED: 221203 (YYMMDD)

addpath(genpath('core_fa'));
addpath(genpath('core_lds'));
addpath(genpath('core_cilds'));
addpath(genpath('core_cifa'));
addpath(genpath('utilities'));
%addpath(genpath('OASIS_matlab-master')); % You can add a line like that
%after downloading oasis to add oasis to the path
% For simulation framework
addpath(genpath('simulation_framework'));