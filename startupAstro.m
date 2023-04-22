% Startup script to use the astrodynamics library
close all; clear; clc

upath = userpath;
restoredefaultpath

% Add folders to path
addpath(genpath('classes'))
addpath(genpath('integrators'))
addpath(genpath('spice'))

%Display welcome message
welcome