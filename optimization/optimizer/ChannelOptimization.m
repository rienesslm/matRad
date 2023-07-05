clear all
close all
clc
%GO = matRad_GeneticOptimizer();
%BAO = beamAngleOptimizer();
results = [];
% for beamNumber=1:5
%     BAO = beamAngleOptimizer('LIVER_PLAN_HEART_&SPINAL.mat', 10, 5, beamNumber);
%     results(beamNumber) = struct(BAO);
% end
%results
COB = matRad_channelOptimizerBT('PROSTATE.mat', 10, 5);