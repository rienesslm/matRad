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
% BAO1 = beamAngleOptimizer('LIVER_PLAN_HEART_&SPINAL.mat', 10, 5, 1);
BAO2 = beamAngleOptimizer('LIVER_PLAN_HEART_&SPINAL.mat', 10, 5, 2);
% BAO3 = beamAngleOptimizer('LIVER_PLAN_HEART_&SPINAL.mat', 10, 5, 3);
% BAO4 = beamAngleOptimizer('LIVER_PLAN_HEART_&SPINAL.mat', 10, 5, 4);
% BAO5 = beamAngleOptimizer('LIVER_PLAN_HEART_&SPINAL.mat', 10, 5, 5);