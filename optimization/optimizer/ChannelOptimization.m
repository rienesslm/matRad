clear all
close all
clc
% 8 Channels
pln.propStf.template.activeNeedles8 =[0 0	0 0	0 0	0 0	0 0	0 0	0
                                      0 0 0 0 0 0 0 0 0 0 1 0 0
                                      0 0 0 0	0 0	0 0	0 0	0 0	0
                                      1 0	0 0	0 0	0 0	0 0	0 0	0
                                      0 0	0 0	0 0	0 0	0 0	0 0	0
                                      0 0 0 0 0 0 0 1 0 0	0 0	0
                                      0 0	0 0	0 0	0 0	0 0	0 0	1
                                      0 0	0 0	0 0	0 0	0 0	0 0	0
                                      0 0	0 0	0 0	0 0	0 0	0 0	0
                                      0 0	0 0	0 0	1 0	0 0	0 0	0
                                      1 0	0 0	0 0	0 0	0 0	0 0	0
                                      0 0	0 0	0 0	0 0	1 0	0 0	0
                                      0 0	0 0	0 0	1 0	0 0	0 0	0];

% 12 Channels
pln.propStf.template.activeNeedles12 =       [0	0	0	0	0	0	0	0	1	0	0	0	0
                                            0	0	0	1	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	0	0	1	0	0	0
                                            0	0	1	0	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	1	0	0	0	0	0
                                            0	0	0	0	0	0	0	0	0	0	1	0	0
                                            0	1	0	0	0	0	0	0	0	0	0	0	1
                                            0	0	0	0	0	0	1	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	0	0	0	0	0	0
                                            0	0	0	1	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	0	0	0	1	0	0
                                            1	0	0	0	0	0	0	0	0	0	0	0	0];

% 16 Channels
pln.propStf.template.activeNeedles16 =       [0	0	0	0	0	0	0	1	0	0	0	0	0
                                            0	1	0	0	0	0	0	0	0	0	1	0	0
                                            0	0	0	1	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	1	0	0	0	0	0
                                            0	0	1	0	0	0	0	0	0	0	0	1	0
                                            0	0	0	0	0	0	0	0	0	0	0	0	0
                                            1	0	0	0	0	0	1	0	0	0	0	0	0
                                            0	0	0	1	0	0	0	0	0	0	1	0	0
                                            0	0	0	0	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	0	1	0	0	0	0
                                            0	0	0	1	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	0	0	0	0	0	1
                                            0	0	0	0	0	1	0	0	0	0	0	1	0];

% 20 Channels
pln.propStf.template.activeNeedles20 =       [0	0	0	0	0	0	1	0	0	0	1	0	0
                                            0	0	1	0	0	1	0	0	0	0	0	0	0
                                            1	0	0	0	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	0	0	1	0	0	0
                                            0	0	0	0	1	0	1	0	0	0	0	0	0
                                            0	1	0	0	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	0	1	0	0	0	0	1
                                            0	0	1	0	0	0	0	0	0	0	1	0	0
                                            0	0	0	0	0	0	0	0	0	0	0	0	0
                                            0	0	0	0	0	0	1	0	0	0	0	0	0
                                            0	0	0	0	1	0	0	0	0	0	1	0	0
                                            0	0	0	0	0	0	0	0	0	1	0	0	1
                                            0	0	0	1	0	0	0	1	0	0	0	0	0];

vector = [1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
pln.propStf.template.activeNeedles8 = reshape(vector, 13, 13);
%GO = matRad_GeneticOptimizer();
%BAO = beamAngleOptimizer();
results = [];
% for beamNumber=1:5
%     BAO = beamAngleOptimizer('LIVER_PLAN_HEART_&SPINAL.mat', 10, 5, beamNumber);
%     results(beamNumber) = struct(BAO);
% end
%results

COB20GS = matRad_channelOptimizerBT('PROSTATE.mat', 4, 5, 8, 8, pln.propStf.template.activeNeedles8);






