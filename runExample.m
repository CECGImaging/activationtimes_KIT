addpath('mesh');

%% Example data: Simulated and inversely reconstructed TMVs in simple AFib

Fs = 500; % sampling frequency in Hz
load('exampleData/geo.mat'); % triangle mesh

% true transmembrane voltages
load('exampleData/trueTmv.mat');
trueLat = latDetection(trueTmv, Fs, geo.vertices, geo.faces);

% reconstructed transmembrane voltages
load('exampleData/recTmv.mat');
recLat = latDetection(recTmv, Fs, geo.vertices, geo.faces);

m = latMetrics(recLat, trueLat);

fprintf('AE:  %.4fms (%.4fms)\n', m.AEmean, m.AEstd)
fprintf('FPR: %.4f%%\n', m.FPR)
fprintf('FNR: %.4f%%\n', m.FNR)
fprintf('DSC: %.4f\n', m.DSC)