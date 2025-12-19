%% Display BRDF
close all;
degree_step = 5;
tpio = combvec(0:degree_step:90, 0, 0:degree_step:90, 0:degree_step:180)';

displayBrdfAll(tpio, spectralon_pipelin.brdf(tpio, 550e-9), figure("Position", [100 600 500 400]));
% displayBrdfAll(tpio, spectralon_pipelin_ggx.brdf(tpio), figure("Position", [100 600 500 400]));

displayBrdfPerPlane(tpio, spectralon_pipelin.brdf(tpio, 550e-9), [], "Spectralon", [180, 0], figure("Position", [100 100 500*2 400]));
% displayBrdfPerPlane(tpio, spectralon_pipelin_ggx.brdf(tpio), [], "Spectralon", [180, 0], figure("Position", [100 100 500*2 400]));

%% Read CSV and compare to measurements
tmp = readmatrix("spectralon_pipelin.csv");
wvl = tmp(1, 5:end);
tpio = tmp(2:end, 1:4);
fr = tmp(2:end, 5:end);
displayBrdfPerPlane(tpio, fr(:, wvl == 550), spectralon_pipelin.brdf(tpio, 550e-9), "Spectralon", [180, 0], figure("Position", [100 100 500*2 400]));
% displayBrdfPerPlane(tpio, fr(:, wvl == 550), spectralon_pipelin_ggx.brdf(tpio), "Spectralon", [180, 0], figure("Position", [100 100 500*2 400]));

%% Compute spectral BRDF
wvl_computed = (380:10:780) * 1e-9;
tpio_lambda = combvec(0:degree_step:90, 0, 0:degree_step:90, 0:degree_step:180, wvl_computed)';

fr_computed = spectralon_pipelin.brdf(tpio_lambda(:, 1:4), tpio_lambda(:, 5));

displayBrdfSpectre(40, 180, tpio_lambda(:, 1:4), tpio_lambda(:, 5), fr_computed, figure("Position", [100 100 500 400]));
displayBrdfSpectre(40, 0, tpio_lambda(:, 1:4), tpio_lambda(:, 5), fr_computed, figure("Position", [600 100 500 400]));
displayBrdfSpectre(90, 180, tpio_lambda(:, 1:4), tpio_lambda(:, 5), fr_computed, figure("Position", [1100 100 500 400]));
displayBrdfSpectre(90, 0, tpio_lambda(:, 1:4), tpio_lambda(:, 5), fr_computed, figure("Position", [1600 100 500 400]));
