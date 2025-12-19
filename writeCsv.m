load("spectralon_pipelin.mat", "tpio", "fr", "wvl");

t = array2table(cat(2, tpio, fr), "VariableNames", ["theta_in", "phi_in", "theta_out", "phi_out", string(380:780)]);
writetable(t, "spectralon_pipelin.csv");