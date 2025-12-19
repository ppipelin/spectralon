function displayBrdfSpectre(theta, phi, tpio, wvl, fr, fig)
	arguments
		theta (1, 1)
		phi (1, 1)
		tpio (:, 4) % Theta_in, Phi_in, Theta_out, Phi_out in degree
		wvl (:, 1)
		fr (:, 1)
		fig = figure()
	end
	idx = tpio(:, 1) == theta & tpio(:, 4) == phi;
	scatter3(tpio(idx, 3), wvl(idx, :) * 1e9, fr(idx, :), 200, ".k")

	xlabel("\theta_{out} (°)");
	ylabel("Wavelength (nm)");
	zlabel("BRDF (sr^{-1})");
	view([-40, 7]);

	title("\Delta_\phi=" + phi + "° with \theta_{in} = " + theta + "°");
end