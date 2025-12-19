function displayBrdfAll(tpio, fr, fig)
	arguments
		tpio (:, 4) % Theta_in, Phi_in, Theta_out, Phi_out in degree
		fr (:, 1)
		fig = figure()
	end

	scatter3(tpio(:, 3), tpio(:, 1), abs(wrapTo180(tpio(:, 4) - tpio(:, 2))), 100, fr, ".");

	set(gca, "ColorScale", "log");

	clim([0.2, 3.5]);
	cb = colorbar("XTick", [0.2, 0.3, 0.5, 0.7, 1, 2, 3, 3.5]);
	cb.Label.String = "BRDF (sr^{-1})";

	view([-15, 20]);

	xlim([0, 90]);
	ylim([0, 90]);
	% zlim([0, 360]);

	xticks([0:20:80, 90]);
	yticks([0:20:80, 90]);
	% zticks(0:30:360);

	xlabel("\theta_{out} (°)");
	ylabel("\theta_{in} (°)");
	zlabel("\Delta_{\phi} (°)");

	drawnow();
end