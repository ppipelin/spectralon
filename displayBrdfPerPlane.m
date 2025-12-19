function displayBrdfPerPlane(tpio, fr, fr_model, name, range, fig)
	arguments
		tpio (:, 4) % Theta_in, Phi_in, Theta_out, Phi_out in degree
		fr (:, 1)
		fr_model (:, 1) = []
		name = ""
		range = 0:45:180 % Plane studied
		fig = figure()
	end
	for deg = range
		nexttile();
		idx_plane = abs(deg - abs(tpio(:, 2) - tpio(:, 4))) < 0.1;

		scatter3(tpio(idx_plane, 3), tpio(idx_plane, 1), fr(idx_plane,:), 200 , ".k"); hold on;

		if (~isempty(fr_model))
			scatter3(tpio(idx_plane, 3), tpio(idx_plane, 1), fr_model(idx_plane, :), 100, "xr"); hold on;
			title(sprintf("Δ_φ=%i°", deg));
			legend(["measure", "model"], "Location", "best");
		else
			title(sprintf("%s with Δ_φ=%i°", name, deg));
		end

		view([-80, 20]);

		xlabel("\theta_{out} (°)");
		ylabel("\theta_{in} (°)");
		zlabel("BRDF (sr^{-1})");
		xticks([0:20:80, 90]);
		yticks([0:20:80, 90]);
		xlim([0, 90]);
		ylim([0, 90]);

		set(gca, "ZScale", "log");
	end
	drawnow();
	hold off;
end