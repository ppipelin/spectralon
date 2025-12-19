classdef spectralon_pipelin
	properties (Constant)
		multiplier = 1e3;

		% alpha_D1, gamma_D1,
		% alpha_D2, beta_D2, gamma_D2, alpha_D3, beta_D3, gamma_D3
		% alpha_F1, beta_F1, gamma_F1, alpha_F2, beta_F2
		% alpha_S1, beta_S1, gamma_S1, alpha_S2, alpha_S3
		% alpha_B1, beta_B1, gamma_B1, alpha_B3, alpha_B4
		% alpha_B2, beta_B2, gamma_B2, alpha_b5, alpha_b6
		% alpha_R1, beta_R1, alpha_R2, beta_R2, gamma_R2
		% alpha_G

		parameters = [
			0.0516378927613567, 2.20549117593996, ...
			0.0419179263204508, 0.629300637825107, 1.46109169536381, 1.00967392110707, 0.00457513407145141, 74.8327340013499, ...
			0.996100780026976, 0.00415387492729943, 2.52733476281762, 0.988456562591691, 0.00557934337130836, ...
			0.989046505426408, 0.00419986853256403, 0.951862172307091, 0.633736041498913, 0.988845900884916, ...
			1.02111733813552, 0.00326156497684221, 0.360428615830561, 0.668181533412052, 0.738325373010976, ...
			1.00547292164243, 0.00601525015516193, 0.478958923244493, 0.0617051816638074, 0.0550056358111914, ...
			0.000108257545500692, 0.000110180421235873, 0.0520104172696285, 0.00961795315901315, 3.85914556039680, ...
			0.332035196771286
			];
	end

	methods (Static)

		function brdf = brdf(tpio, lambda, parameters, delta_phi)
			arguments
				tpio (:, 4) % Theta_in, Phi_in, Theta_out, Phi_out in degree
				lambda (:, 1) % Wavelength in meter
				parameters (:, 34) = spectralon_pipelin.parameters
				delta_phi (:, 1) = abs(wrapTo180(tpio(:,4) - tpio(:,2)))
			end

			brdf = ( ...
				spectralon_pipelin.Dp(parameters, (deg2rad(tpio(:,3) + tpio(:,1)))) + spectralon_pipelin.Df(parameters, deg2rad(tpio(:,1)), deg2rad(tpio(:,3)), deg2rad(delta_phi)) ...
				+ spectralon_pipelin.F(parameters, deg2rad(tpio(:,1)), deg2rad(tpio(:,3)), deg2rad(delta_phi)) .* spectralon_pipelin.R(parameters, lambda*1e6) ...
				+ spectralon_pipelin.S(parameters, deg2rad(tpio(:,1)), deg2rad(tpio(:,3)), deg2rad(delta_phi)) ...
				+ spectralon_pipelin.B(parameters, deg2rad(tpio(:,1)), deg2rad(tpio(:,3)), deg2rad(delta_phi)) .* spectralon_pipelin.R2(parameters, lambda*1e6) ...
				) ...
				.* parameters(end) ...
				.* and(abs(tpio(:,1))<=90,abs(tpio(:,3))<=90);
		end

		function fr = Dp(a, theta)
			alpha_D1 = a(1);
			gamma_D1 = a(2);
			fr = 1 - alpha_D1 .* theta .^ gamma_D1;
		end

		function fr = Dfs(a, theta, delta_phi)
			alpha_D2 = a(3);
			beta_D2 = a(4);
			gamma_D2 = a(5);
			alpha_D3 = a(6);
			beta_D3 = a(7);
			gamma_D3 = a(8);
			fr = (alpha_D2 + beta_D2 .* theta) .^ gamma_D2 .* exp(-(((delta_phi - pi) ./ (alpha_D3 + beta_D3 .* theta) .^ gamma_D3) .^ 2));
		end

		function fr = Df(a, theta_i, theta_o, delta_phi)
			fr = spectralon_pipelin.Dfs(a, theta_i, delta_phi) .* spectralon_pipelin.Dfs(a, theta_o, delta_phi);
		end

		function fr = Fs(a, theta, delta_phi)
			alpha_F1 = a(9);
			beta_F1 = a(10);
			gamma_F1 = a(11)*spectralon_pipelin.multiplier;
			alpha_F2 = a(12);
			beta_F2 = a(13);
			gamma_F2 = 1e3;
			fr = (alpha_F1 + beta_F1 .* theta) .^ gamma_F1 .* (1 + (delta_phi - pi) .^ 2 ./ (alpha_F2 + beta_F2 .* theta) .^ gamma_F2) .^ (-((alpha_F2 + beta_F2 .* theta) .^ gamma_F2 + 1) / 2);
		end

		function fr = F(a, theta_i, theta_o, delta_phi)
			fr = spectralon_pipelin.Fs(a, theta_i, delta_phi) .* spectralon_pipelin.Fs(a, theta_o, delta_phi);
		end

		function fr = S(a, theta_i, theta_o, delta_phi)
			alpha_S1 = a(14);
			beta_S1 = a(15);
			gamma_S1 = a(16)*spectralon_pipelin.multiplier;
			alpha_S2 = a(17);
			alpha_S3 = a(18)*spectralon_pipelin.multiplier;
			fr = (alpha_S1 + beta_S1 .* (theta_i+theta_o)) .^ gamma_S1 .* spectralon_pipelin.gauss(delta_phi - pi, alpha_S2) .* spectralon_pipelin.gauss(theta_i - theta_o, alpha_S3);
		end

		function fr = B(a, theta_i, theta_o, delta_phi)
			alpha_B1 = a(19);
			beta_B1 = a(20);
			gamma_B1 = a(21)*spectralon_pipelin.multiplier;
			alpha_B3 = a(22);
			alpha_B4 = a(23);
			alpha_B2 = a(24);
			beta_B2 = a(25);
			gamma_B2 = a(26)*spectralon_pipelin.multiplier;
			alpha_b5 = a(27);
			alpha_b6 = a(28);
			fr = (alpha_B1 + beta_B1 .* (theta_i+theta_o)) .^ gamma_B1 .* spectralon_pipelin.gauss(delta_phi, alpha_B3) .* spectralon_pipelin.gauss(theta_i - theta_o, alpha_B4) ...
				+ (alpha_B2 + beta_B2 .* (theta_i+theta_o)) .^ gamma_B2 .* spectralon_pipelin.lorentz(delta_phi, alpha_b5) .* spectralon_pipelin.lorentz(theta_i - theta_o, alpha_b6);
		end

		function fr = R(a, lambda)
			alpha_R1 = a(29);
			beta_R1 = a(30);
			fr = alpha_R1 + beta_R1 .* lambda;
		end

		function fr = R2(a, lambda)
			alpha_R2 = a(31);
			beta_R2 = a(32);
			gamma_R2 = a(33);
			fr = (alpha_R2 - beta_R2 .* lambda) .^ gamma_R2;
		end

		function fr = gauss(x, divider)
			fr = exp(-(x ./ divider) .^ 2);
		end

		function fr = lorentz(x, divider)
			fr = 1 ./ (1 + (x ./ divider) .^ 2);
		end
	end
end