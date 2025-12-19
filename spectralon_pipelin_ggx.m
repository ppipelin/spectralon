classdef spectralon_pipelin_ggx
	properties (Constant)
		% roughness, global multiplier, offset, multiplier, exponent, ior

		parameters = [0.9841, 3.8200, 0.1900, 6.4009, 1.5112, 1.1727];
	end

	methods (Static)

		function brdf = brdf(tpio, lambda, parameters, delta_phi)
			arguments
				tpio (:, 4) % Theta_in, Phi_in, Theta_out, Phi_out in degree
				lambda (:, 1) = 550e-9 % Wavelength in meter
				parameters (:, 6) = spectralon_pipelin_ggx.parameters
				delta_phi (:, 1) = abs(wrapTo180(tpio(:,4) - tpio(:,2)))
			end

			roughness = parameters(1);

			wo = spectralon_pipelin_ggx.sphericalToCartesian(tpio(:,4), tpio(:,3));
			wi = spectralon_pipelin_ggx.sphericalToCartesian(tpio(:,2), tpio(:,1));

			wh = wi + wo;
			wh = wh ./ vecnorm(wh, 2, 2);

			dot_wo_wh = dot(wo, wh, 2);
			dot_wi_wh = dot(wi, wh, 2);

			dot_N_wo = wo(:,3);
			dot_N_wi = wi(:,3);
			dot_N_wh = wh(:,3);

			F = spectralon_pipelin_ggx.computeFresnel(wi, wh, parameters(6));

			G1 = spectralon_pipelin_ggx.computeG(roughness, tpio(:,3), dot_wo_wh, dot_N_wo);
			G2 = spectralon_pipelin_ggx.computeG(roughness, tpio(:,1), dot_wi_wh, dot_N_wi);
			G = G1.*G2;

			% tan^2 is sin2 on cos2
			tan2 = (1 - dot_N_wh .^ 2) ./ dot_N_wh .^ 2;
			Xi = double(dot_N_wh > 0);
			D = roughness .^ 2 .* Xi ./ (pi * dot_N_wh .^ 4 .* (roughness .^2 + tan2) .^ 2);

			spec = F .* G .* D ./ (4 * dot_N_wi .* dot_N_wo); % Cook-Torrance use pi but GGX changed to 4

			brdf = parameters(2) .* (parameters(3) + parameters(4) .* spec) .^ parameters(5);
		end

		function o = computeG(roughness, theta, dot_H_w, dot_N_w)
			o = 2 ./ (1 + sqrt(1+roughness.^2.*tand(theta).^2)) .* double((dot_H_w./ dot_N_w)>0);
		end

		function cart = sphericalToCartesian(azimuth, inclination)
			arguments
				azimuth, inclination % In degree
			end
			x = sind(inclination) .* cosd(azimuth);
			y = sind(inclination) .* sind(azimuth);
			z = cosd(inclination);
			cart = cat(2, x, y, z);
		end

		function F = computeFresnel(wi, wh, ior)
			arguments
				wi % Incoming vector
				wh % Half vector
				ior % Index of refraction
			end
			c = abs(dot(wi, wh, 2));
			% g = sqrt(ior .^ 2 / 1 ^ 2 - 1 + c .^ 2);
			g = sqrt(ior .^ 2 - 1 + c.^2);
			F = ((g - c) .^ 2 ./ (2 .* (g + c) .^2 )) ...
				.* (1 + (c .* (g + c) - 1) .^ 2 ./ ((c .* (g - c) + 1) .^ 2));
			F(imag(g) ~= 0) = 1;
		end
	end

end