// Focus particle reference-frame shader.
//
// Draws a lat/long grid plus a contrasting cap on the +Y pole, all in
// object-local space so the markings rotate with the sphere. Without these
// markings a uniformly shaded sphere reads as motionless under axial spin.
//
// Tier color is a uniform so M3 can re-tint per tier (photon / electron /
// baryon) without forking the shader.

shader_type spatial;
render_mode blend_mix, depth_draw_opaque, cull_back, diffuse_burley, specular_schlick_ggx;

uniform vec4 tier_color : source_color = vec4(0.92, 0.86, 0.72, 1.0);
uniform vec4 line_color : source_color = vec4(0.18, 0.14, 0.08, 1.0);
uniform vec4 pole_color : source_color = vec4(0.95, 0.30, 0.20, 1.0);
uniform float metallic : hint_range(0.0, 1.0) = 0.15;
uniform float roughness : hint_range(0.0, 1.0) = 0.35;
uniform int meridian_count : hint_range(1, 24) = 8;
uniform int parallel_count : hint_range(1, 12) = 4;
uniform float line_width : hint_range(0.001, 0.1) = 0.012;
uniform float line_softness : hint_range(0.0, 1.0) = 0.4;
uniform float pole_cap_angle : hint_range(0.0, 1.5) = 0.30;

varying vec3 v_local_pos;

void vertex() {
	v_local_pos = VERTEX;
}

void fragment() {
	vec3 n = normalize(v_local_pos);

	// Polar angle: 0 at +Y, PI at -Y.
	float theta = acos(clamp(n.y, -1.0, 1.0));
	// Longitude around Y: (-PI, PI].
	float phi = atan(n.z, n.x);

	vec3 albedo = tier_color.rgb;

	if (theta < pole_cap_angle) {
		// Pole cap with a soft rim so it doesn't read as a hard sticker.
		float rim = smoothstep(pole_cap_angle, pole_cap_angle - 0.04, theta);
		albedo = mix(tier_color.rgb, pole_color.rgb, rim);
	} else {
		// Meridians: angular distance to nearest, converted to surface arc
		// length via sin(theta) so line width on the surface stays uniform.
		float meridian_step = TAU / float(meridian_count);
		float phi_mod = mod(phi + PI, meridian_step);
		float meridian_angular = min(phi_mod, meridian_step - phi_mod);
		float meridian_dist = meridian_angular * sin(theta);

		// Parallels: evenly spaced in theta between the two poles.
		float parallel_step = PI / float(parallel_count + 1);
		float theta_mod = mod(theta, parallel_step);
		float parallel_dist = min(theta_mod, parallel_step - theta_mod);

		float min_dist = min(meridian_dist, parallel_dist);
		float inner = line_width * (1.0 - line_softness);
		float line_factor = 1.0 - smoothstep(inner, line_width, min_dist);
		albedo = mix(tier_color.rgb, line_color.rgb, line_factor);
	}

	ALBEDO = albedo;
	METALLIC = metallic;
	ROUGHNESS = roughness;
}
