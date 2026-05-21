shader_type spatial;
render_mode blend_mix, depth_draw_never, cull_back, unshaded;

uniform sampler2D heatmap_texture : filter_linear;
uniform float opacity : hint_range(0.0, 1.0) = 0.7;

varying vec3 v_local_pos;

void vertex() {
	v_local_pos = VERTEX;
}

void fragment() {
	vec3 n = normalize(v_local_pos);
	float lat = asin(clamp(n.y, -1.0, 1.0));
	float lon = atan(n.z, n.x);

	float u = (lon + PI) / TAU;
	float v = (lat + PI * 0.5) / PI;

	vec3 tex = texture(heatmap_texture, vec2(u, v)).rgb;
	float cw  = tex.r;
	float ccw = tex.g;
	float total = cw + ccw;

	vec3 base_color = vec3(0.02, 0.02, 0.05);

	if (total < 0.01) {
		ALBEDO = base_color;
		ALPHA = opacity * 0.4;
	} else {
		// High gamma: only peaks near max are bright
		float t = pow(clamp(total, 0.0, 1.0), 3.0);

		// Heat ramp: near-black -> dark red -> orange -> white
		vec3 color;
		if (t < 0.25) {
			color = mix(vec3(0.03, 0.01, 0.03), vec3(0.6, 0.05, 0.0), t * 4.0);
		} else if (t < 0.5) {
			color = mix(vec3(0.6, 0.05, 0.0), vec3(1.0, 0.5, 0.0), (t - 0.25) * 4.0);
		} else if (t < 0.75) {
			color = mix(vec3(1.0, 0.5, 0.0), vec3(1.0, 0.85, 0.3), (t - 0.5) * 4.0);
		} else {
			color = mix(vec3(1.0, 0.85, 0.3), vec3(1.0, 1.0, 0.95), (t - 0.75) * 4.0);
		}

		// Subtle CW/CCW tint shift — gated by signal strength to
		// prevent blue bleed at low intensities near the poles
		float ratio = (cw - ccw) / max(total, 0.001);
		color += vec3(-0.08, -0.04, 0.12) * ratio * smoothstep(0.0, 0.05, t);
		color = max(color, vec3(0.0));

		ALBEDO = color;
		ALPHA = mix(opacity * 0.4, min(opacity, 0.9), sqrt(t));
	}
}
