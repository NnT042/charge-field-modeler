shader_type spatial;
render_mode blend_add, depth_draw_never, cull_disabled, unshaded;

uniform vec4  ghost_color    : source_color = vec4(0.35, 0.85, 1.0, 1.0);
uniform float rim_power      : hint_range(1.0, 8.0)   = 3.0;
uniform float equator_width  : hint_range(0.005, 0.15) = 0.035;
uniform float pole_cap_angle : hint_range(0.02, 0.5)   = 0.13;
uniform int   meridian_count : hint_range(0, 12)       = 4;
uniform float meridian_width : hint_range(0.001, 0.05) = 0.012;

varying vec3 v_local_pos;

void vertex() {
	v_local_pos = VERTEX;
}

void fragment() {
	vec3  n     = normalize(v_local_pos);
	float theta = acos(clamp(n.y, -1.0, 1.0));
	float phi   = atan(n.z, n.x);

	vec3 base = ghost_color.rgb;

	float ndotv = clamp(dot(NORMAL, VIEW), 0.0, 1.0);
	float rim   = pow(1.0 - ndotv, rim_power) * 0.20;

	float eq_dist = abs(theta - PI * 0.5);
	float eq      = (1.0 - smoothstep(0.0, equator_width, eq_dist)) * 0.90;

	float pole_dist = min(theta, PI - theta);
	float pole      = (1.0 - smoothstep(0.0, pole_cap_angle, pole_dist)) * 0.85;

	float meridian = 0.0;
	if (meridian_count > 0 && pole < 0.05) {
		float mstep = TAU / float(meridian_count);
		float phi_m = mod(phi + PI, mstep);
		float mang  = min(phi_m, mstep - phi_m);
		float mdist = mang * sin(theta);
		meridian = (1.0 - smoothstep(0.0, meridian_width, mdist)) * 0.25;
	}

	float alpha = max(rim, max(eq, max(pole, meridian)));

	vec3 rim_col  = base * 0.6;
	vec3 eq_col   = base;
	vec3 pole_col = mix(base, vec3(1.0), 0.4);

	vec3 col = mix(rim_col, eq_col, smoothstep(0.0, 0.05, eq));
	col = mix(col, pole_col, smoothstep(0.0, 0.05, pole));

	ALBEDO = col;
	ALPHA  = alpha;
}
