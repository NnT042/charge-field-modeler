// Ghost sphere: transparent orbital-shell reference. The sphere is oriented
// so its +Y pole aligns with the outermost spin's rotation axis, making the
// equatorial ring sit exactly on that level's orbital plane.
//
// Elements (in priority order for visual hierarchy):
//   pole caps   — orange; mark the rotation axis endpoints
//   equator     — bright cyan; the outermost orbit circle
//   meridians   — faint cyan; reinforce axis orientation
//   fresnel rim — faint blue; shows the sphere boundary at grazing angles

shader_type spatial;
render_mode blend_add, depth_draw_never, cull_disabled, unshaded;

uniform vec4  rim_color      : source_color = vec4(0.2,  0.5,  1.0,  0.20);
uniform vec4  equator_color  : source_color = vec4(0.35, 0.85, 1.0,  0.90);
uniform vec4  pole_color     : source_color = vec4(0.95, 0.5,  0.15, 0.85);
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
	float theta = acos(clamp(n.y, -1.0, 1.0));  // 0 at +Y pole, PI at -Y
	float phi   = atan(n.z, n.x);               // longitude around Y

	// Fresnel rim — faint base glow at grazing angles.
	float ndotv = clamp(dot(NORMAL, VIEW), 0.0, 1.0);
	float rim   = pow(1.0 - ndotv, rim_power) * rim_color.a;

	// Equatorial ring (theta = PI/2) — the outermost level's orbit circle.
	float eq_dist = abs(theta - PI * 0.5);
	float eq      = (1.0 - smoothstep(0.0, equator_width, eq_dist)) * equator_color.a;

	// Pole caps — mark the rotation axis.
	float pole_dist = min(theta, PI - theta);
	float pole      = (1.0 - smoothstep(0.0, pole_cap_angle, pole_dist)) * pole_color.a;

	// Meridians — suppress on pole caps to keep them clean.
	float meridian = 0.0;
	if (meridian_count > 0 && pole < 0.05) {
		float mstep = TAU / float(meridian_count);
		float phi_m = mod(phi + PI, mstep);
		float mang  = min(phi_m, mstep - phi_m);
		float mdist = mang * sin(theta);
		meridian = (1.0 - smoothstep(0.0, meridian_width, mdist)) * 0.25;
	}

	// Compose alpha: each feature overrides the previous.
	float alpha = max(rim, max(eq, max(pole, meridian)));

	// Color: pole > equator > meridian/rim blend.
	vec3 col = mix(rim_color.rgb, equator_color.rgb, smoothstep(0.0, 0.05, eq));
	col = mix(col, pole_color.rgb, smoothstep(0.0, 0.05, pole));

	ALBEDO = col;
	ALPHA  = alpha;
}
