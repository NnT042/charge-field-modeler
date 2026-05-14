shader_type spatial;
render_mode blend_add, depth_draw_never, cull_back, unshaded;

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

	if (total < 0.01) {
		ALPHA = 0.0;
	} else {
		// Dominant chirality sets the hue; ratio sets saturation
		float dominance = abs(cw - ccw) / total; // 0=even split, 1=pure
		vec3 cw_hue  = vec3(0.2, 0.5, 1.0);  // blue
		vec3 ccw_hue = vec3(1.0, 0.35, 0.15); // red-orange
		vec3 hue = (cw >= ccw) ? cw_hue : ccw_hue;
		vec3 col = mix(vec3(0.7, 0.7, 0.7), hue, dominance) * total;

		ALBEDO = col;
		ALPHA = total * opacity;
	}
}
