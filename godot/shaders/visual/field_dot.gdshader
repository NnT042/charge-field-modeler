shader_type spatial;
render_mode unshaded, cull_disabled;

void vertex() {
	mat4 mv = VIEW_MATRIX * MODEL_MATRIX;
	float sx = length(MODEL_MATRIX[0].xyz);
	float sy = length(MODEL_MATRIX[1].xyz);
	float sz = length(MODEL_MATRIX[2].xyz);
	mv[0] = vec4(sx, 0.0, 0.0, 0.0);
	mv[1] = vec4(0.0, sy, 0.0, 0.0);
	mv[2] = vec4(0.0, 0.0, sz, 0.0);
	MODELVIEW_MATRIX = mv;
}

void fragment() {
	vec2 uv_c = UV * 2.0 - 1.0;
	if (dot(uv_c, uv_c) > 1.0) {
		discard;
	}
	ALBEDO = COLOR.rgb;
	ALPHA = COLOR.a;
}
