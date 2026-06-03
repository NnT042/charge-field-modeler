extends Control
## Real-time line graph of exit-angle distribution.
##
## Shows percentage of all deflections in each 1-degree bin.
## Unsigned mode: 0° (equator) to 90° (pole), 90 bins.
## Signed mode: -90° to +90°, 180 bins.

const MARGIN_L: float = 32.0
const MARGIN_B: float = 20.0
const MARGIN_T: float = 18.0
const MARGIN_R: float = 8.0

var _histogram: PackedFloat32Array = PackedFloat32Array()
var _sample_count: int = 0
var _mean_angle: float = 0.0
var _y_ceiling: float = 1.0
var _hover_pos: Vector2 = Vector2.ZERO
var _hovering: bool = false
var _signed_mode: bool = false

var _line_color: Color = Color(1.0, 0.6, 0.12, 0.9)
var _grid_color: Color = Color(0.25, 0.25, 0.30, 0.5)
var _label_color: Color = Color(0.55, 0.55, 0.6, 1.0)
var _peak_color: Color = Color(1.0, 0.85, 0.3, 0.9)
var _mean_color: Color = Color(0.4, 0.75, 1.0, 0.6)
var _bg_color: Color = Color(0.04, 0.04, 0.07, 0.88)
var _crosshair_color: Color = Color(0.9, 0.9, 0.95, 0.7)

signal export_requested


func _ready() -> void:
	mouse_filter = Control.MOUSE_FILTER_STOP


func update_data(data: PackedFloat32Array, sample_n: int, mean_deg: float) -> void:
	_histogram = data
	_sample_count = sample_n
	_mean_angle = mean_deg
	queue_redraw()


func set_signed_mode(enabled: bool) -> void:
	_signed_mode = enabled
	_y_ceiling = 1.0
	queue_redraw()


func is_signed_mode() -> bool:
	return _signed_mode


func _gui_input(event: InputEvent) -> void:
	if event is InputEventMouseMotion:
		_hover_pos = (event as InputEventMouseMotion).position
		_hovering = true
		queue_redraw()


func _notification(what: int) -> void:
	if what == NOTIFICATION_MOUSE_EXIT:
		_hovering = false
		queue_redraw()


func _snap_ceiling(raw_max: float) -> float:
	var steps: Array = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0]
	for s: float in steps:
		if raw_max <= s:
			return s
	return ceil(raw_max / 10.0) * 10.0


func _draw() -> void:
	draw_rect(Rect2(Vector2.ZERO, size), _bg_color)

	var plot_x: float = MARGIN_L
	var plot_y: float = MARGIN_T
	var plot_w: float = size.x - MARGIN_L - MARGIN_R
	var plot_h: float = size.y - MARGIN_T - MARGIN_B

	var font: Font = ThemeDB.fallback_font
	var fs_sm: int = 9
	var fs_title: int = 10

	var n: int = _histogram.size()
	var expected_n: int = 180 if _signed_mode else 90
	var deg_min: float = -90.0 if _signed_mode else 0.0
	var deg_max: float = 90.0
	var deg_range: float = deg_max - deg_min

	# Title + sample count
	var mode_tag: String = "±90°" if _signed_mode else "|lat|"
	var title: String = "Exit Angle (%s)" % mode_tag
	if _sample_count > 0:
		if _sample_count >= 1000000:
			title += "  (n=%.1fM)" % (float(_sample_count) / 1000000.0)
		elif _sample_count >= 1000:
			title += "  (n=%dK)" % (_sample_count / 1000)
		else:
			title += "  (n=%d)" % _sample_count
	draw_string(font, Vector2(plot_x, plot_y - 5), title,
		HORIZONTAL_ALIGNMENT_LEFT, -1, fs_title, _label_color)

	if n < 2 or n != expected_n or _sample_count < 10:
		draw_string(font, Vector2(plot_x + plot_w * 0.25, plot_y + plot_h * 0.5),
			"Accumulating...", HORIZONTAL_ALIGNMENT_LEFT, -1, fs_sm,
			Color(0.4, 0.4, 0.4))
		return

	# Find actual max
	var raw_max: float = 0.01
	for v: float in _histogram:
		if v > raw_max:
			raw_max = v

	var target: float = _snap_ceiling(raw_max)
	if target > _y_ceiling:
		_y_ceiling = target
	elif target < _y_ceiling * 0.7:
		_y_ceiling = target
	var max_pct: float = _y_ceiling

	# Vertical grid lines at degree marks
	var deg_marks: Array
	if _signed_mode:
		deg_marks = [-90, -60, -30, 0, 30, 60, 90]
	else:
		deg_marks = [0, 15, 30, 45, 60, 75, 90]
	for deg: int in deg_marks:
		var x: float = plot_x + float(deg - int(deg_min)) / deg_range * plot_w
		draw_line(Vector2(x, plot_y), Vector2(x, plot_y + plot_h),
			_grid_color, 1.0)
		var show_label: bool
		if _signed_mode:
			show_label = deg == -90 or deg == -30 or deg == 0 or deg == 30 or deg == 90
		else:
			show_label = deg == 0 or deg == 30 or deg == 60 or deg == 90
		if show_label:
			draw_string(font, Vector2(x - 5, plot_y + plot_h + 14),
				str(deg), HORIZONTAL_ALIGNMENT_LEFT, -1, fs_sm, _label_color)

	# Horizontal grid lines
	var y_steps: int = 4
	for i in range(0, y_steps + 1):
		var frac: float = float(i) / float(y_steps)
		var y: float = plot_y + plot_h * (1.0 - frac)
		draw_line(Vector2(plot_x, y), Vector2(plot_x + plot_w, y),
			_grid_color, 1.0)
		if i > 0:
			var val: float = max_pct * frac
			var lbl: String
			if val >= 1.0:
				lbl = "%.0f%%" % val
			else:
				lbl = "%.1f%%" % val
			draw_string(font, Vector2(1, y + 3),
				lbl, HORIZONTAL_ALIGNMENT_LEFT, -1, fs_sm - 1, _label_color)

	# Mean angle vertical marker
	if _signed_mode:
		# Mean is always positive (absolute), show on both sides
		if _mean_angle > 0.5 and _mean_angle < 89.5:
			var mx_pos: float = plot_x + (_mean_angle + 90.0) / deg_range * plot_w
			var mx_neg: float = plot_x + (-_mean_angle + 90.0) / deg_range * plot_w
			draw_dashed_line(Vector2(mx_pos, plot_y), Vector2(mx_pos, plot_y + plot_h),
				_mean_color, 1.0, 4.0)
			draw_dashed_line(Vector2(mx_neg, plot_y), Vector2(mx_neg, plot_y + plot_h),
				_mean_color, 1.0, 4.0)
	else:
		if _mean_angle > 0.5 and _mean_angle < 89.5:
			var mx: float = plot_x + _mean_angle / deg_range * plot_w
			draw_dashed_line(Vector2(mx, plot_y), Vector2(mx, plot_y + plot_h),
				_mean_color, 1.0, 4.0)

	# Line graph
	var points: PackedVector2Array = PackedVector2Array()
	points.resize(n)
	for i in n:
		var x: float = plot_x + (float(i) + 0.5) / float(n) * plot_w
		var y: float = plot_y + plot_h * (1.0 - _histogram[i] / max_pct)
		points[i] = Vector2(x, clampf(y, plot_y, plot_y + plot_h))

	if points.size() >= 2:
		draw_polyline(points, _line_color, 1.5, true)

	# Peaks
	var sorted_vals: Array = []
	for v: float in _histogram:
		sorted_vals.append(v)
	sorted_vals.sort()
	var median: float = sorted_vals[sorted_vals.size() / 2] if sorted_vals.size() > 0 else 0.0
	var threshold: float = maxf(median * 2.0, max_pct * 0.15)

	var peak_bins: Array = []
	for i in range(2, n - 2):
		if _histogram[i] > threshold:
			if _histogram[i] >= _histogram[i - 1] and _histogram[i] >= _histogram[i + 1]:
				if _histogram[i] >= _histogram[i - 2] and _histogram[i] >= _histogram[i + 2]:
					var dominated: bool = false
					for p: int in peak_bins:
						if absi(i - p) < 6 and _histogram[p] >= _histogram[i]:
							dominated = true
							break
					if not dominated:
						peak_bins.append(i)

	for peak_bin: int in peak_bins:
		var px: float = plot_x + (float(peak_bin) + 0.5) / float(n) * plot_w
		var py: float = plot_y + plot_h * (1.0 - _histogram[peak_bin] / max_pct)
		draw_circle(Vector2(px, py), 3.0, _peak_color)
		var peak_deg: int = peak_bin + int(deg_min) if _signed_mode else peak_bin
		var plbl: String = "%d" % peak_deg
		draw_string(font, Vector2(px + 5, py - 1),
			plbl, HORIZONTAL_ALIGNMENT_LEFT, -1, fs_sm, _peak_color)

	# Plot border
	draw_rect(Rect2(plot_x, plot_y, plot_w, plot_h),
		Color(0.3, 0.3, 0.35, 0.5), false, 1.0)

	# Mouse-over crosshair + readout
	if _hovering and n > 0:
		var mx: float = _hover_pos.x
		var my: float = _hover_pos.y
		if mx >= plot_x and mx <= plot_x + plot_w and my >= plot_y and my <= plot_y + plot_h:
			var frac: float = (mx - plot_x) / plot_w
			var bin_idx: int = clampi(int(frac * float(n)), 0, n - 1)
			var pct_val: float = _histogram[bin_idx]
			var deg_val: int = bin_idx + int(deg_min) if _signed_mode else bin_idx

			draw_line(Vector2(mx, plot_y), Vector2(mx, plot_y + plot_h),
				_crosshair_color, 1.0)

			var val_y: float = plot_y + plot_h * (1.0 - pct_val / max_pct)
			val_y = clampf(val_y, plot_y, plot_y + plot_h)
			draw_line(Vector2(plot_x, val_y), Vector2(plot_x + plot_w, val_y),
				Color(_crosshair_color.r, _crosshair_color.g, _crosshair_color.b, 0.3), 1.0)

			draw_circle(Vector2(mx, val_y), 4.0, _crosshair_color)

			var tip: String = "%d°  %.2f%%" % [deg_val, pct_val]
			var tip_x: float = mx + 8
			if tip_x + 70 > plot_x + plot_w:
				tip_x = mx - 75
			var tip_y: float = val_y - 4
			if tip_y < plot_y + 10:
				tip_y = val_y + 14
			draw_string(font, Vector2(tip_x, tip_y), tip,
				HORIZONTAL_ALIGNMENT_LEFT, -1, fs_sm + 1,
				Color(1.0, 1.0, 1.0, 0.95))
