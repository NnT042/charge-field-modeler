extends Control
## Real-time line graph of exit-angle distribution.
##
## Shows percentage of all deflections in each 1-degree bin from
## 0 degrees (equatorial exit) to 90 degrees (polar exit). Based on exit
## POSITION latitude (same projection as the heatmap sphere), so histogram
## peaks correspond directly to the bright bands on the heatmap.

const MARGIN_L: float = 32.0
const MARGIN_B: float = 20.0
const MARGIN_T: float = 18.0
const MARGIN_R: float = 8.0

var _histogram: PackedFloat32Array = PackedFloat32Array()
var _sample_count: int = 0
var _mean_angle: float = 0.0
var _y_ceiling: float = 1.0  # smoothed max — only ratchets up, slow decay
var _hover_pos: Vector2 = Vector2.ZERO
var _hovering: bool = false

var _line_color: Color = Color(1.0, 0.6, 0.12, 0.9)
var _grid_color: Color = Color(0.25, 0.25, 0.30, 0.5)
var _label_color: Color = Color(0.55, 0.55, 0.6, 1.0)
var _peak_color: Color = Color(1.0, 0.85, 0.3, 0.9)
var _mean_color: Color = Color(0.4, 0.75, 1.0, 0.6)
var _bg_color: Color = Color(0.04, 0.04, 0.07, 0.88)
var _crosshair_color: Color = Color(0.9, 0.9, 0.95, 0.7)


func _ready() -> void:
	mouse_filter = Control.MOUSE_FILTER_STOP


func update_data(data: PackedFloat32Array, sample_n: int, mean_deg: float) -> void:
	_histogram = data
	_sample_count = sample_n
	_mean_angle = mean_deg
	queue_redraw()


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
	# Snap to a clean value: 0.5, 1, 2, 3, 5, 8, 10, 15, 20, ...
	var steps: Array = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0]
	for s: float in steps:
		if raw_max <= s:
			return s
	return ceil(raw_max / 10.0) * 10.0


func _draw() -> void:
	# Background
	draw_rect(Rect2(Vector2.ZERO, size), _bg_color)

	var plot_x: float = MARGIN_L
	var plot_y: float = MARGIN_T
	var plot_w: float = size.x - MARGIN_L - MARGIN_R
	var plot_h: float = size.y - MARGIN_T - MARGIN_B

	var font: Font = ThemeDB.fallback_font
	var fs_sm: int = 9
	var fs_title: int = 10

	# Title + sample count
	var title: String = "Exit Angle"
	if _sample_count > 0:
		if _sample_count >= 1000000:
			title += "  (n=%.1fM)" % (float(_sample_count) / 1000000.0)
		elif _sample_count >= 1000:
			title += "  (n=%dK)" % (_sample_count / 1000)
		else:
			title += "  (n=%d)" % _sample_count
	draw_string(font, Vector2(plot_x, plot_y - 5), title,
		HORIZONTAL_ALIGNMENT_LEFT, -1, fs_title, _label_color)

	if _histogram.size() < 2 or _sample_count < 10:
		draw_string(font, Vector2(plot_x + plot_w * 0.25, plot_y + plot_h * 0.5),
			"Accumulating...", HORIZONTAL_ALIGNMENT_LEFT, -1, fs_sm,
			Color(0.4, 0.4, 0.4))
		return

	# Find actual max
	var raw_max: float = 0.01
	for v: float in _histogram:
		if v > raw_max:
			raw_max = v

	# Smooth Y ceiling: ratchet up instantly, decay very slowly
	var target: float = _snap_ceiling(raw_max)
	if target > _y_ceiling:
		_y_ceiling = target
	elif target < _y_ceiling * 0.7:
		# Only shrink if actual max is well below ceiling
		_y_ceiling = target
	var max_pct: float = _y_ceiling

	# Vertical grid lines at degree marks
	var deg_marks: Array = [0, 15, 30, 45, 60, 75, 90]
	for deg: int in deg_marks:
		var x: float = plot_x + float(deg) / 90.0 * plot_w
		draw_line(Vector2(x, plot_y), Vector2(x, plot_y + plot_h),
			_grid_color, 1.0)
		if deg == 0 or deg == 30 or deg == 60 or deg == 90:
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
	if _mean_angle > 0.5 and _mean_angle < 89.5:
		var mx: float = plot_x + _mean_angle / 90.0 * plot_w
		draw_dashed_line(Vector2(mx, plot_y), Vector2(mx, plot_y + plot_h),
			_mean_color, 1.0, 4.0)

	# Line graph
	var n: int = _histogram.size()
	var points: PackedVector2Array = PackedVector2Array()
	points.resize(n)
	for i in n:
		var x: float = plot_x + (float(i) + 0.5) / 90.0 * plot_w
		var y: float = plot_y + plot_h * (1.0 - _histogram[i] / max_pct)
		points[i] = Vector2(x, clampf(y, plot_y, plot_y + plot_h))

	if points.size() >= 2:
		draw_polyline(points, _line_color, 1.5, true)

	# Find and annotate peaks (local maxima well above noise floor)
	var sorted_vals: Array = []
	for v: float in _histogram:
		sorted_vals.append(v)
	sorted_vals.sort()
	var median: float = sorted_vals[sorted_vals.size() / 2] if sorted_vals.size() > 0 else 0.0
	var threshold: float = maxf(median * 2.0, max_pct * 0.15)

	var peak_bins: Array = []
	for i in range(2, n - 2):
		if _histogram[i] > threshold:
			# Must be >= both immediate and next neighbors (3-wide local max)
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
		var px: float = plot_x + (float(peak_bin) + 0.5) / 90.0 * plot_w
		var py: float = plot_y + plot_h * (1.0 - _histogram[peak_bin] / max_pct)
		draw_circle(Vector2(px, py), 3.0, _peak_color)
		var plbl: String = "%d" % peak_bin
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
			# Map mouse X to degree bin
			var deg_f: float = (mx - plot_x) / plot_w * 90.0
			var bin_idx: int = clampi(int(deg_f), 0, n - 1)
			var pct_val: float = _histogram[bin_idx]

			# Vertical crosshair line
			draw_line(Vector2(mx, plot_y), Vector2(mx, plot_y + plot_h),
				_crosshair_color, 1.0)

			# Horizontal line at the bin's value
			var val_y: float = plot_y + plot_h * (1.0 - pct_val / max_pct)
			val_y = clampf(val_y, plot_y, plot_y + plot_h)
			draw_line(Vector2(plot_x, val_y), Vector2(plot_x + plot_w, val_y),
				Color(_crosshair_color.r, _crosshair_color.g, _crosshair_color.b, 0.3), 1.0)

			# Intersection dot
			draw_circle(Vector2(mx, val_y), 4.0, _crosshair_color)

			# Tooltip label
			var tip: String = "%d°  %.2f%%" % [bin_idx, pct_val]
			var tip_x: float = mx + 8
			# Flip to left side if too close to right edge
			if tip_x + 70 > plot_x + plot_w:
				tip_x = mx - 75
			var tip_y: float = val_y - 4
			if tip_y < plot_y + 10:
				tip_y = val_y + 14
			draw_string(font, Vector2(tip_x, tip_y), tip,
				HORIZONTAL_ALIGNMENT_LEFT, -1, fs_sm + 1,
				Color(1.0, 1.0, 1.0, 0.95))
