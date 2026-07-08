function findLocalPeaks(xs, ys) {
  var ymax = 0;
  for (var i = 0; i < ys.length; i++) if (ys[i] > ymax) ymax = ys[i];
  var peaks = [];
  if (ymax <= 0) return peaks;
  for (var j = 1; j < ys.length - 1; j++) {
    if (xs[j] < 0.05 || xs[j] > 12) continue;
    if (ys[j] < 0.04 * ymax) continue;
    if (ys[j] > ys[j - 1] && ys[j] >= ys[j + 1]) peaks.push({ omega: xs[j], rel: ys[j] / ymax });
  }
  peaks.sort(function(a, b) { return b.rel - a.rel; });
  return peaks.slice(0, 8);
}
