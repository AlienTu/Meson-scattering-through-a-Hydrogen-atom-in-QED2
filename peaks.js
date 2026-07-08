function findLocalPeaks(xs, ys) {
  var ymax = 0;
  for (var i = 0; i < ys.length; i++) if (ys[i] > ymax) ymax = ys[i];
  var peaks = [];
  if (ymax <= 0) return peaks;
  for (var j = 1; j < ys.length - 1; j++) {
    if (xs[j] < 0.05 || xs[j] > 12) continue;
    if (ys[j] < 0.04 * ymax) continue;
    if (ys[j] > ys[j - 1] && ys[j] >= ys[j + 1]) {
      var omega = xs[j];
      var denom = ys[j - 1] - 2 * ys[j] + ys[j + 1];
      if (Math.abs(denom) > 1e-300) {
        omega = xs[j] + 0.5 * (ys[j - 1] - ys[j + 1]) / denom * (xs[j] - xs[j - 1]);
      }
      peaks.push({ omega: omega, rel: ys[j] / ymax, power: ys[j] });
    }
  }
  peaks.sort(function(a, b) { return b.power - a.power; });
  var kept = [];
  for (var k = 0; k < peaks.length; k++) {
    var ok = true;
    for (var m = 0; m < kept.length; m++) {
      if (Math.abs(peaks[k].omega - kept[m].omega) < 0.1) ok = false;
    }
    if (ok) kept.push(peaks[k]);
    if (kept.length >= 8) break;
  }
  kept.sort(function(a, b) { return a.omega - b.omega; });
  return kept;
}

function peakListText(name, xs, ys) {
  var ps = findLocalPeaks(xs, ys);
  if (ps.length === 0) return name + ': no clear peak';
  return name + ': ' + ps.map(function(p) {
    return p.omega.toFixed(4) + ' [' + (100 * p.rel).toFixed(1) + '%]';
  }).join(', ');
}
