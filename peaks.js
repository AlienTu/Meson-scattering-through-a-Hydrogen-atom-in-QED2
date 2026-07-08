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
      if (Math.abs(denom) > 1e-300) omega = xs[j] + 0.5 * (ys[j - 1] - ys[j + 1]) / denom * (xs[j] - xs[j - 1]);
      peaks.push({ omega: omega, rel: ys[j] / ymax, power: ys[j] });
    }
  }
  peaks.sort(function(a, b) { return b.power - a.power; });
  return peaks.slice(0, 8).sort(function(a, b) { return a.omega - b.omega; });
}

function peakListText(name, xs, ys) {
  var ps = findLocalPeaks(xs, ys);
  if (ps.length === 0) return name + ': no clear peak';
  return name + ': ' + ps.map(function(p) { return p.omega.toFixed(4); }).join(', ');
}

function spectrumPeakReport(spec) {
  if (!spec) return 'no spectrum';
  return [
    peakListText('left raw', spec.omega, spec.leftPower),
    peakListText('right raw', spec.omega, spec.rightPower),
    peakListText('incident at left', spec.omega, spec.rightMovingLeft),
    peakListText('reflected at left', spec.omega, spec.leftMovingLeft),
    peakListText('transmitted at right', spec.omega, spec.rightMovingRight)
  ].join('\n');
}
