function atomLinearOperator(vec, out) {
  var p = st.p;
  var n = p.N;
  var dx2 = p.dx * p.dx;
  for (var i = 0; i < n; i++) { out[i] = 0; out[n + i] = 0; }
  for (var j = 1; j < n - 1; j++) {
    var u1 = vec[j];
    var u2 = vec[n + j];
    var lap1 = (-vec[j - 1] + 2 * u1 - vec[j + 1]) / dx2;
    var lap2 = (-vec[n + j - 1] + 2 * u2 - vec[n + j + 1]) / dx2;
    var h11 = 0.5 * p.g * p.g + 2 * Math.PI * p.k1 * Math.cos(TWO_SQRT_PI * st.phi1Static[j]);
    var h22 = 0.5 * p.g * p.g + 2 * Math.PI * p.k2 * Math.cos(TWO_SQRT_PI * st.phi2Static[j]);
    var h12 = 0.5 * p.g * p.g;
    out[j] = lap1 + h11 * u1 + h12 * u2;
    out[n + j] = lap2 + h12 * u1 + h22 * u2;
  }
}

function atomDot(a, b) {
  var s = 0;
  var dx = st.p.dx;
  for (var i = 0; i < a.length; i++) s += a[i] * b[i] * dx;
  return s;
}

function atomNormalize(v) {
  var s = Math.sqrt(Math.max(atomDot(v, v), 1e-300));
  for (var i = 0; i < v.length; i++) v[i] /= s;
}

function atomProjectOut(v, z) {
  var c = atomDot(v, z);
  for (var i = 0; i < v.length; i++) v[i] -= c * z[i];
}

function translationMode() {
  var p = st.p;
  var n = p.N;
  var z = new Float64Array(2 * n);
  for (var i = 1; i < n - 1; i++) {
    z[i] = (st.phi1Static[i + 1] - st.phi1Static[i - 1]) / (2 * p.dx);
    z[n + i] = (st.phi2Static[i + 1] - st.phi2Static[i - 1]) / (2 * p.dx);
  }
  atomNormalize(z);
  return z;
}

function atomLocalization(v) {
  var p = st.p;
  var n = p.N;
  var near = 0;
  var total = 0;
  var R = Math.min(18, 0.25 * p.L);
  for (var i = 0; i < n; i++) {
    var den = (v[i] * v[i] + v[n + i] * v[n + i]) * p.dx;
    total += den;
    if (Math.abs(st.x[i]) < R) near += den;
  }
  return total > 0 ? near / total : 0;
}

function negCount2x2(a, b, d) {
  var tr = a + d;
  var disc = Math.sqrt(Math.max((a - d) * (a - d) + 4 * b * b, 0));
  var e1 = 0.5 * (tr - disc);
  var e2 = 0.5 * (tr + disc);
  return (e1 < 0 ? 1 : 0) + (e2 < 0 ? 1 : 0);
}

function blockSturmCount(lambda) {
  var p = st.p;
  var n = p.N;
  var m = n - 2;
  var aoff = 1 / (p.dx * p.dx);
  var prev00 = 0, prev01 = 0, prev11 = 0;
  var havePrev = false;
  var count = 0;
  for (var r = 0; r < m; r++) {
    var j = r + 1;
    var h11 = 0.5 * p.g * p.g + 2 * Math.PI * p.k1 * Math.cos(TWO_SQRT_PI * st.phi1Static[j]);
    var h22 = 0.5 * p.g * p.g + 2 * Math.PI * p.k2 * Math.cos(TWO_SQRT_PI * st.phi2Static[j]);
    var h12 = 0.5 * p.g * p.g;
    var s00 = 2 * aoff + h11 - lambda;
    var s01 = h12;
    var s11 = 2 * aoff + h22 - lambda;
    if (havePrev) {
      var det = prev00 * prev11 - prev01 * prev01;
      if (Math.abs(det) < 1e-18) det = det >= 0 ? 1e-18 : -1e-18;
      var inv00 = prev11 / det;
      var inv01 = -prev01 / det;
      var inv11 = prev00 / det;
      var aa = aoff * aoff;
      s00 -= aa * inv00;
      s01 -= aa * inv01;
      s11 -= aa * inv11;
    }
    count += negCount2x2(s00, s01, s11);
    prev00 = s00;
    prev01 = s01;
    prev11 = s11;
    havePrev = true;
  }
  return count;
}

function bisectEigenvalue(k, lo, hi) {
  for (var it = 0; it < 70; it++) {
    var mid = 0.5 * (lo + hi);
    if (blockSturmCount(mid) >= k) hi = mid;
    else lo = mid;
  }
  return 0.5 * (lo + hi);
}

function smallVectorOf2x2(a, b, c, d) {
  var aa = a * a + c * c;
  var bb = a * b + c * d;
  var dd = b * b + d * d;
  var tr = aa + dd;
  var disc = Math.sqrt(Math.max((aa - dd) * (aa - dd) + 4 * bb * bb, 0));
  var lam = 0.5 * (tr - disc);
  var x, y;
  if (Math.abs(bb) + Math.abs(aa - lam) > Math.abs(bb) + Math.abs(dd - lam)) {
    x = -bb;
    y = aa - lam;
  } else {
    x = dd - lam;
    y = -bb;
  }
  var norm = Math.sqrt(x * x + y * y) || 1;
  return [x / norm, y / norm];
}

function eigenvectorByShooting(lambda) {
  var p = st.p;
  var n = p.N;
  var m = n - 2;
  var aoff = 1 / (p.dx * p.dx);
  var yA0 = new Float64Array(m + 2);
  var yA1 = new Float64Array(m + 2);
  var yB0 = new Float64Array(m + 2);
  var yB1 = new Float64Array(m + 2);
  yA0[1] = 1;
  yA1[1] = 0;
  yB0[1] = 0;
  yB1[1] = 1;
  for (var r = 1; r <= m; r++) {
    var j = r;
    var site = j;
    var h11 = 0.5 * p.g * p.g + 2 * Math.PI * p.k1 * Math.cos(TWO_SQRT_PI * st.phi1Static[site]);
    var h22 = 0.5 * p.g * p.g + 2 * Math.PI * p.k2 * Math.cos(TWO_SQRT_PI * st.phi2Static[site]);
    var h12 = 0.5 * p.g * p.g;
    var d11 = 2 * aoff + h11 - lambda;
    var d22 = 2 * aoff + h22 - lambda;
    yA0[r + 1] = (d11 * yA0[r] + h12 * yA1[r]) / aoff - yA0[r - 1];
    yA1[r + 1] = (h12 * yA0[r] + d22 * yA1[r]) / aoff - yA1[r - 1];
    yB0[r + 1] = (d11 * yB0[r] + h12 * yB1[r]) / aoff - yB0[r - 1];
    yB1[r + 1] = (h12 * yB0[r] + d22 * yB1[r]) / aoff - yB1[r - 1];
    var scale = Math.max(Math.abs(yA0[r + 1]), Math.abs(yA1[r + 1]), Math.abs(yB0[r + 1]), Math.abs(yB1[r + 1]), 1);
    if (scale > 1e80) {
      for (var q = 0; q <= r + 1; q++) {
        yA0[q] /= scale;
        yA1[q] /= scale;
        yB0[q] /= scale;
        yB1[q] /= scale;
      }
    }
  }
  var c = smallVectorOf2x2(yA0[m + 1], yB0[m + 1], yA1[m + 1], yB1[m + 1]);
  var v = new Float64Array(2 * n);
  for (var rr = 1; rr <= m; rr++) {
    v[rr] = c[0] * yA0[rr] + c[1] * yB0[rr];
    v[n + rr] = c[0] * yA1[rr] + c[1] * yB1[rr];
  }
  atomNormalize(v);
  return v;
}

function solveInternalMode() {
  if (!st) buildAtom();
  var p = st.p;
  var n = p.N;
  var len = 2 * n;
  var zero = translationMode();
  var low = -10;
  while (blockSturmCount(low) > 0) low *= 2;
  var high = st.md.mLight2 * 0.999;
  var nBelow = blockSturmCount(high);
  var candidates = [];
  for (var k = 1; k <= Math.min(nBelow, 8); k++) {
    var lambda = bisectEigenvalue(k, low, high);
    var v = eigenvectorByShooting(lambda);
    var Lv = new Float64Array(len);
    atomLinearOperator(v, Lv);
    var ray = atomDot(v, Lv);
    var omega = Math.sqrt(Math.max(ray, 0));
    var loc = atomLocalization(v);
    var zeroOverlap = Math.abs(atomDot(v, zero));
    var residual = 0;
    for (var i = 0; i < len; i++) residual += (Lv[i] - ray * v[i]) * (Lv[i] - ray * v[i]) * p.dx;
    residual = Math.sqrt(residual);
    candidates.push({ k: k, lambda: ray, omega: omega, loc: loc, zero: zeroOverlap, residual: residual, u: v });
  }
  var chosen = null;
  var zeroCut = 0.25;
  for (var cidx = 0; cidx < candidates.length; cidx++) {
    var x = candidates[cidx];
    if (x.omega > zeroCut && x.omega < st.md.mLight && x.loc > 0.55 && x.zero < 0.50) { chosen = x; break; }
  }
  if (!chosen) {
    for (var cidx2 = 0; cidx2 < candidates.length; cidx2++) {
      var y = candidates[cidx2];
      if (y.omega > zeroCut && y.loc > 0.50) { chosen = y; break; }
    }
  }
  if (!chosen && candidates.length) chosen = candidates[candidates.length - 1];
  if (!chosen) {
    var box0 = document.getElementById('internalModeReadout');
    if (box0) box0.textContent = 'No sub-threshold eigenvalue found. Check static background.';
    return null;
  }
  st.internalMode = { u: chosen.u.slice(), omega: chosen.omega, lambda: chosen.lambda, candidates: candidates };
  var lines = [];
  lines.push('block-tridiagonal Sturm spectrum, no target frequency');
  lines.push('sub-threshold eigenvalue count = ' + nBelow);
  lines.push('chosen mode after zero-mode filtering: omega = ' + chosen.omega.toFixed(6) + ', omega^2 = ' + chosen.lambda.toFixed(6));
  lines.push('continuum estimate: m_light = ' + st.md.mLight.toFixed(6));
  lines.push('zero cutoff: omega > ' + zeroCut.toFixed(3));
  lines.push('localization = ' + chosen.loc.toFixed(4) + ', zero overlap = ' + chosen.zero.toExponential(2));
  lines.push('residual norm = ' + chosen.residual.toExponential(3));
  lines.push('sub-threshold modes:');
  for (var rline = 0; rline < candidates.length; rline++) {
    var cc = candidates[rline];
    lines.push(cc.k + ': omega=' + cc.omega.toFixed(6) + ', loc=' + cc.loc.toFixed(3) + ', zero=' + cc.zero.toExponential(1) + ', res=' + cc.residual.toExponential(1));
  }
  var box = document.getElementById('internalModeReadout');
  if (box) box.textContent = lines.join('\n');
  if (typeof message === 'function') message('bound mode solved without target shift: omega = ' + chosen.omega.toFixed(4));
  return st.internalMode;
}

function applyInternalMode() {
  if (!st || !st.internalMode) solveInternalMode();
  if (!st || !st.internalMode) return;
  var ampEl = document.getElementById('internalAmp');
  var phaseEl = document.getElementById('internalPhase');
  var amp = ampEl ? parseFloat(ampEl.value) : 0.08;
  var phase = phaseEl ? parseFloat(phaseEl.value) : 0;
  var p = st.p;
  var n = p.N;
  var u = st.internalMode.u;
  var omega = st.internalMode.omega;
  st.running = false;
  st.phi1 = st.phi1Static.slice();
  st.phi2 = st.phi2Static.slice();
  st.pi1 = new Float64Array(n);
  st.pi2 = new Float64Array(n);
  st.records = [];
  st.lastSpectrum = null;
  st.time = 0;
  st.step = 0;
  for (var i = 1; i < n - 1; i++) {
    st.phi1[i] += amp * Math.cos(phase) * u[i];
    st.phi2[i] += amp * Math.cos(phase) * u[n + i];
    st.pi1[i] += -amp * omega * Math.sin(phase) * u[i];
    st.pi2[i] += -amp * omega * Math.sin(phase) * u[n + i];
  }
  if (p.driveMode === 'packet') addIncomingPacket(st);
  st.initialEnergy = totalEnergy(st);
  drawAll();
  updateReadouts();
  if (typeof message === 'function') message('internal vibration applied');
}

function atomOnlyRun() {
  var dm = document.getElementById('driveMode');
  if (dm) dm.value = 'none';
  applyInternalMode();
  if (st) st.running = true;
}

function setupInternalButtons() {
  var solveBtn = document.getElementById('solveInternalBtn');
  var applyBtn = document.getElementById('applyInternalBtn');
  var atomBtn = document.getElementById('atomOnlyBtn');
  if (solveBtn) solveBtn.addEventListener('click', solveInternalMode);
  if (applyBtn) applyBtn.addEventListener('click', applyInternalMode);
  if (atomBtn) atomBtn.addEventListener('click', atomOnlyRun);
}

setupInternalButtons();
