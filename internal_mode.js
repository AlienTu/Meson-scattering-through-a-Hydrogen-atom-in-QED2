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

function solve2x2(a, b, c, d, r0, r1) {
  var det = a * d - b * c;
  if (Math.abs(det) < 1e-14) det = det >= 0 ? 1e-14 : -1e-14;
  return [(d * r0 - b * r1) / det, (-c * r0 + a * r1) / det];
}

function solveShiftedBlockTridiagonal(rhs, sigma) {
  var p = st.p;
  var n = p.N;
  var m = n - 2;
  var aoff = 1 / (p.dx * p.dx);
  var P00 = new Float64Array(m);
  var P01 = new Float64Array(m);
  var P10 = new Float64Array(m);
  var P11 = new Float64Array(m);
  var G0 = new Float64Array(m);
  var G1 = new Float64Array(m);

  for (var r = 0; r < m; r++) {
    var j = r + 1;
    var h11 = 0.5 * p.g * p.g + 2 * Math.PI * p.k1 * Math.cos(TWO_SQRT_PI * st.phi1Static[j]);
    var h22 = 0.5 * p.g * p.g + 2 * Math.PI * p.k2 * Math.cos(TWO_SQRT_PI * st.phi2Static[j]);
    var h12 = 0.5 * p.g * p.g;

    var m00 = 2 * aoff + h11 - sigma;
    var m01 = h12;
    var m10 = h12;
    var m11 = 2 * aoff + h22 - sigma;

    var b0 = rhs[j];
    var b1 = rhs[n + j];
    if (r > 0) {
      m00 += aoff * P00[r - 1];
      m01 += aoff * P01[r - 1];
      m10 += aoff * P10[r - 1];
      m11 += aoff * P11[r - 1];
      b0 += aoff * G0[r - 1];
      b1 += aoff * G1[r - 1];
    }

    var g = solve2x2(m00, m01, m10, m11, b0, b1);
    G0[r] = g[0];
    G1[r] = g[1];

    if (r < m - 1) {
      var pcol0 = solve2x2(m00, m01, m10, m11, -aoff, 0);
      var pcol1 = solve2x2(m00, m01, m10, m11, 0, -aoff);
      P00[r] = pcol0[0];
      P10[r] = pcol0[1];
      P01[r] = pcol1[0];
      P11[r] = pcol1[1];
    }
  }

  var y = new Float64Array(2 * n);
  for (var r2 = m - 1; r2 >= 0; r2--) {
    var y0 = G0[r2];
    var y1 = G1[r2];
    if (r2 < m - 1) {
      var yn0 = y[r2 + 2];
      var yn1 = y[n + r2 + 2];
      y0 -= P00[r2] * yn0 + P01[r2] * yn1;
      y1 -= P10[r2] * yn0 + P11[r2] * yn1;
    }
    y[r2 + 1] = y0;
    y[n + r2 + 1] = y1;
  }
  return y;
}

function solveInternalMode() {
  if (!st) buildAtom();
  var p = st.p;
  var n = p.N;
  var len = 2 * n;
  var z = translationMode();
  var targetOmega = Math.min(1.79, 0.75 * st.md.mLight);
  var sigma = targetOmega * targetOmega;
  var q = new Float64Array(len);
  var width = 7.0;
  for (var i = 1; i < n - 1; i++) {
    var env = Math.exp(-0.5 * st.x[i] * st.x[i] / (width * width));
    q[i] = env;
    q[n + i] = -env;
  }
  atomProjectOut(q, z);
  atomNormalize(q);

  var iterations = 28;
  var lambda = 0;
  var Lv = new Float64Array(len);
  for (var it = 0; it < iterations; it++) {
    var y = solveShiftedBlockTridiagonal(q, sigma);
    atomProjectOut(y, z);
    atomNormalize(y);
    q = y;
    atomLinearOperator(q, Lv);
    lambda = atomDot(q, Lv);
  }
  atomLinearOperator(q, Lv);
  lambda = atomDot(q, Lv);
  var omega = Math.sqrt(Math.max(lambda, 0));
  var loc = atomLocalization(q);
  var zero = Math.abs(atomDot(q, z));
  var residual = 0;
  for (var r = 0; r < len; r++) residual += (Lv[r] - lambda * q[r]) * (Lv[r] - lambda * q[r]) * p.dx;
  residual = Math.sqrt(residual);

  st.internalMode = { u: q.slice(), omega: omega, lambda: lambda, sigma: sigma, residual: residual };
  var lines = [];
  lines.push('shift-invert block-tridiagonal solve');
  lines.push('target omega = ' + targetOmega.toFixed(6) + ', sigma = ' + sigma.toFixed(6));
  lines.push('omega = ' + omega.toFixed(6) + ', omega^2 = ' + lambda.toFixed(6));
  lines.push('continuum estimate: m_light = ' + st.md.mLight.toFixed(6));
  lines.push('benchmark: omega_B about 1.79, threshold about 2.92');
  lines.push('localization = ' + loc.toFixed(4) + ', zero overlap = ' + zero.toExponential(2));
  lines.push('residual norm = ' + residual.toExponential(3));
  lines.push('cost: O(N) per inverse iteration, no dense diagonalization');
  var box = document.getElementById('internalModeReadout');
  if (box) box.textContent = lines.join('\n');
  if (typeof message === 'function') message('fast shift-invert bound mode solved: omega = ' + omega.toFixed(4));
  return st.internalMode;
}

function applyInternalMode() {
  if (!st || !st.internalMode) solveInternalMode();
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
  st.running = true;
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
