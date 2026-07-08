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

function jacobiEigenSym(A, maxIter) {
  var n = A.length;
  var V = [];
  for (var i = 0; i < n; i++) {
    V[i] = new Float64Array(n);
    V[i][i] = 1;
  }
  for (var it = 0; it < maxIter; it++) {
    var p = 0, q = 1, max = 0;
    for (var i2 = 0; i2 < n; i2++) {
      for (var j2 = i2 + 1; j2 < n; j2++) {
        var a = Math.abs(A[i2][j2]);
        if (a > max) { max = a; p = i2; q = j2; }
      }
    }
    if (max < 1e-9) break;
    var app = A[p][p], aqq = A[q][q], apq = A[p][q];
    var tau = (aqq - app) / (2 * apq);
    var t = Math.sign(tau || 1) / (Math.abs(tau) + Math.sqrt(1 + tau * tau));
    var c = 1 / Math.sqrt(1 + t * t);
    var s = t * c;
    for (var k = 0; k < n; k++) {
      if (k !== p && k !== q) {
        var akp = A[k][p], akq = A[k][q];
        A[k][p] = A[p][k] = c * akp - s * akq;
        A[k][q] = A[q][k] = s * akp + c * akq;
      }
    }
    A[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
    A[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
    A[p][q] = A[q][p] = 0;
    for (var k2 = 0; k2 < n; k2++) {
      var vkp = V[k2][p], vkq = V[k2][q];
      V[k2][p] = c * vkp - s * vkq;
      V[k2][q] = s * vkp + c * vkq;
    }
  }
  var vals = [];
  for (var i3 = 0; i3 < n; i3++) vals.push({ value: A[i3][i3], index: i3 });
  vals.sort(function(a, b) { return a.value - b.value; });
  return { values: vals, vectors: V };
}

function buildDirectFluctuationMatrix(windowHalfWidth) {
  var p = st.p;
  var dx2 = p.dx * p.dx;
  var idx = [];
  for (var i = 1; i < p.N - 1; i++) {
    if (Math.abs(st.x[i]) <= windowHalfWidth) idx.push(i);
  }
  var m = idx.length;
  var dim = 2 * m;
  var A = [];
  for (var r = 0; r < dim; r++) A[r] = new Float64Array(dim);
  for (var a = 0; a < m; a++) {
    var j = idx[a];
    var h11 = 0.5 * p.g * p.g + 2 * Math.PI * p.k1 * Math.cos(TWO_SQRT_PI * st.phi1Static[j]);
    var h22 = 0.5 * p.g * p.g + 2 * Math.PI * p.k2 * Math.cos(TWO_SQRT_PI * st.phi2Static[j]);
    var h12 = 0.5 * p.g * p.g;
    var r1 = a;
    var r2 = m + a;
    A[r1][r1] = 2 / dx2 + h11;
    A[r2][r2] = 2 / dx2 + h22;
    A[r1][r2] = h12;
    A[r2][r1] = h12;
    if (a > 0) {
      A[r1][r1 - 1] = -1 / dx2;
      A[r2][r2 - 1] = -1 / dx2;
    }
    if (a < m - 1) {
      A[r1][r1 + 1] = -1 / dx2;
      A[r2][r2 + 1] = -1 / dx2;
    }
  }
  return { A: A, idx: idx, m: m, dim: dim };
}

function fullVectorFromWindow(eigVec, idx, m) {
  var n = st.p.N;
  var v = new Float64Array(2 * n);
  for (var a = 0; a < m; a++) {
    var j = idx[a];
    v[j] = eigVec[a];
    v[n + j] = eigVec[m + a];
  }
  atomNormalize(v);
  return v;
}

function solveInternalMode() {
  if (!st) buildAtom();
  var p = st.p;
  var windowHalfWidth = Math.min(24, 0.35 * p.L);
  var pack = buildDirectFluctuationMatrix(windowHalfWidth);
  var eig = jacobiEigenSym(pack.A, 60 * pack.dim * pack.dim);
  var z = translationMode();
  var candidates = [];
  for (var c = 0; c < Math.min(30, eig.values.length); c++) {
    var col = eig.values[c].index;
    var lambda = eig.values[c].value;
    var ev = new Float64Array(pack.dim);
    for (var r = 0; r < pack.dim; r++) ev[r] = eig.vectors[r][col];
    var v = fullVectorFromWindow(ev, pack.idx, pack.m);
    var Lv = new Float64Array(2 * p.N);
    atomLinearOperator(v, Lv);
    lambda = atomDot(v, Lv);
    var omega = Math.sqrt(Math.max(lambda, 0));
    var loc = atomLocalization(v);
    var zero = Math.abs(atomDot(v, z));
    candidates.push({ lambda: lambda, omega: omega, loc: loc, zero: zero, u: v });
  }
  candidates.sort(function(a, b) { return a.lambda - b.lambda; });
  var chosen = null;
  for (var i = 0; i < candidates.length; i++) {
    var x = candidates[i];
    if (x.omega > 0.25 && x.omega < st.md.mLight && x.loc > 0.65) { chosen = x; break; }
  }
  if (!chosen) {
    for (var j2 = 0; j2 < candidates.length; j2++) {
      var y = candidates[j2];
      if (y.omega > 0.25 && y.loc > 0.55) { chosen = y; break; }
    }
  }
  if (!chosen) chosen = candidates[0];
  st.internalMode = { u: chosen.u.slice(), omega: chosen.omega, lambda: chosen.lambda, candidates: candidates };
  var lines = [];
  lines.push('direct matrix diagonalization on central window |x| < ' + windowHalfWidth.toFixed(1));
  lines.push('matrix dimension = ' + pack.dim + ' x ' + pack.dim);
  lines.push('chosen bound mode: omega = ' + chosen.omega.toFixed(6) + ', omega^2 = ' + chosen.lambda.toFixed(6));
  lines.push('continuum estimate: m_light = ' + st.md.mLight.toFixed(6));
  lines.push('expected benchmark: omega_B about 1.79, threshold about 2.92');
  lines.push('localization = ' + chosen.loc.toFixed(4) + ', zero overlap = ' + chosen.zero.toExponential(2));
  lines.push('lowest direct matrix modes:');
  for (var q = 0; q < Math.min(10, candidates.length); q++) {
    var cc = candidates[q];
    lines.push(q + ': omega=' + cc.omega.toFixed(6) + ', loc=' + cc.loc.toFixed(3) + ', zero=' + cc.zero.toExponential(1));
  }
  var box = document.getElementById('internalModeReadout');
  if (box) box.textContent = lines.join('\n');
  if (typeof message === 'function') message('direct matrix bound mode solved: omega = ' + chosen.omega.toFixed(4));
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
