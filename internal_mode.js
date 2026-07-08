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

function atomNorm(a) { return Math.sqrt(Math.max(atomDot(a, a), 1e-300)); }

function atomNormalize(v) {
  var s = atomNorm(v);
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

function jacobiEigenSym(A, maxIter) {
  var n = A.length;
  var V = [];
  for (var i = 0; i < n; i++) {
    V[i] = new Float64Array(n);
    V[i][i] = 1;
  }
  for (var it = 0; it < maxIter; it++) {
    var p = 0, q = 1, max = 0;
    for (var i = 0; i < n; i++) {
      for (var j = i + 1; j < n; j++) {
        var a = Math.abs(A[i][j]);
        if (a > max) { max = a; p = i; q = j; }
      }
    }
    if (max < 1e-10) break;
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
    for (var k = 0; k < n; k++) {
      var vkp = V[k][p], vkq = V[k][q];
      V[k][p] = c * vkp - s * vkq;
      V[k][q] = s * vkp + c * vkq;
    }
  }
  var vals = [];
  for (var i = 0; i < n; i++) vals.push({ value: A[i][i], index: i });
  vals.sort(function(a, b) { return a.value - b.value; });
  return { values: vals, vectors: V };
}

function solveInternalMode() {
  if (!st) buildAtom();
  var p = st.p;
  var n = p.N;
  var len = 2 * n;
  var z = translationMode();
  var m = Math.min(90, len - 2);
  var Q = [];
  var alpha = new Float64Array(m);
  var beta = new Float64Array(m);
  var q = new Float64Array(len);
  for (var i = 1; i < n - 1; i++) {
    var env = Math.exp(-0.5 * st.x[i] * st.x[i] / 64);
    q[i] = env * (0.7 + 0.3 * Math.sin(0.37 * st.x[i]));
    q[n + i] = -env * (0.9 + 0.2 * Math.cos(0.51 * st.x[i]));
  }
  atomProjectOut(q, z);
  atomNormalize(q);
  var qPrev = new Float64Array(len);
  var w = new Float64Array(len);
  var actualM = 0;
  for (var it = 0; it < m; it++) {
    Q.push(q.slice());
    atomLinearOperator(q, w);
    alpha[it] = atomDot(q, w);
    for (var k = 0; k < len; k++) w[k] -= alpha[it] * q[k] + (it > 0 ? beta[it - 1] * qPrev[k] : 0);
    atomProjectOut(w, z);
    for (var r = 0; r < Q.length; r++) atomProjectOut(w, Q[r]);
    beta[it] = atomNorm(w);
    actualM = it + 1;
    if (beta[it] < 1e-9 || it === m - 1) break;
    qPrev = q;
    q = w.slice();
    for (var k2 = 0; k2 < len; k2++) q[k2] /= beta[it];
  }
  var T = [];
  for (var a = 0; a < actualM; a++) {
    T[a] = new Float64Array(actualM);
    T[a][a] = alpha[a];
    if (a + 1 < actualM) T[a][a + 1] = T[a + 1][a] = beta[a];
  }
  var eig = jacobiEigenSym(T, 30 * actualM * actualM);
  var candidates = [];
  for (var c = 0; c < Math.min(20, eig.values.length); c++) {
    var ridx = eig.values[c].index;
    var lambda = eig.values[c].value;
    var v = new Float64Array(len);
    for (var a2 = 0; a2 < actualM; a2++) {
      var coeff = eig.vectors[a2][ridx];
      for (var b2 = 0; b2 < len; b2++) v[b2] += coeff * Q[a2][b2];
    }
    atomProjectOut(v, z);
    atomNormalize(v);
    var Lv = new Float64Array(len);
    atomLinearOperator(v, Lv);
    lambda = atomDot(v, Lv);
    var omega = Math.sqrt(Math.max(lambda, 0));
    var loc = atomLocalization(v);
    var ovZero = Math.abs(atomDot(v, z));
    candidates.push({ lambda: lambda, omega: omega, loc: loc, zero: ovZero, u: v });
  }
  candidates.sort(function(a, b) { return a.lambda - b.lambda; });
  var chosen = null;
  for (var c2 = 0; c2 < candidates.length; c2++) {
    var x = candidates[c2];
    if (x.lambda > 1e-6 && x.omega < st.md.mLight && x.loc > 0.55 && x.zero < 1e-3) { chosen = x; break; }
  }
  if (!chosen) {
    for (var c3 = 0; c3 < candidates.length; c3++) {
      var y = candidates[c3];
      if (y.lambda > 1e-6 && y.loc > 0.45 && y.zero < 1e-3) { chosen = y; break; }
    }
  }
  if (!chosen) chosen = candidates[0];
  st.internalMode = { u: chosen.u.slice(), omega: chosen.omega, lambda: chosen.lambda, candidates: candidates };
  var lines = [];
  lines.push('chosen bound mode: omega = ' + chosen.omega.toFixed(6) + ', omega^2 = ' + chosen.lambda.toFixed(6));
  lines.push('continuum estimate: m_light = ' + st.md.mLight.toFixed(6));
  lines.push('localization = ' + chosen.loc.toFixed(4) + ', zero overlap = ' + chosen.zero.toExponential(2));
  lines.push('lowest Ritz modes:');
  for (var c4 = 0; c4 < Math.min(8, candidates.length); c4++) {
    var cc = candidates[c4];
    lines.push(c4 + ': omega=' + cc.omega.toFixed(6) + ', loc=' + cc.loc.toFixed(3) + ', zero=' + cc.zero.toExponential(1));
  }
  var box = document.getElementById('internalModeReadout');
  if (box) box.textContent = lines.join('\n');
  if (typeof message === 'function') message('linear fluctuation mode solved: omega = ' + chosen.omega.toFixed(4));
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
