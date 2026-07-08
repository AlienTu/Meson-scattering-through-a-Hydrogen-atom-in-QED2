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

function solveInternalMode() {
  if (!st) buildAtom();
  var p = st.p;
  var n = p.N;
  var len = 2 * n;
  var z = translationMode();
  var v = new Float64Array(len);
  var Lv = new Float64Array(len);
  var tmp = new Float64Array(len);
  var width = 8.0;
  for (var i = 1; i < n - 1; i++) {
    var env = Math.exp(-0.5 * st.x[i] * st.x[i] / (width * width));
    v[i] = env * Math.sin(0.7 * st.x[i]);
    v[n + i] = -env * Math.cos(0.5 * st.x[i]);
  }
  atomProjectOut(v, z);
  atomNormalize(v);
  var tau = 0.18 * p.dx * p.dx;
  var lambda = 0;
  for (var it = 0; it < 2600; it++) {
    atomLinearOperator(v, Lv);
    lambda = atomDot(v, Lv);
    for (var k = 0; k < len; k++) tmp[k] = v[k] - tau * (Lv[k] - lambda * v[k]);
    atomProjectOut(tmp, z);
    atomNormalize(tmp);
    var swap = v; v = tmp; tmp = swap;
  }
  atomLinearOperator(v, Lv);
  lambda = atomDot(v, Lv);
  var omega = Math.sqrt(Math.max(lambda, 0));
  st.internalMode = { u: v.slice(), omega: omega, lambda: lambda };
  var box = document.getElementById('internalModeReadout');
  if (box) box.textContent = 'omega_internal = ' + omega.toFixed(6) + ', omega^2 = ' + lambda.toFixed(6) + ', m_light = ' + st.md.mLight.toFixed(6);
  if (typeof message === 'function') message('internal mode solved: omega = ' + omega.toFixed(4));
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
