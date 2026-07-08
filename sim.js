"use strict";

const SQRT_PI = Math.sqrt(Math.PI);
const TWO_SQRT_PI = 2 * SQRT_PI;

const $ = (id) => document.getElementById(id);

const ui = {
  buildBtn: $("buildBtn"), startBtn: $("startBtn"), pauseBtn: $("pauseBtn"), resetBtn: $("resetBtn"),
  spectrumBtn: $("spectrumBtn"), downloadBtn: $("downloadBtn"), clearDataBtn: $("clearDataBtn"),
  message: $("message"), time: $("time"), Nreadout: $("Nreadout"), dxReadout: $("dxReadout"),
  massReadout: $("massReadout"), kReadout: $("kReadout"), energyReadout: $("energyReadout"),
  fieldCanvas: $("fieldCanvas"), leftTraceCanvas: $("leftTraceCanvas"), rightTraceCanvas: $("rightTraceCanvas"),
  spectrumCanvas: $("spectrumCanvas"), directionCanvas: $("directionCanvas")
};

let st = null;

function num(id) { return parseFloat($(id).value); }
function intNum(id) { return parseInt($(id).value, 10); }
function clamp(v, lo, hi) { return Math.max(lo, Math.min(hi, v)); }
function sqr(v) { return v * v; }

function message(text) { ui.message.textContent = text; }

function readParams() {
  const L = num("L");
  const requestedDx = num("dx");
  const N = Math.round(L / requestedDx) + 1;
  const dx = L / (N - 1);
  const p = {
    L, N, dx,
    g: num("g"), k1: num("k1"), k2: num("k2"),
    dt: num("dt"), substeps: intNum("substeps"), relaxSteps: intNum("relaxSteps"),
    driveMode: $("driveMode").value,
    omegaIn: num("omegaIn"), ampIn: num("ampIn"), packetWidth: num("packetWidth"),
    sourceX: num("sourceX"), sourceWidth: num("sourceWidth"),
    spongeStart: num("spongeStart"), spongeStrength: num("spongeStrength"),
    xLeftProbe: num("xLeftProbe"), xRightProbe: num("xRightProbe"), probeSep: num("probeSep"),
    recordEvery: intNum("recordEvery"), maxRecords: intNum("maxRecords")
  };
  p.dt = Math.min(p.dt, 0.45 * p.dx);
  return p;
}

function makeGrid(p) {
  const x = new Float64Array(p.N);
  for (let i = 0; i < p.N; i++) x[i] = -0.5 * p.L + i * p.dx;
  return x;
}

function massData(p) {
  const a = 0.5 * p.g * p.g + 2 * Math.PI * p.k1;
  const d = 0.5 * p.g * p.g + 2 * Math.PI * p.k2;
  const b = 0.5 * p.g * p.g;
  const tr = a + d;
  const disc = Math.sqrt((a - d) * (a - d) + 4 * b * b);
  const lamLight = 0.5 * (tr - disc);
  const lamHeavy = 0.5 * (tr + disc);
  let v1, v2;
  if (Math.abs(b) > 1e-14) {
    v1 = b;
    v2 = lamLight - a;
  } else {
    v1 = 0; v2 = 1;
  }
  const norm = Math.sqrt(v1 * v1 + v2 * v2) || 1;
  v1 /= norm; v2 /= norm;
  if (v2 < 0) { v1 = -v1; v2 = -v2; }
  return {
    mLight2: Math.max(lamLight, 0), mHeavy2: Math.max(lamHeavy, 0),
    mLight: Math.sqrt(Math.max(lamLight, 0)), mHeavy: Math.sqrt(Math.max(lamHeavy, 0)),
    vLight: [v1, v2]
  };
}

function indexAt(x, p) {
  return clamp(Math.round((x + 0.5 * p.L) / p.dx), 0, p.N - 1);
}

function potentialGradient(phi1, phi2, p) {
  const common = 0.5 * p.g * p.g * (phi1 + phi2);
  return [
    common + SQRT_PI * p.k1 * Math.sin(TWO_SQRT_PI * phi1),
    common + SQRT_PI * p.k2 * Math.sin(TWO_SQRT_PI * phi2)
  ];
}

function buildAtom() {
  const p = readParams();
  const x = makeGrid(p);
  const phi1 = new Float64Array(p.N);
  const phi2 = new Float64Array(p.N);
  const v1 = new Float64Array(p.N);
  const v2 = new Float64Array(p.N);
  const left1 = 0.0, left2 = 0.0;
  const right1 = SQRT_PI, right2 = -SQRT_PI;
  const width = 3.0;

  for (let i = 0; i < p.N; i++) {
    const s = 0.5 * (1 + Math.tanh(x[i] / width));
    phi1[i] = left1 + (right1 - left1) * s;
    phi2[i] = left2 + (right2 - left2) * s;
  }
  phi1[0] = left1; phi2[0] = left2; phi1[p.N - 1] = right1; phi2[p.N - 1] = right2;

  const relaxDt = Math.min(0.035, 0.25 * p.dx);
  const damp = 1.8;
  let maxResidual = Infinity;
  for (let it = 0; it < p.relaxSteps; it++) {
    maxResidual = 0;
    for (let i = 1; i < p.N - 1; i++) {
      const lap1 = (phi1[i - 1] - 2 * phi1[i] + phi1[i + 1]) / (p.dx * p.dx);
      const lap2 = (phi2[i - 1] - 2 * phi2[i] + phi2[i + 1]) / (p.dx * p.dx);
      const [d1, d2] = potentialGradient(phi1[i], phi2[i], p);
      const f1 = lap1 - d1;
      const f2 = lap2 - d2;
      v1[i] += relaxDt * (f1 - damp * v1[i]);
      v2[i] += relaxDt * (f2 - damp * v2[i]);
      maxResidual = Math.max(maxResidual, Math.abs(f1), Math.abs(f2));
    }
    for (let i = 1; i < p.N - 1; i++) {
      phi1[i] += relaxDt * v1[i];
      phi2[i] += relaxDt * v2[i];
    }
    phi1[0] = left1; phi2[0] = left2; phi1[p.N - 1] = right1; phi2[p.N - 1] = right2;
    v1[0] = 0; v2[0] = 0; v1[p.N - 1] = 0; v2[p.N - 1] = 0;
  }

  const md = massData(p);
  st = {
    p, x, md,
    phi1Static: phi1.slice(), phi2Static: phi2.slice(),
    phi1: phi1.slice(), phi2: phi2.slice(),
    pi1: new Float64Array(p.N), pi2: new Float64Array(p.N),
    force1: new Float64Array(p.N), force2: new Float64Array(p.N),
    gamma: makeSponge(p, x),
    time: 0, step: 0, running: false, built: true,
    records: [], initialEnergy: null, lastSpectrum: null,
    probe: updateProbeIndices(p)
  };
  if (p.driveMode === "packet") addIncomingPacket(st);
  st.initialEnergy = totalEnergy(st);
  updateReadouts(maxResidual);
  drawAll();
  message(`atom built; relaxation residual ${maxResidual.toExponential(2)}`);
}

function makeSponge(p, x) {
  const gamma = new Float64Array(p.N);
  const denom = Math.max(0.5 * p.L - p.spongeStart, p.dx);
  for (let i = 0; i < p.N; i++) {
    const ax = Math.abs(x[i]);
    if (ax > p.spongeStart) {
      const s = (ax - p.spongeStart) / denom;
      gamma[i] = p.spongeStrength * s * s;
    }
  }
  return gamma;
}

function updateProbeIndices(p) {
  return {
    L1: indexAt(p.xLeftProbe, p),
    L2: indexAt(p.xLeftProbe + p.probeSep, p),
    R1: indexAt(p.xRightProbe, p),
    R2: indexAt(p.xRightProbe + p.probeSep, p)
  };
}

function addIncomingPacket(s) {
  const p = s.p;
  const omega = p.omegaIn;
  const k = Math.sqrt(Math.max(omega * omega - s.md.mLight2, 0));
  const x0 = p.sourceX;
  const sig = p.packetWidth;
  const [v1, v2] = s.md.vLight;
  if (k <= 0) {
    message("omega_in is below the light meson mass; packet is evanescent.");
    return;
  }
  for (let i = 1; i < p.N - 1; i++) {
    const y = s.x[i] - x0;
    const env = Math.exp(-0.5 * y * y / (sig * sig));
    const phase = k * y;
    const val = p.ampIn * env * Math.cos(phase);
    const vel = p.ampIn * omega * env * Math.sin(phase);
    s.phi1[i] += v1 * val;
    s.phi2[i] += v2 * val;
    s.pi1[i] += v1 * vel;
    s.pi2[i] += v2 * vel;
  }
}

function computeForces(s) {
  const p = s.p;
  const f1 = s.force1, f2 = s.force2;
  const [mv1, mv2] = s.md.vLight;
  const continuous = p.driveMode === "continuous";
  const on = continuous ? (1 - Math.exp(-s.time * s.time / 25.0)) : 0;
  const drivePhase = Math.sin(p.omegaIn * s.time);
  for (let i = 1; i < p.N - 1; i++) {
    const lap1 = (s.phi1[i - 1] - 2 * s.phi1[i] + s.phi1[i + 1]) / (p.dx * p.dx);
    const lap2 = (s.phi2[i - 1] - 2 * s.phi2[i] + s.phi2[i + 1]) / (p.dx * p.dx);
    const [d1, d2] = potentialGradient(s.phi1[i], s.phi2[i], p);
    let src1 = 0, src2 = 0;
    if (continuous) {
      const y = s.x[i] - p.sourceX;
      const mask = Math.exp(-0.5 * y * y / (p.sourceWidth * p.sourceWidth));
      const amp = p.ampIn * on * drivePhase * mask;
      src1 = mv1 * amp;
      src2 = mv2 * amp;
    }
    f1[i] = lap1 - d1 - s.gamma[i] * s.pi1[i] + src1;
    f2[i] = lap2 - d2 - s.gamma[i] * s.pi2[i] + src2;
  }
  f1[0] = 0; f2[0] = 0; f1[p.N - 1] = 0; f2[p.N - 1] = 0;
}

function stepSimulation(s) {
  const p = s.p;
  computeForces(s);
  for (let i = 1; i < p.N - 1; i++) {
    s.pi1[i] += p.dt * s.force1[i];
    s.pi2[i] += p.dt * s.force2[i];
    s.phi1[i] += p.dt * s.pi1[i];
    s.phi2[i] += p.dt * s.pi2[i];
  }
  s.phi1[0] = s.phi1Static[0]; s.phi2[0] = s.phi2Static[0];
  s.phi1[p.N - 1] = s.phi1Static[p.N - 1]; s.phi2[p.N - 1] = s.phi2Static[p.N - 1];
  s.pi1[0] = 0; s.pi2[0] = 0; s.pi1[p.N - 1] = 0; s.pi2[p.N - 1] = 0;
  s.time += p.dt;
  s.step++;
  if (s.step % p.recordEvery === 0) recordProbes(s);
}

function recordProbes(s) {
  const pr = s.probe;
  const rec = { t: s.time };
  for (const [name, idx] of Object.entries(pr)) {
    rec[`x_${name}`] = s.x[idx];
    rec[`phi1_${name}`] = s.phi1[idx];
    rec[`phi2_${name}`] = s.phi2[idx];
    rec[`dphi1_${name}`] = s.phi1[idx] - s.phi1Static[idx];
    rec[`dphi2_${name}`] = s.phi2[idx] - s.phi2Static[idx];
    rec[`pi1_${name}`] = s.pi1[idx];
    rec[`pi2_${name}`] = s.pi2[idx];
  }
  s.records.push(rec);
  if (s.records.length > s.p.maxRecords) s.records.shift();
}

function totalEnergy(s) {
  const p = s.p;
  let E = 0;
  for (let i = 1; i < p.N - 1; i++) {
    const d1 = (s.phi1[i + 1] - s.phi1[i - 1]) / (2 * p.dx);
    const d2 = (s.phi2[i + 1] - s.phi2[i - 1]) / (2 * p.dx);
    const kin = 0.5 * (s.pi1[i] * s.pi1[i] + s.pi2[i] * s.pi2[i]);
    const grad = 0.5 * (d1 * d1 + d2 * d2);
    const pot = 0.25 * p.g * p.g * sqr(s.phi1[i] + s.phi2[i])
      - 0.5 * p.k1 * Math.cos(TWO_SQRT_PI * s.phi1[i])
      - 0.5 * p.k2 * Math.cos(TWO_SQRT_PI * s.phi2[i]);
    E += (kin + grad + pot) * p.dx;
  }
  return E;
}

function updateReadouts(residual = null) {
  if (!st) return;
  const p = st.p;
  const omega = p.omegaIn;
  const k = Math.sqrt(Math.max(omega * omega - st.md.mLight2, 0));
  ui.time.textContent = st.time.toFixed(3);
  ui.Nreadout.textContent = p.N.toString();
  ui.dxReadout.textContent = p.dx.toFixed(4);
  ui.massReadout.textContent = `m_light=${st.md.mLight.toFixed(3)}, m_heavy=${st.md.mHeavy.toFixed(3)}`;
  ui.kReadout.textContent = k > 0 ? k.toFixed(3) : "evanescent";
  const E = totalEnergy(st);
  const drift = st.initialEnergy === null ? 0 : E - st.initialEnergy;
  ui.energyReadout.textContent = `${drift.toExponential(3)}${residual !== null ? `, residual=${residual.toExponential(2)}` : ""}`;
}

function drawAll() {
  if (!st) return;
  drawFields();
  drawTrace(ui.leftTraceCanvas, "L");
  drawTrace(ui.rightTraceCanvas, "R");
  if (st.lastSpectrum) {
    drawSpectrum(st.lastSpectrum);
    drawDirection(st.lastSpectrum);
  }
}

function panel(ctx, x, y, w, h, title) {
  ctx.strokeStyle = "#30363d";
  ctx.strokeRect(x, y, w, h);
  ctx.fillStyle = "#c9d1d9";
  ctx.font = "14px sans-serif";
  ctx.fillText(title, x + 8, y + 18);
}

function drawCurve(ctx, arr, x0, y0, w, h, ymin, ymax, color, lineWidth = 1.5, dashed = false) {
  const N = arr.length;
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = lineWidth;
  if (dashed) ctx.setLineDash([5, 4]);
  ctx.beginPath();
  for (let i = 0; i < N; i++) {
    const px = x0 + w * i / (N - 1);
    const py = y0 + h * (1 - (arr[i] - ymin) / (ymax - ymin || 1));
    if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
  }
  ctx.stroke();
  ctx.restore();
}

function drawFields() {
  const c = ui.fieldCanvas;
  const ctx = c.getContext("2d");
  const W = c.width, H = c.height;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = "#080b10";
  ctx.fillRect(0, 0, W, H);
  const pad = 42;
  const h = (H - 3 * pad) / 2;
  const w = W - 2 * pad;
  panel(ctx, pad, pad, w, h, "full fields: phi1, phi2, static background dashed");
  panel(ctx, pad, 2 * pad + h, w, h, "vibrations: phi - phi_static");

  let ymin = Infinity, ymax = -Infinity, dymin = Infinity, dymax = -Infinity;
  const d1 = new Float64Array(st.p.N), d2 = new Float64Array(st.p.N);
  for (let i = 0; i < st.p.N; i++) {
    ymin = Math.min(ymin, st.phi1[i], st.phi2[i], st.phi1Static[i], st.phi2Static[i]);
    ymax = Math.max(ymax, st.phi1[i], st.phi2[i], st.phi1Static[i], st.phi2Static[i]);
    d1[i] = st.phi1[i] - st.phi1Static[i];
    d2[i] = st.phi2[i] - st.phi2Static[i];
    dymin = Math.min(dymin, d1[i], d2[i]);
    dymax = Math.max(dymax, d1[i], d2[i]);
  }
  const ypad = 0.08 * (ymax - ymin || 1); ymin -= ypad; ymax += ypad;
  const dabs = Math.max(Math.abs(dymin), Math.abs(dymax), 1e-4); dymin = -dabs; dymax = dabs;
  drawCurve(ctx, st.phi1Static, pad, pad, w, h, ymin, ymax, "#8b949e", 1, true);
  drawCurve(ctx, st.phi2Static, pad, pad, w, h, ymin, ymax, "#8b949e", 1, true);
  drawCurve(ctx, st.phi1, pad, pad, w, h, ymin, ymax, "#58a6ff", 1.8);
  drawCurve(ctx, st.phi2, pad, pad, w, h, ymin, ymax, "#f2cc60", 1.8);
  drawCurve(ctx, d1, pad, 2 * pad + h, w, h, dymin, dymax, "#58a6ff", 1.6);
  drawCurve(ctx, d2, pad, 2 * pad + h, w, h, dymin, dymax, "#f2cc60", 1.6);
  markX(ctx, st.p.xLeftProbe, pad, pad, w, 2 * h + pad, "#3fb950", "L");
  markX(ctx, st.p.xRightProbe, pad, pad, w, 2 * h + pad, "#f85149", "R");
  ctx.fillStyle = "#58a6ff"; ctx.fillText("phi1", pad + w - 90, pad + 20);
  ctx.fillStyle = "#f2cc60"; ctx.fillText("phi2", pad + w - 45, pad + 20);
}

function markX(ctx, xval, x0, y0, w, h, color, label) {
  const px = x0 + w * (xval + 0.5 * st.p.L) / st.p.L;
  ctx.save();
  ctx.strokeStyle = color;
  ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(px, y0); ctx.lineTo(px, y0 + h); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = color; ctx.fillText(label, px + 4, y0 + 14);
  ctx.restore();
}

function drawTrace(canvas, side) {
  const ctx = canvas.getContext("2d");
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = "#080b10"; ctx.fillRect(0, 0, W, H);
  const recs = st.records.slice(-1200);
  if (recs.length < 2) return;
  const key1 = `dphi1_${side}1`, key2 = `dphi2_${side}1`;
  let ymax = 1e-5;
  for (const r of recs) ymax = Math.max(ymax, Math.abs(r[key1]), Math.abs(r[key2]));
  const arr1 = recs.map(r => r[key1]);
  const arr2 = recs.map(r => r[key2]);
  drawAxes(ctx, W, H, -ymax, ymax);
  drawCurve(ctx, arr1, 38, 18, W - 58, H - 48, -ymax, ymax, "#58a6ff", 1.5);
  drawCurve(ctx, arr2, 38, 18, W - 58, H - 48, -ymax, ymax, "#f2cc60", 1.5);
  ctx.fillStyle = "#8b949e"; ctx.font = "12px sans-serif";
  ctx.fillText(`last ${recs.length} records, t=${st.time.toFixed(2)}`, 48, H - 12);
}

function drawAxes(ctx, W, H, ymin, ymax) {
  ctx.strokeStyle = "#30363d";
  ctx.strokeRect(38, 18, W - 58, H - 48);
  const y0 = 18 + (H - 48) * (1 - (0 - ymin) / (ymax - ymin || 1));
  ctx.beginPath(); ctx.moveTo(38, y0); ctx.lineTo(W - 20, y0); ctx.stroke();
}

function lightProjectionRecord(r, side, second = false) {
  const tag = second ? `${side}2` : `${side}1`;
  const [v1, v2] = st.md.vLight;
  return v1 * r[`dphi1_${tag}`] + v2 * r[`dphi2_${tag}`];
}

function hann(n, N) { return 0.5 * (1 - Math.cos(2 * Math.PI * n / Math.max(1, N - 1))); }

function fft(re, im) {
  const n = re.length;
  let j = 0;
  for (let i = 1; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = -2 * Math.PI / len;
    const wlenRe = Math.cos(ang), wlenIm = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let wRe = 1, wIm = 0;
      for (let k = 0; k < len / 2; k++) {
        const uRe = re[i + k], uIm = im[i + k];
        const vRe = re[i + k + len / 2] * wRe - im[i + k + len / 2] * wIm;
        const vIm = re[i + k + len / 2] * wIm + im[i + k + len / 2] * wRe;
        re[i + k] = uRe + vRe;
        im[i + k] = uIm + vIm;
        re[i + k + len / 2] = uRe - vRe;
        im[i + k + len / 2] = uIm - vIm;
        const nwRe = wRe * wlenRe - wIm * wlenIm;
        const nwIm = wRe * wlenIm + wIm * wlenRe;
        wRe = nwRe; wIm = nwIm;
      }
    }
  }
}

function pow2Floor(n) {
  let p = 1;
  while ((p << 1) <= n) p <<= 1;
  return p;
}

function complexExp(theta) { return { re: Math.cos(theta), im: Math.sin(theta) }; }
function cMul(a, b) { return { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }; }
function cSub(a, b) { return { re: a.re - b.re, im: a.im - b.im }; }
function cDiv(a, b) { const den = b.re * b.re + b.im * b.im || 1e-300; return { re: (a.re * b.re + a.im * b.im) / den, im: (a.im * b.re - a.re * b.im) / den }; }
function cAbs2(a) { return a.re * a.re + a.im * a.im; }

function analyzeSpectrum() {
  if (!st || st.records.length < 32) { message("not enough probe records for spectrum"); return; }
  const n = pow2Floor(st.records.length);
  const recs = st.records.slice(st.records.length - n);
  const dtS = st.p.dt * st.p.recordEvery;
  const signals = {};
  for (const key of ["L1", "L2", "R1", "R2"]) {
    const re = new Float64Array(n), im = new Float64Array(n);
    const side = key[0], second = key[1] === "2";
    let mean = 0;
    for (let i = 0; i < n; i++) mean += lightProjectionRecord(recs[i], side, second);
    mean /= n;
    for (let i = 0; i < n; i++) re[i] = (lightProjectionRecord(recs[i], side, second) - mean) * hann(i, n);
    fft(re, im);
    signals[key] = { re, im };
  }
  const omega = [], leftPower = [], rightPower = [];
  const leftMovingLeft = [], rightMovingRight = [], rightMovingLeft = [], leftMovingRight = [];
  const xL1 = recs[n - 1].x_L1, xL2 = recs[n - 1].x_L2;
  const xR1 = recs[n - 1].x_R1, xR2 = recs[n - 1].x_R2;

  for (let m = 1; m < n / 2; m++) {
    const om = 2 * Math.PI * m / (n * dtS);
    omega.push(om);
    const CL = { re: signals.L1.re[m] / n, im: signals.L1.im[m] / n };
    const CR = { re: signals.R1.re[m] / n, im: signals.R1.im[m] / n };
    leftPower.push(cAbs2(CL)); rightPower.push(cAbs2(CR));

    const dirL = directionalAtFrequency(om, xL1, xL2,
      { re: signals.L1.re[m] / n, im: signals.L1.im[m] / n },
      { re: signals.L2.re[m] / n, im: signals.L2.im[m] / n });
    const dirR = directionalAtFrequency(om, xR1, xR2,
      { re: signals.R1.re[m] / n, im: signals.R1.im[m] / n },
      { re: signals.R2.re[m] / n, im: signals.R2.im[m] / n });
    rightMovingLeft.push(dirL.right);
    leftMovingLeft.push(dirL.left);
    rightMovingRight.push(dirR.right);
    leftMovingRight.push(dirR.left);
  }
  st.lastSpectrum = { omega, leftPower, rightPower, leftMovingLeft, rightMovingLeft, leftMovingRight, rightMovingRight };
  drawSpectrum(st.lastSpectrum);
  drawDirection(st.lastSpectrum);
  message(`spectrum computed with ${n} samples`);
}

function directionalAtFrequency(omega, x1, x2, C1, C2) {
  if (omega <= st.md.mLight) return { right: 0, left: 0 };
  const k = Math.sqrt(Math.max(omega * omega - st.md.mLight2, 0));
  if (k * Math.abs(x2 - x1) < 1e-6) return { right: 0, left: 0 };
  const a1 = complexExp(k * x1), b1 = complexExp(-k * x1);
  const a2 = complexExp(k * x2), b2 = complexExp(-k * x2);
  const den = cSub(cMul(a1, b2), cMul(a2, b1));
  const R = cDiv(cSub(cMul(C1, b2), cMul(C2, b1)), den);
  const L = cDiv(cSub(cMul(a1, C2), cMul(a2, C1)), den);
  return { right: cAbs2(R), left: cAbs2(L) };
}

function drawSpectrum(spec) {
  const ctx = ui.spectrumCanvas.getContext("2d");
  const W = ui.spectrumCanvas.width, H = ui.spectrumCanvas.height;
  ctx.clearRect(0, 0, W, H); ctx.fillStyle = "#080b10"; ctx.fillRect(0, 0, W, H);
  if (!spec.omega.length) return;
  const omMax = Math.min(12, spec.omega[spec.omega.length - 1]);
  const vals = [];
  for (let i = 0; i < spec.omega.length; i++) if (spec.omega[i] <= omMax) vals.push(spec.leftPower[i], spec.rightPower[i]);
  const ymax = Math.max(...vals, 1e-16);
  drawSpectralAxes(ctx, W, H, omMax, ymax, "raw light-mode power: left and right probes");
  drawXY(ctx, spec.omega, spec.leftPower, omMax, ymax, "#3fb950");
  drawXY(ctx, spec.omega, spec.rightPower, omMax, ymax, "#f85149");
  ctx.fillStyle = "#3fb950"; ctx.fillText("left probe", 52, 42);
  ctx.fillStyle = "#f85149"; ctx.fillText("right probe", 140, 42);
}

function drawDirection(spec) {
  const ctx = ui.directionCanvas.getContext("2d");
  const W = ui.directionCanvas.width, H = ui.directionCanvas.height;
  ctx.clearRect(0, 0, W, H); ctx.fillStyle = "#080b10"; ctx.fillRect(0, 0, W, H);
  if (!spec.omega.length) return;
  const omMax = Math.min(12, spec.omega[spec.omega.length - 1]);
  const vals = [];
  for (let i = 0; i < spec.omega.length; i++) if (spec.omega[i] <= omMax) vals.push(spec.leftMovingLeft[i], spec.rightMovingRight[i], spec.rightMovingLeft[i]);
  const ymax = Math.max(...vals, 1e-16);
  drawSpectralAxes(ctx, W, H, omMax, ymax, "directional power, two-probe estimate");
  drawXY(ctx, spec.omega, spec.leftMovingLeft, omMax, ymax, "#f2cc60");
  drawXY(ctx, spec.omega, spec.rightMovingRight, omMax, ymax, "#58a6ff");
  drawXY(ctx, spec.omega, spec.rightMovingLeft, omMax, ymax, "#8b949e");
  ctx.fillStyle = "#f2cc60"; ctx.fillText("reflected at left", 52, 42);
  ctx.fillStyle = "#58a6ff"; ctx.fillText("transmitted at right", 170, 42);
  ctx.fillStyle = "#8b949e"; ctx.fillText("incident at left", 330, 42);
}

function drawSpectralAxes(ctx, W, H, omMax, ymax, title) {
  ctx.strokeStyle = "#30363d"; ctx.strokeRect(42, 24, W - 64, H - 64);
  ctx.fillStyle = "#c9d1d9"; ctx.font = "12px sans-serif"; ctx.fillText(title, 52, 18);
  ctx.fillStyle = "#8b949e";
  ctx.fillText("omega", W - 58, H - 14);
  ctx.fillText("0", 34, H - 38);
  ctx.fillText(omMax.toFixed(1), W - 72, H - 38);
  ctx.fillText(ymax.toExponential(1), 45, 38);
}

function drawXY(ctx, xs, ys, xMax, yMax, color) {
  const W = ctx.canvas.width, H = ctx.canvas.height;
  const x0 = 42, y0 = 24, w = W - 64, h = H - 64;
  ctx.strokeStyle = color; ctx.lineWidth = 1.4; ctx.beginPath();
  let started = false;
  for (let i = 0; i < xs.length; i++) {
    if (xs[i] > xMax) break;
    const px = x0 + w * xs[i] / xMax;
    const py = y0 + h * (1 - ys[i] / yMax);
    if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py);
  }
  ctx.stroke();
}

function downloadCSV() {
  if (!st || !st.records.length) { message("no records to download"); return; }
  const keys = ["t", "x_L1", "phi1_L1", "phi2_L1", "dphi1_L1", "dphi2_L1", "pi1_L1", "pi2_L1",
    "x_L2", "phi1_L2", "phi2_L2", "dphi1_L2", "dphi2_L2", "pi1_L2", "pi2_L2",
    "x_R1", "phi1_R1", "phi2_R1", "dphi1_R1", "dphi2_R1", "pi1_R1", "pi2_R1",
    "x_R2", "phi1_R2", "phi2_R2", "dphi1_R2", "dphi2_R2", "pi1_R2", "pi2_R2"];
  let csv = keys.join(",") + "\n";
  for (const r of st.records) csv += keys.map(k => Number.isFinite(r[k]) ? r[k] : "").join(",") + "\n";
  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = `qed2_meson_scattering_records_t${st.time.toFixed(1)}.csv`;
  document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(url);
}

function resetDynamic() {
  if (!st) { buildAtom(); return; }
  const p = readParams();
  const sameGrid = p.N === st.p.N && Math.abs(p.dx - st.p.dx) < 1e-12 && Math.abs(p.L - st.p.L) < 1e-12 && p.g === st.p.g && p.k1 === st.p.k1 && p.k2 === st.p.k2;
  if (!sameGrid) { buildAtom(); return; }
  st.p = p;
  st.md = massData(p);
  st.gamma = makeSponge(p, st.x);
  st.probe = updateProbeIndices(p);
  st.phi1 = st.phi1Static.slice(); st.phi2 = st.phi2Static.slice();
  st.pi1 = new Float64Array(p.N); st.pi2 = new Float64Array(p.N);
  st.time = 0; st.step = 0; st.records = []; st.lastSpectrum = null; st.running = false;
  if (p.driveMode === "packet") addIncomingPacket(st);
  st.initialEnergy = totalEnergy(st);
  drawAll(); updateReadouts(); message("dynamic fields reset");
}

function loop() {
  if (st && st.running) {
    for (let n = 0; n < st.p.substeps; n++) stepSimulation(st);
    updateReadouts();
    drawAll();
  }
  requestAnimationFrame(loop);
}

ui.buildBtn.addEventListener("click", buildAtom);
ui.startBtn.addEventListener("click", () => { if (!st) buildAtom(); st.running = true; message("running"); });
ui.pauseBtn.addEventListener("click", () => { if (st) st.running = false; message("paused"); });
ui.resetBtn.addEventListener("click", resetDynamic);
ui.spectrumBtn.addEventListener("click", analyzeSpectrum);
ui.downloadBtn.addEventListener("click", downloadCSV);
ui.clearDataBtn.addEventListener("click", () => { if (st) { st.records = []; st.lastSpectrum = null; drawAll(); message("records cleared"); } });

buildAtom();
loop();
