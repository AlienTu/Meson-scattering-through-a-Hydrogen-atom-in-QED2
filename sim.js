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
  if (Math.abs(b) > 1e-14) { v1 = b; v2 = lamLight - a; }
  else { v1 = 0; v2 = 1; }
  const norm = Math.sqrt(v1 * v1 + v2 * v2) || 1;
  v1 /= norm; v2 /= norm;
  if (v2 < 0) { v1 = -v1; v2 = -v2; }
  return { mLight2: Math.max(lamLight, 0), mHeavy2: Math.max(lamHeavy, 0), mLight: Math.sqrt(Math.max(lamLight, 0)), mHeavy: Math.sqrt(Math.max(lamHeavy, 0)), vLight: [v1, v2] };
}

function indexAt(x, p) { return clamp(Math.round((x + 0.5 * p.L) / p.dx), 0, p.N - 1); }

function potentialGradient(phi1, phi2, p) {
  const common = 0.5 * p.g * p.g * (phi1 + phi2);
  return [common + SQRT_PI * p.k1 * Math.sin(TWO_SQRT_PI * phi1), common + SQRT_PI * p.k2 * Math.sin(TWO_SQRT_PI * phi2)];
}

function buildAtom() {
  const p = readParams();
  const x = makeGrid(p);
  const phi1 = new Float64Array(p.N), phi2 = new Float64Array(p.N), v1 = new Float64Array(p.N), v2 = new Float64Array(p.N);
  const left1 = 0.0, left2 = 0.0, right1 = SQRT_PI, right2 = -SQRT_PI, width = 3.0;
  for (let i = 0; i < p.N; i++) { const s = 0.5 * (1 + Math.tanh(x[i] / width)); phi1[i] = left1 + (right1 - left1) * s; phi2[i] = left2 + (right2 - left2) * s; }
  phi1[0] = left1; phi2[0] = left2; phi1[p.N - 1] = right1; phi2[p.N - 1] = right2;
  const relaxDt = Math.min(0.035, 0.25 * p.dx), damp = 1.8;
  let maxResidual = Infinity;
  for (let it = 0; it < p.relaxSteps; it++) {
    maxResidual = 0;
    for (let i = 1; i < p.N - 1; i++) {
      const lap1 = (phi1[i - 1] - 2 * phi1[i] + phi1[i + 1]) / (p.dx * p.dx);
      const lap2 = (phi2[i - 1] - 2 * phi2[i] + phi2[i + 1]) / (p.dx * p.dx);
      const [d1, d2] = potentialGradient(phi1[i], phi2[i], p);
      const f1 = lap1 - d1, f2 = lap2 - d2;
      v1[i] += relaxDt * (f1 - damp * v1[i]); v2[i] += relaxDt * (f2 - damp * v2[i]);
      maxResidual = Math.max(maxResidual, Math.abs(f1), Math.abs(f2));
    }
    for (let i = 1; i < p.N - 1; i++) { phi1[i] += relaxDt * v1[i]; phi2[i] += relaxDt * v2[i]; }
    phi1[0] = left1; phi2[0] = left2; phi1[p.N - 1] = right1; phi2[p.N - 1] = right2;
    v1[0] = 0; v2[0] = 0; v1[p.N - 1] = 0; v2[p.N - 1] = 0;
  }
  const md = massData(p);
  st = { p, x, md, phi1Static: phi1.slice(), phi2Static: phi2.slice(), phi1: phi1.slice(), phi2: phi2.slice(), pi1: new Float64Array(p.N), pi2: new Float64Array(p.N), force1: new Float64Array(p.N), force2: new Float64Array(p.N), gamma: makeSponge(p, x), time: 0, step: 0, running: false, built: true, records: [], initialEnergy: null, lastSpectrum: null, probe: updateProbeIndices(p) };
  if (p.driveMode === "packet") addIncomingPacket(st);
  st.initialEnergy = totalEnergy(st);
  updateReadouts(maxResidual); drawAll(); message(`atom built; relaxation residual ${maxResidual.toExponential(2)}`);
}

function makeSponge(p, x) {
  const gamma = new Float64Array(p.N);
  const denom = Math.max(0.5 * p.L - p.spongeStart, p.dx);
  for (let i = 0; i < p.N; i++) { const ax = Math.abs(x[i]); if (ax > p.spongeStart) { const s = (ax - p.spongeStart) / denom; gamma[i] = p.spongeStrength * s * s; } }
  return gamma;
}
function updateProbeIndices(p) { return { L1: indexAt(p.xLeftProbe, p), L2: indexAt(p.xLeftProbe + p.probeSep, p), R1: indexAt(p.xRightProbe, p), R2: indexAt(p.xRightProbe + p.probeSep, p) }; }
function addIncomingPacket(s) { const p = s.p, omega = p.omegaIn, k = Math.sqrt(Math.max(omega * omega - s.md.mLight2, 0)), x0 = p.sourceX, sig = p.packetWidth, [v1, v2] = s.md.vLight; if (k <= 0) { message("omega_in is below the light meson mass; packet is evanescent."); return; } for (let i = 1; i < p.N - 1; i++) { const y = s.x[i] - x0, env = Math.exp(-0.5 * y * y / (sig * sig)), phase = k * y, val = p.ampIn * env * Math.cos(phase), vel = p.ampIn * omega * env * Math.sin(phase); s.phi1[i] += v1 * val; s.phi2[i] += v2 * val; s.pi1[i] += v1 * vel; s.pi2[i] += v2 * vel; } }
function computeForces(s) { const p = s.p, f1 = s.force1, f2 = s.force2, [mv1, mv2] = s.md.vLight, continuous = p.driveMode === "continuous", on = continuous ? (1 - Math.exp(-s.time * s.time / 25.0)) : 0, drivePhase = Math.sin(p.omegaIn * s.time); for (let i = 1; i < p.N - 1; i++) { const lap1 = (s.phi1[i - 1] - 2 * s.phi1[i] + s.phi1[i + 1]) / (p.dx * p.dx), lap2 = (s.phi2[i - 1] - 2 * s.phi2[i] + s.phi2[i + 1]) / (p.dx * p.dx), [d1, d2] = potentialGradient(s.phi1[i], s.phi2[i], p); let src1 = 0, src2 = 0; if (continuous) { const y = s.x[i] - p.sourceX, mask = Math.exp(-0.5 * y * y / (p.sourceWidth * p.sourceWidth)), amp = p.ampIn * on * drivePhase * mask; src1 = mv1 * amp; src2 = mv2 * amp; } f1[i] = lap1 - d1 - s.gamma[i] * s.pi1[i] + src1; f2[i] = lap2 - d2 - s.gamma[i] * s.pi2[i] + src2; } f1[0] = 0; f2[0] = 0; f1[p.N - 1] = 0; f2[p.N - 1] = 0; }
function stepSimulation(s) { const p = s.p; computeForces(s); for (let i = 1; i < p.N - 1; i++) { s.pi1[i] += p.dt * s.force1[i]; s.pi2[i] += p.dt * s.force2[i]; s.phi1[i] += p.dt * s.pi1[i]; s.phi2[i] += p.dt * s.pi2[i]; } s.phi1[0] = s.phi1Static[0]; s.phi2[0] = s.phi2Static[0]; s.phi1[p.N - 1] = s.phi1Static[p.N - 1]; s.phi2[p.N - 1] = s.phi2Static[p.N - 1]; s.pi1[0] = 0; s.pi2[0] = 0; s.pi1[p.N - 1] = 0; s.pi2[p.N - 1] = 0; s.time += p.dt; s.step++; if (s.step % p.recordEvery === 0) recordProbes(s); }
function recordProbes(s) { const pr = s.probe, rec = { t: s.time }; for (const [name, idx] of Object.entries(pr)) { rec[`x_${name}`] = s.x[idx]; rec[`phi1_${name}`] = s.phi1[idx]; rec[`phi2_${name}`] = s.phi2[idx]; rec[`dphi1_${name}`] = s.phi1[idx] - s.phi1Static[idx]; rec[`dphi2_${name}`] = s.phi2[idx] - s.phi2Static[idx]; rec[`pi1_${name}`] = s.pi1[idx]; rec[`pi2_${name}`] = s.pi2[idx]; } s.records.push(rec); if (s.records.length > s.p.maxRecords) s.records.shift(); }
function totalEnergy(s) { const p = s.p; let E = 0; for (let i = 1; i < p.N - 1; i++) { const d1 = (s.phi1[i + 1] - s.phi1[i - 1]) / (2 * p.dx), d2 = (s.phi2[i + 1] - s.phi2[i - 1]) / (2 * p.dx), kin = 0.5 * (s.pi1[i] * s.pi1[i] + s.pi2[i] * s.pi2[i]), grad = 0.5 * (d1 * d1 + d2 * d2), pot = 0.25 * p.g * p.g * sqr(s.phi1[i] + s.phi2[i]) - 0.5 * p.k1 * Math.cos(TWO_SQRT_PI * s.phi1[i]) - 0.5 * p.k2 * Math.cos(TWO_SQRT_PI * s.phi2[i]); E += (kin + grad + pot) * p.dx; } return E; }
function updateReadouts(residual = null) { if (!st) return; ui.time.textContent = st.time.toFixed(3); ui.Nreadout.textContent = st.p.N; ui.dxReadout.textContent = st.p.dx.toFixed(5); ui.massReadout.textContent = `m_light=${st.md.mLight.toFixed(4)}, m_heavy=${st.md.mHeavy.toFixed(4)}`; const k = Math.sqrt(Math.max(st.p.omegaIn * st.p.omegaIn - st.md.mLight2, 0)); ui.kReadout.textContent = k.toFixed(4); const E = totalEnergy(st); ui.energyReadout.textContent = st.initialEnergy ? `${((E - st.initialEnergy) / Math.max(Math.abs(st.initialEnergy), 1)).toExponential(2)}` : "-"; if (residual !== null) ui.energyReadout.textContent += `, bg res ${residual.toExponential(1)}`; }
function drawAxes(ctx, W, H, ymin, ymax) { ctx.strokeStyle = "#30363d"; ctx.strokeRect(38, 18, W - 58, H - 48); ctx.fillStyle = "#8b949e"; ctx.font = "12px sans-serif"; ctx.fillText(ymax.toExponential(1), 42, 32); ctx.fillText(ymin.toExponential(1), 42, H - 18); }
function drawAll() { if (!st) return; drawFields(); drawTrace(ui.leftTraceCanvas, "L"); drawTrace(ui.rightTraceCanvas, "R"); if (st.lastSpectrum) { drawSpectrum(st.lastSpectrum); drawDirection(st.lastSpectrum); } }
function panel(ctx, x, y, w, h, title) { ctx.strokeStyle = "#30363d"; ctx.strokeRect(x, y, w, h); ctx.fillStyle = "#c9d1d9"; ctx.font = "13px sans-serif"; ctx.fillText(title, x + 8, y + 18); }
function drawCurve(ctx, arr, x0, y0, w, h, ymin, ymax, color, lineWidth = 1.5, dashed = false) { const N = arr.length; ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = lineWidth; if (dashed) ctx.setLineDash([5, 4]); ctx.beginPath(); for (let i = 0; i < N; i++) { const px = x0 + w * i / (N - 1), py = y0 + h * (1 - (arr[i] - ymin) / (ymax - ymin || 1)); if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py); } ctx.stroke(); ctx.restore(); }
function drawFields() { const c = ui.fieldCanvas, ctx = c.getContext("2d"), W = c.width, H = c.height; ctx.clearRect(0, 0, W, H); ctx.fillStyle = "#080b10"; ctx.fillRect(0, 0, W, H); const pad = 42, h = (H - 3 * pad) / 2, w = W - 2 * pad; panel(ctx, pad, pad, w, h, "full fields: phi1, phi2, static background dashed"); panel(ctx, pad, 2 * pad + h, w, h, "vibrations: phi - phi_static"); let ymin = Infinity, ymax = -Infinity, dymin = Infinity, dymax = -Infinity; const d1 = new Float64Array(st.p.N), d2 = new Float64Array(st.p.N); for (let i = 0; i < st.p.N; i++) { ymin = Math.min(ymin, st.phi1[i], st.phi2[i], st.phi1Static[i], st.phi2Static[i]); ymax = Math.max(ymax, st.phi1[i], st.phi2[i], st.phi1Static[i], st.phi2Static[i]); d1[i] = st.phi1[i] - st.phi1Static[i]; d2[i] = st.phi2[i] - st.phi2Static[i]; dymin = Math.min(dymin, d1[i], d2[i]); dymax = Math.max(dymax, d1[i], d2[i]); } const ypad = 0.08 * (ymax - ymin || 1); ymin -= ypad; ymax += ypad; const dabs = Math.max(Math.abs(dymin), Math.abs(dymax), 1e-4); dymin = -dabs; dymax = dabs; drawCurve(ctx, st.phi1Static, pad, pad, w, h, ymin, ymax, "#8b949e", 1, true); drawCurve(ctx, st.phi2Static, pad, pad, w, h, ymin, ymax, "#8b949e", 1, true); drawCurve(ctx, st.phi1, pad, pad, w, h, ymin, ymax, "#58a6ff", 1.8); drawCurve(ctx, st.phi2, pad, pad, w, h, ymin, ymax, "#f2cc60", 1.8); drawCurve(ctx, d1, pad, 2 * pad + h, w, h, dymin, dymax, "#58a6ff", 1.6); drawCurve(ctx, d2, pad, 2 * pad + h, w, h, dymin, dymax, "#f2cc60", 1.6); markX(ctx, st.p.xLeftProbe, pad, pad, w, 2 * h + pad, "#3fb950", "L"); markX(ctx, st.p.xRightProbe, pad, pad, w, 2 * h + pad, "#f85149", "R"); ctx.fillStyle = "#58a6ff"; ctx.fillText("phi1", pad + w - 90, pad + 20); ctx.fillStyle = "#f2cc60"; ctx.fillText("phi2", pad + w - 45, pad + 20); }
function markX(ctx, xval, x0, y0, w, h, color, label) { const px = x0 + w * (xval + 0.5 * st.p.L) / st.p.L; ctx.save(); ctx.strokeStyle = color; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(px, y0); ctx.lineTo(px, y0 + h); ctx.stroke(); ctx.setLineDash([]); ctx.fillStyle = color; ctx.fillText(label, px + 4, y0 + 14); ctx.restore(); }
function drawTrace(canvas, side) { const ctx = canvas.getContext("2d"), W = canvas.width, H = canvas.height; ctx.clearRect(0, 0, W, H); ctx.fillStyle = "#080b10"; ctx.fillRect(0, 0, W, H); const recs = st.records.slice(-1200); if (recs.length < 2) return; const key1 = `dphi1_${side}1`, key2 = `dphi2_${side}1`; let ymax = 1e-5; for (const r of recs) ymax = Math.max(ymax, Math.abs(r[key1]), Math.abs(r[key2])); const arr1 = recs.map(r => r[key1]), arr2 = recs.map(r => r[key2]); drawAxes(ctx, W, H, -ymax, ymax); drawCurve(ctx, arr1, 38, 18, W - 58, H - 48, -ymax, ymax, "#58a6ff", 1.5); drawCurve(ctx, arr2, 38, 18, W - 58, H - 48, -ymax, ymax, "#f2cc60", 1.5); }
function fftComplex(re, im) { const n = re.length; if (n <= 1) return; const evenRe = new Float64Array(n / 2), evenIm = new Float64Array(n / 2), oddRe = new Float64Array(n / 2), oddIm = new Float64Array(n / 2); for (let i = 0; i < n / 2; i++) { evenRe[i] = re[2 * i]; evenIm[i] = im[2 * i]; oddRe[i] = re[2 * i + 1]; oddIm[i] = im[2 * i + 1]; } fftComplex(evenRe, evenIm); fftComplex(oddRe, oddIm); for (let k = 0; k < n / 2; k++) { const ang = -2 * Math.PI * k / n, cr = Math.cos(ang), ci = Math.sin(ang), tr = cr * oddRe[k] - ci * oddIm[k], ti = cr * oddIm[k] + ci * oddRe[k]; re[k] = evenRe[k] + tr; im[k] = evenIm[k] + ti; re[k + n / 2] = evenRe[k] - tr; im[k + n / 2] = evenIm[k] - ti; } }
function nextPow2(n) { let p = 1; while (p < n) p <<= 1; return p; }
function complexAt(re, im, omega, dt) { const n = re.length; const k = Math.floor((omega * dt * n) / (2 * Math.PI)); if (k < 0 || k >= n) return [0, 0]; return [re[k], im[k]]; }
function analyzeSpectrum() { if (!st || st.records.length < 32) { message("not enough records for FFT"); return; } const recs = st.records, n0 = nextPow2(recs.length), n = Math.min(n0, 8192), start = recs.length - n, dt = (recs[1].t - recs[0].t), [v1, v2] = st.md.vLight; const left = new Float64Array(n), right = new Float64Array(n), l2 = new Float64Array(n), r2 = new Float64Array(n); for (let i = 0; i < n; i++) { const r = recs[start + i], win = 0.5 * (1 - Math.cos(2 * Math.PI * i / (n - 1))); left[i] = win * (v1 * r.dphi1_L1 + v2 * r.dphi2_L1); right[i] = win * (v1 * r.dphi1_R1 + v2 * r.dphi2_R1); l2[i] = win * (v1 * r.dphi1_L2 + v2 * r.dphi2_L2); r2[i] = win * (v1 * r.dphi1_R2 + v2 * r.dphi2_R2); } const imL = new Float64Array(n), imR = new Float64Array(n), imL2 = new Float64Array(n), imR2 = new Float64Array(n); fftComplex(left, imL); fftComplex(right, imR); fftComplex(l2, imL2); fftComplex(r2, imR2); const half = n / 2, omega = [], leftPower = [], rightPower = [], leftMovingLeft = [], rightMovingLeft = [], leftMovingRight = [], rightMovingRight = []; const dxProbe = st.p.probeSep; for (let k = 1; k < half; k++) { const om = 2 * Math.PI * k / (n * dt); omega.push(om); const pL = left[k] * left[k] + imL[k] * imL[k], pR = right[k] * right[k] + imR[k] * imR[k]; leftPower.push(pL); rightPower.push(pR); const kk = Math.sqrt(Math.max(om * om - st.md.mLight2, 0)); if (kk * dxProbe < 1e-8) { leftMovingLeft.push(0); rightMovingLeft.push(pL); leftMovingRight.push(0); rightMovingRight.push(pR); continue; } const c = Math.cos(kk * dxProbe), s = Math.sin(kk * dxProbe); const det = 2 * s; function split(Are, Aim, Bre, Bim) { const Rre = (Bim - Aim * c - Are * s) / det; const Rim = (-Bre + Are * c - Aim * s) / det; const Lre = Are - Rre, Lim = Aim - Rim; return [Lre * Lre + Lim * Lim, Rre * Rre + Rim * Rim]; } const [LpL, RpL] = split(left[k], imL[k], l2[k], imL2[k]), [LpR, RpR] = split(right[k], imR[k], r2[k], imR2[k]); leftMovingLeft.push(LpL); rightMovingLeft.push(RpL); leftMovingRight.push(LpR); rightMovingRight.push(RpR); } st.lastSpectrum = { omega, leftPower, rightPower, leftMovingLeft, rightMovingLeft, leftMovingRight, rightMovingRight }; drawSpectrum(st.lastSpectrum); drawDirection(st.lastSpectrum); message("spectrum computed"); }
function drawSpectrum(spec) { const ctx = ui.spectrumCanvas.getContext("2d"), W = ui.spectrumCanvas.width, H = ui.spectrumCanvas.height; ctx.clearRect(0, 0, W, H); ctx.fillStyle = "#080b10"; ctx.fillRect(0, 0, W, H); if (!spec || !spec.omega.length) return; const omMax = Math.min(12, spec.omega[spec.omega.length - 1]); const vals = []; for (let i = 0; i < spec.omega.length; i++) if (spec.omega[i] <= omMax) vals.push(spec.leftPower[i], spec.rightPower[i]); const ymax = Math.max(...vals, 1e-16); drawSpectralAxes(ctx, W, H, omMax, ymax, "raw projected spectra"); drawXY(ctx, spec.omega, spec.leftPower, omMax, ymax, "#3fb950"); drawXY(ctx, spec.omega, spec.rightPower, omMax, ymax, "#f85149"); ctx.fillStyle = "#3fb950"; ctx.fillText("left probe", 52, 42); ctx.fillStyle = "#f85149"; ctx.fillText("right probe", 150, 42); }
function drawDirection(spec) { const ctx = ui.directionCanvas.getContext("2d"), W = ui.directionCanvas.width, H = ui.directionCanvas.height; ctx.clearRect(0, 0, W, H); ctx.fillStyle = "#080b10"; ctx.fillRect(0, 0, W, H); if (!spec || !spec.omega.length) return; const omMax = Math.min(12, spec.omega[spec.omega.length - 1]); const vals = []; for (let i = 0; i < spec.omega.length; i++) if (spec.omega[i] <= omMax) vals.push(spec.leftMovingLeft[i], spec.rightMovingRight[i], spec.rightMovingLeft[i]); const ymax = Math.max(...vals, 1e-16); drawSpectralAxes(ctx, W, H, omMax, ymax, "directional power, two-probe estimate"); drawXY(ctx, spec.omega, spec.leftMovingLeft, omMax, ymax, "#f2cc60"); drawXY(ctx, spec.omega, spec.rightMovingRight, omMax, ymax, "#58a6ff"); drawXY(ctx, spec.omega, spec.rightMovingLeft, omMax, ymax, "#8b949e"); ctx.fillStyle = "#f2cc60"; ctx.fillText("reflected at left", 52, 42); ctx.fillStyle = "#58a6ff"; ctx.fillText("transmitted at right", 170, 42); ctx.fillStyle = "#8b949e"; ctx.fillText("incident at left", 330, 42); }
function drawSpectralAxes(ctx, W, H, omMax, ymax, title) { ctx.strokeStyle = "#30363d"; ctx.strokeRect(42, 24, W - 64, H - 64); ctx.fillStyle = "#c9d1d9"; ctx.font = "12px sans-serif"; ctx.fillText(title, 52, 18); ctx.fillStyle = "#8b949e"; ctx.fillText("omega", W - 58, H - 14); ctx.fillText("0", 34, H - 38); ctx.fillText(omMax.toFixed(1), W - 72, H - 38); ctx.fillText(ymax.toExponential(1), 45, 38); }
function drawXY(ctx, xs, ys, xMax, yMax, color) { const W = ctx.canvas.width, H = ctx.canvas.height, x0 = 42, y0 = 24, w = W - 64, h = H - 64; ctx.strokeStyle = color; ctx.lineWidth = 1.4; ctx.beginPath(); let started = false; for (let i = 0; i < xs.length; i++) { if (xs[i] > xMax) break; const px = x0 + w * xs[i] / xMax, py = y0 + h * (1 - ys[i] / yMax); if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py); } ctx.stroke(); }
function downloadCSV() { if (!st || !st.records.length) { message("no records to download"); return; } const keys = ["t", "x_L1", "phi1_L1", "phi2_L1", "dphi1_L1", "dphi2_L1", "pi1_L1", "pi2_L1", "x_L2", "phi1_L2", "phi2_L2", "dphi1_L2", "dphi2_L2", "pi1_L2", "pi2_L2", "x_R1", "phi1_R1", "phi2_R1", "dphi1_R1", "dphi2_R1", "pi1_R1", "pi2_R1", "x_R2", "phi1_R2", "phi2_R2", "dphi1_R2", "dphi2_R2", "pi1_R2", "pi2_R2"]; let csv = keys.join(",") + "\n"; for (const r of st.records) csv += keys.map(k => Number.isFinite(r[k]) ? r[k] : "").join(",") + "\n"; const blob = new Blob([csv], { type: "text/csv" }), url = URL.createObjectURL(blob), a = document.createElement("a"); a.href = url; a.download = `qed2_meson_scattering_records_t${st.time.toFixed(1)}.csv`; document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(url); }
function resetDynamic() { if (!st) { buildAtom(); return; } const p = readParams(), sameGrid = p.N === st.p.N && Math.abs(p.dx - st.p.dx) < 1e-12 && Math.abs(p.L - st.p.L) < 1e-12 && p.g === st.p.g && p.k1 === st.p.k1 && p.k2 === st.p.k2; if (!sameGrid) { buildAtom(); return; } st.p = p; st.md = massData(p); st.gamma = makeSponge(p, st.x); st.probe = updateProbeIndices(p); st.phi1 = st.phi1Static.slice(); st.phi2 = st.phi2Static.slice(); st.pi1 = new Float64Array(p.N); st.pi2 = new Float64Array(p.N); st.time = 0; st.step = 0; st.records = []; st.lastSpectrum = null; st.running = false; if (p.driveMode === "packet") addIncomingPacket(st); st.initialEnergy = totalEnergy(st); drawAll(); updateReadouts(); message("dynamic fields reset"); }
function loop() { if (st && st.running) { for (let n = 0; n < st.p.substeps; n++) stepSimulation(st); updateReadouts(); drawAll(); } requestAnimationFrame(loop); }
ui.buildBtn.addEventListener("click", buildAtom); ui.startBtn.addEventListener("click", () => { if (!st) buildAtom(); st.running = true; message("running"); }); ui.pauseBtn.addEventListener("click", () => { if (st) st.running = false; message("paused"); }); ui.resetBtn.addEventListener("click", resetDynamic); ui.spectrumBtn.addEventListener("click", analyzeSpectrum); ui.downloadBtn.addEventListener("click", downloadCSV); ui.clearDataBtn.addEventListener("click", () => { if (st) { st.records = []; st.lastSpectrum = null; drawAll(); message("records cleared"); } });
message("ready; click Build atom first, then Solve bound mode");
loop();
