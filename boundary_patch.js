"use strict";

// Open-boundary patch for sim.js.
// Static background is still found with fixed topological end values.
// Real-time evolution then applies an outgoing Sommerfeld condition to
// perturbations delta phi = phi - phi_static.

(function () {
  const oldReadParams = readParams;
  readParams = function () {
    const p = oldReadParams();
    const el = document.getElementById("boundaryMode");
    p.boundaryMode = el ? el.value : "radiation";
    return p;
  };

  const oldMakeSponge = makeSponge;
  makeSponge = function (p, x) {
    if (p.boundaryMode === "radiation") return new Float64Array(p.N);
    return oldMakeSponge(p, x);
  };

  function applyDirichletBoundary(s) {
    const N = s.p.N;
    s.phi1[0] = s.phi1Static[0];
    s.phi2[0] = s.phi2Static[0];
    s.phi1[N - 1] = s.phi1Static[N - 1];
    s.phi2[N - 1] = s.phi2Static[N - 1];
    s.pi1[0] = 0;
    s.pi2[0] = 0;
    s.pi1[N - 1] = 0;
    s.pi2[N - 1] = 0;
  }

  function applyRadiationBoundary(s, oldL1, oldL2, oldR1, oldR2) {
    const p = s.p;
    const N = p.N;
    const c = p.dt / p.dx;

    // Left boundary: outgoing-to-left perturbation, (dt - dx) delta phi = 0.
    const dL1 = (oldL1 - s.phi1Static[0]) + c * ((s.phi1[1] - s.phi1Static[1]) - (oldL1 - s.phi1Static[0]));
    const dL2 = (oldL2 - s.phi2Static[0]) + c * ((s.phi2[1] - s.phi2Static[1]) - (oldL2 - s.phi2Static[0]));

    // Right boundary: outgoing-to-right perturbation, (dt + dx) delta phi = 0.
    const dR1 = (oldR1 - s.phi1Static[N - 1]) - c * ((oldR1 - s.phi1Static[N - 1]) - (s.phi1[N - 2] - s.phi1Static[N - 2]));
    const dR2 = (oldR2 - s.phi2Static[N - 1]) - c * ((oldR2 - s.phi2Static[N - 1]) - (s.phi2[N - 2] - s.phi2Static[N - 2]));

    s.phi1[0] = s.phi1Static[0] + dL1;
    s.phi2[0] = s.phi2Static[0] + dL2;
    s.phi1[N - 1] = s.phi1Static[N - 1] + dR1;
    s.phi2[N - 1] = s.phi2Static[N - 1] + dR2;

    s.pi1[0] = (s.phi1[0] - oldL1) / p.dt;
    s.pi2[0] = (s.phi2[0] - oldL2) / p.dt;
    s.pi1[N - 1] = (s.phi1[N - 1] - oldR1) / p.dt;
    s.pi2[N - 1] = (s.phi2[N - 1] - oldR2) / p.dt;
  }

  stepSimulation = function (s) {
    const p = s.p;
    const N = p.N;
    const oldL1 = s.phi1[0];
    const oldL2 = s.phi2[0];
    const oldR1 = s.phi1[N - 1];
    const oldR2 = s.phi2[N - 1];

    computeForces(s);
    for (let i = 1; i < N - 1; i++) {
      s.pi1[i] += p.dt * s.force1[i];
      s.pi2[i] += p.dt * s.force2[i];
      s.phi1[i] += p.dt * s.pi1[i];
      s.phi2[i] += p.dt * s.pi2[i];
    }

    if (p.boundaryMode === "dirichlet") {
      applyDirichletBoundary(s);
    } else {
      applyRadiationBoundary(s, oldL1, oldL2, oldR1, oldR2);
    }

    s.time += p.dt;
    s.step++;
    if (s.step % p.recordEvery === 0) recordProbes(s);
  };

  if (typeof st !== "undefined" && st) {
    st.p = readParams();
    st.gamma = makeSponge(st.p, st.x);
    st.probe = updateProbeIndices(st.p);
    if (typeof message === "function") message("open radiation boundary patch loaded");
  }
})();
