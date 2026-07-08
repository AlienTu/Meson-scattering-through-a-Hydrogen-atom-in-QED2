# Meson scattering through a hydrogen atom in QED2

A browser-based semiclassical simulation for meson scattering through a one-atom background in the bosonized two-flavour massive Schwinger model.

The model evolved in flavour variables is

```math
\partial_t^2\phi_i-\partial_x^2\phi_i+\frac{g^2}{2}(\phi_1+\phi_2)+\sqrt{\pi}\,\kappa_i\sin(2\sqrt{\pi}\phi_i)=J_i(t,x),\quad i=1,2.
```

The app:

1. relaxes a static one-atom background near `x=0` with neutral topological boundary values,
2. injects either a continuous left source or a finite right-moving packet,
3. evolves the full nonlinear classical field equations with absorbing sponge layers,
4. records `phi1`, `phi2` near `x=-40` and `x=+40`,
5. computes Fourier spectra and an approximate right/left directional decomposition using two nearby probes.

Default parameters are `L=100`, `dx=0.1`, `g=1.5`, `kappa2=1`, and a heavy `kappa1=20`.

## Notes

- Dirichlet values are only used to hold the static one-atom sector. The simulation uses sponge layers near the boundaries to suppress reflections.
- The continuous source is not a perfectly one-way boundary condition. It emits both directions locally, but the left-moving component is damped by the left sponge. For a cleaner single-shot scattering experiment, use the finite packet mode.
- Directional spectra are asymptotic linear diagnostics. They are most meaningful when the probe points are outside the nonlinear atom core and before the sponge region.
