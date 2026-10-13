// Same normalized frames as XFOIL BLPRV/BLKIN; the pinned full reference uses HVRAT=0.
export function edgeState(ui, mach, theta, re, h) {
  const msq = mach * mach;
  const beta = Math.sqrt(1 - msq);
  const tk = msq / (1 + beta) ** 2;
  const u = ui * (1 - tk) / (1 - tk * ui * ui);
  const hstinv = .4 * msq / (1 + .2 * msq);
  const he = 1 - .5 * u * u * hstinv;
  const localMsq = u * u * hstinv / (.4 * he);
  const density = ((1 + .2 * msq) / (1 + .2 * localMsq)) ** 2.5;
  const heInf = 1 - .5 * hstinv;
  const reBl = re * Math.sqrt(heInf ** 3) / heInf;
  const viscosity = Math.sqrt(he ** 3) / he / reBl;
  return { u, localMsq, density, rtheta: density * u * theta / viscosity,
    hk: (h - .29 * localMsq) / (1 + .113 * localMsq) };
}

