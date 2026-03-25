/*
 * Ported from http://svn.r-project.org/R/trunk/src/nmath/qnorm.c
 *
 * Mathlib : A C Library of Special Functions
 * Copyright (C) 1998       Ross Ihaka
 * Copyright (C) 2000--2005 The R Core Team
 * based on AS 111 (C) 1977 Royal Statistical Society
 * and   on AS 241 (C) 1988 Royal Statistical Society
 *
 * Ported to JavaScript by Ryosuke Niwa (https://github.com/rniwa/js-shapiro-wilk)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

function normalQuantile(p: number): number {
  let r: number, val: number;

  const q = p - 0.5;

  if (0.075 <= p && p <= 0.925) {
    r = 0.180625 - q * q;
    val =
      (q *
        (((((((r * 2509.0809287301226727 +
          33430.575583588128105) *
          r +
          67265.770927008700853) *
          r +
          45921.953931549871457) *
          r +
          13731.693765509461125) *
          r +
          1971.5909503065514427) *
          r +
          133.14166789178437745) *
          r +
          3.387132872796366608)) /
      (((((((r * 5226.495278852854561 +
        28729.085735721942674) *
        r +
        39307.89580009271061) *
        r +
        21213.794301586595867) *
        r +
        5394.1960214247511077) *
        r +
        687.1870074920579083) *
        r +
        42.313330701600911252) *
        r +
        1);
  } else {
    r = q > 0 ? 1 - p : p;
    r = Math.sqrt(-Math.log(r));

    if (r <= 5.0) {
      r -= 1.6;
      val =
        (((((((r * 7.7454501427834140764e-4 +
          0.0227238449892691845833) *
          r +
          0.24178072517745061177) *
          r +
          1.27045825245236838258) *
          r +
          3.64784832476320460504) *
          r +
          5.7694972214606914055) *
          r +
          4.6303378461565452959) *
          r +
          1.42343711074968357734) /
        (((((((r * 1.05075007164441684324e-9 +
          5.475938084995344946e-4) *
          r +
          0.0151986665636164571966) *
          r +
          0.14810397642748007459) *
          r +
          0.68976733498510000455) *
          r +
          1.6763848301838038494) *
          r +
          2.05319162663775882187) *
          r +
          1);
    } else {
      r -= 5.0;
      val =
        (((((((r * 2.01033439929228813265e-7 +
          2.71155556874348757815e-5) *
          r +
          0.0012426609473880784386) *
          r +
          0.026532189526576123093) *
          r +
          0.29656057182850489123) *
          r +
          1.7848265399172913358) *
          r +
          5.4637849111641143699) *
          r +
          6.6579046435011037772) /
        (((((((r * 2.04426310338993978564e-15 +
          1.4215117583164458887e-7) *
          r +
          1.8463183175100546818e-5) *
          r +
          7.868691311456132591e-4) *
          r +
          0.0148753612908506148525) *
          r +
          0.13692988092273580531) *
          r +
          0.59983220655588793769) *
          r +
          1.0);
    }

    if (q < 0.0) val = -val;
  }
  return val;
}

/*
 * Ported from http://svn.r-project.org/R/trunk/src/library/stats/src/swilk.c
 *
 * R : A Computer Language for Statistical Data Analysis
 * Copyright (C) 2000-12   The R Core Team.
 *
 * Based on Applied Statistics algorithms AS181, R94
 *   (C) Royal Statistical Society 1982, 1995
 *
 * Ported to JavaScript by Ryosuke Niwa (https://github.com/rniwa/js-shapiro-wilk)
 */

function poly(cc: number[], nord: number, x: number): number {
  let ret_val = cc[0];
  if (nord > 1) {
    let p = x * cc[nord - 1];
    for (let j = nord - 2; j > 0; j--) p = (p + cc[j]) * x;
    ret_val += p;
  }
  return ret_val;
}

function sign(x: number): number {
  if (x === 0) return 0;
  return x > 0 ? 1 : -1;
}

// Upper-tail normal CDF: P(X > x) for X ~ N(m, s^2)
// Uses erfc approximation (Horner form, max error ~1.5e-7)
function pnormUpper(x: number, m: number, s: number): number {
  const z = (x - m) / (s * Math.SQRT2);
  if (z < 0) {
    const t = 1 / (1 - 0.3275911 * z);
    const y =
      t *
      Math.exp(-z * z - 1.26551223 +
        t * (1.00002368 +
          t * (0.37409196 +
            t * (0.09678418 +
              t * (-0.18628806 +
                t * (0.27886807 +
                  t * (-1.13520398 +
                    t * (1.48851587 +
                      t * (-0.82215223 + t * 0.17087294)))))))));
    return 1 - y / 2;
  }
  const t = 1 / (1 + 0.3275911 * z);
  const y =
    t *
    Math.exp(-z * z - 1.26551223 +
      t * (1.00002368 +
        t * (0.37409196 +
          t * (0.09678418 +
            t * (-0.18628806 +
              t * (0.27886807 +
                t * (-1.13520398 +
                  t * (1.48851587 +
                    t * (-0.82215223 + t * 0.17087294)))))))));
  return y / 2;
}

function shapiroWilkW(x: number[]): { w: number; pvalue: number } | null {
  x = [...x].sort((a, b) => a - b);
  const n = x.length;
  if (n < 3) return null;

  const nn2 = Math.floor(n / 2);
  const a: number[] = new Array(nn2 + 1);

  const small = 1e-19;
  const g = [-2.273, 0.459];
  const c1 = [0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056];
  const c2 = [0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633];
  const c3 = [0.544, -0.39978, 0.025054, -6.714e-4];
  const c4 = [1.3822, -0.77857, 0.062767, -0.0020322];
  const c5 = [-1.5861, -0.31082, -0.083751, 0.0038915];
  const c6 = [-0.4803, -0.082676, 0.0030302];

  let pw = 1;
  const an = n;

  if (n === 3) {
    a[1] = 0.70710678;
  } else {
    const an25 = an + 0.25;
    let summ2 = 0.0;
    for (let i = 1; i <= nn2; i++) {
      a[i] = normalQuantile((i - 0.375) / an25);
      summ2 += a[i] * a[i];
    }
    summ2 *= 2;
    const ssumm2 = Math.sqrt(summ2);
    const rsn = 1 / Math.sqrt(an);
    const a1 = poly(c1, 6, rsn) - a[1] / ssumm2;

    let i1: number;
    let fac: number;
    if (n > 5) {
      i1 = 3;
      const a2 = -a[2] / ssumm2 + poly(c2, 6, rsn);
      fac = Math.sqrt(
        (summ2 - 2 * a[1] * a[1] - 2 * a[2] * a[2]) /
          (1 - 2 * a1 * a1 - 2 * a2 * a2),
      );
      a[2] = a2;
    } else {
      i1 = 2;
      fac = Math.sqrt(
        (summ2 - 2 * a[1] * a[1]) / (1 - 2 * a1 * a1),
      );
    }
    a[1] = a1;
    for (let i = i1; i <= nn2; i++) a[i] /= -fac;
  }

  const range = x[n - 1] - x[0];
  if (range < small) return null;

  let xx = x[0] / range;
  let sx = xx;
  let sa = -a[1];
  for (let i = 1, j = n - 1; i < n; j--) {
    const xi = x[i] / range;
    if (xx - xi > small) return null;
    sx += xi;
    i++;
    if (i !== j) sa += sign(i - j) * a[Math.min(i, j)];
    xx = xi;
  }
  if (n > 5000) return null;

  sa /= n;
  sx /= n;
  let ssa = 0, ssx = 0, sax = 0;
  for (let i = 0, j = n - 1; i < n; i++, j--) {
    const asa = i !== j ? sign(i - j) * a[1 + Math.min(i, j)] - sa : -sa;
    const xsx = x[i] / range - sx;
    ssa += asa * asa;
    ssx += xsx * xsx;
    sax += asa * xsx;
  }

  const ssassx = Math.sqrt(ssa * ssx);
  const w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
  const w = 1 - w1;

  if (n === 3) {
    const pi6 = 1.90985931710274;
    const stqr = 1.04719755119660;
    pw = pi6 * (Math.asin(Math.sqrt(w)) - stqr);
    if (pw < 0) pw = 0;
    return { w, pvalue: pw };
  }

  const y = Math.log(w1);
  const logn = Math.log(an);
  let m: number, s: number;

  if (n <= 11) {
    const gamma = poly(g, 2, an);
    if (y >= gamma) return { w, pvalue: 1e-99 };
    const y2 = -Math.log(gamma - y);
    m = poly(c3, 4, an);
    s = Math.exp(poly(c4, 4, an));
    pw = pnormUpper(y2, m, s);
  } else {
    m = poly(c5, 4, logn);
    s = Math.exp(poly(c6, 3, logn));
    pw = pnormUpper(y, m, s);
  }

  return { w, pvalue: pw };
}

function iqrFilter(data: number[]): number[] {
  if (data.length < 4) return data;
  const sorted = [...data].sort((a, b) => a - b);
  const n = sorted.length;
  const q1 = sorted[Math.floor(n * 0.25)];
  const q3 = sorted[Math.floor(n * 0.75)];
  const iqr = q3 - q1;
  const low = q1 - 1.5 * iqr;
  const high = q3 + 1.5 * iqr;
  return sorted.filter((x) => x >= low && x <= high);
}

export function shapiroWilkTest(
  data: number[],
): { w: number; pvalue: number } | null {
  return shapiroWilkW(iqrFilter(data));
}
