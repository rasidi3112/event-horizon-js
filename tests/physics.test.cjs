// Physics regression tests. Run with: node --test tests/*.test.cjs
const test = require('node:test');
const assert = require('node:assert/strict');
const path = require('node:path');

const {
    CONSTANTS, BlackHolePhysics, GravitationalWaveSimulator,
    formatTemperature, formatTime
} = require(path.join(__dirname, '..', 'physics.js'));

const relClose = (actual, expected, tol, msg) => {
    const rel = Math.abs(actual - expected) / Math.abs(expected);
    assert.ok(rel <= tol, `${msg}: got ${actual}, expected ${expected} (rel. error ${rel})`);
};

test('Schwarzschild radius of 1 M☉ ≈ 2953.25 m', () => {
    const bh = new BlackHolePhysics(1, 0);
    assert.ok(Math.abs(bh.rs - 2953.25) <= 0.1, `rs = ${bh.rs}`);
});

test('Hawking temperature of 1 M☉ ≈ 6.17e-8 K', () => {
    const bh = new BlackHolePhysics(1, 0);
    relClose(bh.hawkingTemperature, 6.17e-8, 0.005, 'T_H');
});

test('surface gravity for a* = 0 equals c⁴/(4GM)', () => {
    const bh = new BlackHolePhysics(10, 0);
    const expected = CONSTANTS.c ** 4 / (4 * CONSTANTS.G * bh.M);
    relClose(bh.surfaceGravity, expected, 1e-9, 'κ');
});

test('Kerr entropy ratio, equatorial ergosphere and ISCO at a* = 0.998', () => {
    const a = 0.998;
    const kerr = new BlackHolePhysics(1, a);
    const schw = new BlackHolePhysics(1, 0);
    relClose(kerr.entropy / schw.entropy, (1 + Math.sqrt(1 - a ** 2)) / 2, 1e-12, 'S_Kerr/S_Schw');
    relClose(kerr.rErgo, kerr.rs, 1e-12, 'r_ergo(equator)');
    relClose(kerr.ergosphereRadius(0), kerr.rPlus, 1e-12, 'r_ergo(pole) = r₊');
    relClose(schw.isco, 6 * schw.rg, 1e-12, 'ISCO(a*=0)');
    relClose(kerr.isco / kerr.rg, 1.237, 1e-3, 'ISCO(a*=0.998)');
});

test('Hawking temperature equals ħκ/(2πck_B) for several spins', () => {
    for (const a of [0, 0.5, 0.9, 0.998]) {
        const bh = new BlackHolePhysics(5, a);
        const fromKappa = CONSTANTS.hbar * bh.surfaceGravity / (2 * Math.PI * CONSTANTS.c * CONSTANTS.kB);
        relClose(bh.hawkingTemperature, fromKappa, 1e-12, `T_H(a*=${a})`);
    }
});

test('OTOC is in [0, 1], non-increasing and deterministic', () => {
    const bh = new BlackHolePhysics(10, 0.5);
    const tMax = 3 * bh.scramblingTime;
    let previous = Infinity;
    for (let i = 0; i <= 200; i++) {
        const t = (i / 200) * tMax;
        const F = bh.getOTOCDecay(t);
        assert.ok(F >= 0 && F <= 1, `F(${t}) = ${F}`);
        assert.ok(F <= previous, `F not monotone at t = ${t}`);
        assert.equal(F, bh.getOTOCDecay(t));
        previous = F;
    }
});

test('naive Page curve peaks at ≈ S₀/2 near x = 1 − 2^(−3/2)', () => {
    const bh = new BlackHolePhysics(1, 0);
    const tEvap = bh.evaporationTime;
    let best = { x: 0, s: -Infinity };
    for (let i = 0; i <= 10000; i++) {
        const x = i / 10000;
        const s = bh.getPageCurveEntropy(x * tEvap, tEvap);
        if (s > best.s) best = { x, s };
    }
    relClose(best.s, bh.entropy / 2, 1e-3, 'S_Page,max');
    assert.ok(Math.abs(best.x - (1 - 2 ** -1.5)) < 1e-3, `peak at x = ${best.x}`);
    relClose(bh.pageTime / tEvap, 0.646, 1e-3, 'pageTime / t_evap');
});

test('GW frequency for 36 + 29 M☉ at τ = 0.1 s ≈ 44.5 Hz', () => {
    const gw = new GravitationalWaveSimulator(36, 29, 1.4e9);
    relClose(gw.getGWFrequency(0, 0.1), 44.5, 0.01, 'f_GW');
});

test('eikonal QNM l = 2, n = 0 gives Mω_R = 2.5/√27', () => {
    const bh = new BlackHolePhysics(1, 0);
    const omega = bh.getQNMFrequency(2, 0);
    const geometricUnit = CONSTANTS.G * bh.M / CONSTANTS.c ** 3; // GM/c³ (s)
    relClose(omega.real * geometricUnit, 2.5 / Math.sqrt(27), 1e-12, 'Mω_R');
    relClose(omega.imag * geometricUnit, -0.5 / Math.sqrt(27), 1e-12, 'Mω_I');
});

test('formatters stay readable for extreme values', () => {
    const t = formatTemperature(5.76e-18);
    assert.ok(!t.includes('0.00'), `formatTemperature gave "${t}"`);
    const s = formatTime(1e80);
    assert.ok(s.length <= 30, `formatTime gave "${s}" (${s.length} chars)`);
});
