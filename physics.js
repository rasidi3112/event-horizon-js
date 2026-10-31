// === Pure physics engine (no DOM) ===
// Notation: a* = dimensionless spin (this.a), rg = GM/c² (m), aLen = a*·rg (m).

// === Physical Constants (SI Units) ===
// CODATA 2018 exact/recommended values; GM_sun is the IAU 2015 nominal value (Resolution B3).
const CONSTANTS = (() => {
    const c = 299792458;              // Speed of light (m/s), exact
    const G = 6.67430e-11;            // Gravitational constant (m³/kg/s²), CODATA 2018
    const hbar = 1.054571817e-34;     // Reduced Planck constant (J·s), CODATA 2018
    const kB = 1.380649e-23;          // Boltzmann constant (J/K), exact
    const GM_sun = 1.3271244e20;      // Nominal solar mass parameter (m³/s²), IAU 2015
    return {
        c, G, hbar, kB, GM_sun,
        M_sun: GM_sun / G,                  // Solar mass (kg) = GM_sun / G
        lP: Math.sqrt(hbar * G / c ** 3),   // Planck length (m) = √(ħG/c³)
        tP: Math.sqrt(hbar * G / c ** 5),   // Planck time (s) = √(ħG/c⁵)
        mP: Math.sqrt(hbar * c / G),        // Planck mass (kg) = √(ħc/G)
        ly: 9460730472580800,               // Light year (m), IAU definition (Julian year × c)
        year: 31557600,                     // Julian year (s) = 365.25 d
    };
})();

// === Known Black Holes Database (from observations) ===
const KNOWN_BLACK_HOLES = [
    { name: "Sagittarius A*", mass: 4.154e6, distance: 26673, type: "Supermassive", spin: 0.9, source: "EHT 2022" },
    { name: "M87*", mass: 6.5e9, distance: 53.49e6, type: "Supermassive", spin: 0.9, source: "EHT 2019" },
    { name: "Cygnus X-1", mass: 21.2, distance: 6070, type: "Stellar", spin: 0.998, source: "LIGO" },
    { name: "GW150914 Primary", mass: 35.6, distance: 1.4e9, type: "Stellar", spin: 0.31, source: "LIGO 2015" },
    { name: "GW150914 Secondary", mass: 30.6, distance: 1.4e9, type: "Stellar", spin: 0.46, source: "LIGO 2015" },
    { name: "GW150914 Remnant", mass: 63.1, distance: 1.4e9, type: "Stellar", spin: 0.69, source: "LIGO 2015" },
    { name: "TON 618", mass: 6.6e10, distance: 1.04e10, type: "Ultramassive", spin: 0.7, source: "Optical" },
    { name: "NGC 1277", mass: 1.7e10, distance: 2.2e8, type: "Supermassive", spin: 0.85, source: "Hubble" },
    { name: "Holmberg 15A*", mass: 4e10, distance: 7e8, type: "Ultramassive", spin: 0.6, source: "VLT" },
    { name: "Phoenix A*", mass: 1e11, distance: 5.8e9, type: "Ultramassive", spin: 0.5, source: "Radio" },
    { name: "GRS 1915+105", mass: 12.4, distance: 28000, type: "Stellar", spin: 0.98, source: "X-ray" },
    { name: "LB-1 Primary", mass: 68, distance: 15000, type: "Stellar", spin: 0.7, source: "Spectroscopy" },
];

// === Quantum Chaos Black Hole Physics Engine ===
class BlackHolePhysics {
    constructor(massInSolarMasses, spin = 0) {
        this.M_solar = massInSolarMasses;
        this.M = massInSolarMasses * CONSTANTS.M_sun;
        this.a = Math.min(Math.max(spin, 0), 0.998); // Kerr spin parameter a* (0 to 0.998, Thorne limit)
    }

    // Gravitational radius rg = GM/c² (m)
    get rg() {
        return (CONSTANTS.G * this.M) / (CONSTANTS.c ** 2);
    }

    // Schwarzschild radius rs = 2GM/c² (m)
    get rs() {
        return 2 * this.rg;
    }

    // Kerr parameter as a length aLen = a*·rg = J/(Mc) (m)
    get aLen() {
        return this.a * this.rg;
    }

    // Kerr horizons r± = rg(1 ± √(1 − a*²)) (m)
    get rPlus() {
        return this.rg * (1 + Math.sqrt(1 - this.a ** 2));
    }

    get rMinus() {
        return this.rg * (1 - Math.sqrt(1 - this.a ** 2));
    }

    // Ergosphere (static limit) at the equator: r_E(π/2) = 2rg = rs (m)
    get rErgo() {
        return this.ergosphereRadius(Math.PI / 2);
    }

    // Static limit r_E(θ) = rg(1 + √(1 − a*²cos²θ)) (m), Boyer–Lindquist
    ergosphereRadius(theta) {
        return this.rg * (1 + Math.sqrt(1 - this.a ** 2 * Math.cos(theta) ** 2));
    }

    // Hawking temperature T_H = ħc(r₊ − r₋) / (4π k_B (r₊² + aLen²)) (K)
    // Equivalent to ħκ/(2πck_B); for a* = 0 it reduces to ħc³/(8πGMk_B).
    get hawkingTemperature() {
        return (CONSTANTS.hbar * CONSTANTS.c * (this.rPlus - this.rMinus)) /
            (4 * Math.PI * CONSTANTS.kB * (this.rPlus ** 2 + this.aLen ** 2));
    }

    // Bekenstein–Hawking entropy S/k_B = A/(4 l_P²), Kerr area A = 4π(r₊² + aLen²)
    get entropy() {
        const A = 4 * Math.PI * (this.rPlus ** 2 + this.aLen ** 2);
        return A / (4 * CONSTANTS.lP ** 2);
    }

    // Surface gravity κ = c²(r₊ − r₋) / (2(r₊² + aLen²)) (m/s²); a* = 0 → c⁴/(4GM)
    get surfaceGravity() {
        return (CONSTANTS.c ** 2 * (this.rPlus - this.rMinus)) /
            (2 * (this.rPlus ** 2 + this.aLen ** 2));
    }

    // MSS Lyapunov bound λ_L ≤ 2πk_B T/ħ (s⁻¹), Maldacena–Shenker–Stanford 2016
    get lyapunovBound() {
        return (2 * Math.PI * CONSTANTS.kB * this.hawkingTemperature) / CONSTANTS.hbar;
    }

    // Scrambling time t* = (ħ/(2πk_B T)) ln S (s)
    get scramblingTime() {
        const S = this.entropy;
        const T = this.hawkingTemperature;
        return (CONSTANTS.hbar / (2 * Math.PI * CONSTANTS.kB * T)) * Math.log(S);
    }

    // Page time of the naive model below: S_rad,thermal = S_BH at x = 1 − 2^(−3/2) ≈ 0.646 (s).
    // Naive model (no greybody factors, radiation entropy = BH entropy lost). Page (2013,
    // JCAP 09, 028; arXiv:1301.4995) includes the entropy produced by the emission
    // (β ≈ 1.4847) and finds t_Page ≈ 0.5381 t_evap.
    get pageTime() {
        return (1 - Math.pow(2, -1.5)) * this.evaporationTime;
    }

    // Evaporation time t_evap = 5120πG²M³/(ħc⁴) (s)
    // Schwarzschild, photon-only geometric-optics estimate: no greybody factors, no spin.
    get evaporationTime() {
        return (5120 * Math.PI * CONSTANTS.G ** 2 * this.M ** 3) /
            (CONSTANTS.hbar * CONSTANTS.c ** 4);
    }

    // Innermost Stable Circular Orbit (prograde), Bardeen–Press–Teukolsky 1972 (m)
    get isco() {
        if (this.a === 0) return 6 * this.rg;

        const Z1 = 1 + Math.cbrt(1 - this.a ** 2) * (Math.cbrt(1 + this.a) + Math.cbrt(1 - this.a));
        const Z2 = Math.sqrt(3 * this.a ** 2 + Z1 ** 2);

        return this.rg * (3 + Z2 - Math.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)));
    }

    // Prograde equatorial photon orbit r_ph = 2rg(1 + cos[(2/3)arccos(−a*)]) (m), BPT 1972
    get photonSphere() {
        return this.rs * (1 + Math.cos((2 / 3) * Math.acos(-this.a)));
    }

    // Maximum extractable rotational energy fraction 1 − M_irr/M,
    // M_irr/M = √((1 + √(1 − a*²))/2) (Christodoulou 1970). ≤ 1 − 1/√2 ≈ 29% at a* = 1.
    // Not the efficiency of a single Penrose event (max (√2 − 1)/2 ≈ 20.7%, Chandrasekhar 1983).
    get maxExtractableRotationalFraction() {
        return 1 - Math.sqrt(0.5 * (1 + Math.sqrt(1 - this.a ** 2)));
    }

    // Horizon angular velocity Ω_H = aLen·c/(r₊² + aLen²) = a*c/(2r₊) (rad/s)
    get horizonAngularVelocity() {
        return (this.a * CONSTANTS.c) / (2 * this.rPlus);
    }

    // OTOC F(t) = max(0, 1 − ε e^{λ_L t}) with ε = 1/S, written as exp(λ_L t − ln S)
    // so the exponent never overflows before subtraction. Dimensionless, deterministic.
    getOTOCDecay(time) {
        const exponent = this.lyapunovBound * time - Math.log(this.entropy);
        return Math.max(0, 1 - Math.exp(exponent));
    }

    // Naive Page curve components at x = t/t_evap (no greybody factor):
    // M(t) = M₀(1 − x)^{1/3}, S_BH ∝ M² → S_BH = S₀(1 − x)^{2/3},
    // S_rad,thermal = S₀ − S_BH, S_Page = min(S_rad,thermal, S_BH).
    getPageCurveComponents(t, tTotal) {
        const x = Math.min(Math.max(t / tTotal, 0), 1);
        const S0 = this.entropy;
        const sBH = S0 * Math.pow(1 - x, 2 / 3);
        const sRadThermal = S0 - sBH;
        return { x, sBH, sRadThermal, sPage: Math.min(sRadThermal, sBH) };
    }

    // Page curve entropy S_Page(t) (units of k_B)
    getPageCurveEntropy(t, tTotal) {
        return this.getPageCurveComponents(t, tTotal).sPage;
    }

    // Schwarzschild QNM, eikonal approximation (Cardoso et al. 2009, arXiv:0812.1806):
    // ω = Ω_c(l + 1/2) − i(n + 1/2)|λ_geo|, Ω_c = |λ_geo| = c³/(3√3 GM) (rad/s).
    // Convention e^{−iωt}, Im ω < 0. Error ~29% in Re ω at l = 2 (see SCHWARZSCHILD_QNM_REFERENCE).
    getQNMFrequency(l, n) {
        const omegaC = CONSTANTS.c / (3 * Math.sqrt(3) * this.rg);
        const omega_R = (l + 0.5) * omegaC;
        const omega_I = -(n + 0.5) * omegaC;
        return { real: omega_R, imag: omega_I };
    }

    // Lloyd's bound dC/dt ≤ 2Mc²/(πħ) (s⁻¹)
    get complexityGrowthRate() {
        return (2 * this.M * CONSTANTS.c ** 2) / (Math.PI * CONSTANTS.hbar);
    }

    // Format for display
    toDisplayObject() {
        return {
            mass: {
                solar: this.M_solar,
                kg: this.M,
                planck: this.M / CONSTANTS.mP
            },
            spin: this.a,
            radii: {
                schwarzschild: this.rs,
                eventHorizon: this.rPlus,
                cauchyHorizon: this.rMinus,
                ergosphere: this.rErgo,
                photonSphere: this.photonSphere,
                isco: this.isco
            },
            thermodynamics: {
                temperature: this.hawkingTemperature,
                entropy: this.entropy,
                surfaceGravity: this.surfaceGravity,
                evaporationTime: this.evaporationTime,
                evaporationTimeYears: this.evaporationTime / CONSTANTS.year
            },
            quantumChaos: {
                lyapunovBound: this.lyapunovBound,
                scramblingTime: this.scramblingTime,
                pageTime: this.pageTime,
                complexityRate: this.complexityGrowthRate
            },
            kerrProperties: {
                horizonAngularVelocity: this.horizonAngularVelocity,
                maxExtractableRotationalPercent: this.maxExtractableRotationalFraction * 100
            }
        };
    }
}

// Schwarzschild gravitational QNMs, Mω in G = c = 1, convention e^{−iωt} (Leaver 1985,
// Proc. R. Soc. A 402, 285; tabulated in Kokkotas & Schmidt 1999, Living Rev. Rel. 2, 2, Table 1).
BlackHolePhysics.SCHWARZSCHILD_QNM_REFERENCE = {
    2: [
        { n: 0, real: 0.37367, imag: -0.08896 },
        { n: 1, real: 0.34671, imag: -0.27391 },
        { n: 2, real: 0.30105, imag: -0.47828 }
    ],
    3: [{ n: 0, real: 0.59944, imag: -0.09270 }],
    4: [{ n: 0, real: 0.80918, imag: -0.09416 }]
};

// === Quantum Chaos Analyzer ===
class QuantumChaosAnalyzer {
    constructor() {
        this.SYK_N = 20; // Number of Majorana fermions
        this.SYK_J = 1;  // Coupling strength
    }

    // Wigner surmise for GUE level spacing P(s) = (32/π²) s² e^{−4s²/π}
    getGUESpacing(s) {
        return (32 / (Math.PI ** 2)) * s ** 2 * Math.exp(-(4 / Math.PI) * s ** 2);
    }

    // OTOC for SYK-like systems with maximal λ = 2π/β: F = max(0, 1 − e^{λt − ln N})
    getSYKOTOC(t, beta) {
        const lambda = 2 * Math.PI / beta; // MSS saturating
        return Math.max(0, 1 - Math.exp(lambda * t - Math.log(this.SYK_N)));
    }
}

// === Gravitational Wave Signal Generator ===
class GravitationalWaveSimulator {
    constructor(m1, m2, distance) {
        this.m1 = m1 * CONSTANTS.M_sun;
        this.m2 = m2 * CONSTANTS.M_sun;
        this.M_total = this.m1 + this.m2;
        // Chirp mass M_c = (m₁m₂)^{3/5} / (m₁ + m₂)^{1/5} (kg)
        this.M_chirp = Math.pow(this.m1 * this.m2, 3 / 5) / Math.pow(this.M_total, 1 / 5);
        this.distance = distance * CONSTANTS.ly;
    }

    // Leading-order (Newtonian) GW frequency at time-to-merger τ (Hz):
    // f_GW(τ) = (1/π)·(5/(256τ))^{3/8}·(GM_c/c³)^{−5/8}; orbital frequency f_orb = f_GW/2.
    getGWFrequency(t, t_merge) {
        const tau = t_merge - t;
        if (tau <= 0) return Infinity;

        const tChirp = CONSTANTS.G * this.M_chirp / CONSTANTS.c ** 3; // GM_c/c³ (s)
        return Math.pow(5 / (256 * tau), 3 / 8) * Math.pow(tChirp, -5 / 8) / Math.PI;
    }

    // Strain amplitude h = 4(GM_c)^{5/3}(πf_GW)^{2/3} / (c⁴D) (dimensionless)
    getStrainAmplitude(f) {
        const h = (4 * CONSTANTS.G ** (5 / 3) * this.M_chirp ** (5 / 3) *
            (Math.PI * f) ** (2 / 3)) / (CONSTANTS.c ** 4 * this.distance);
        return h;
    }

    // Generate waveform data; GW phase Φ += 2π f_GW dt
    generateWaveform(duration, sampleRate) {
        const samples = [];
        const dt = 1 / sampleRate;
        let phase = 0;

        for (let t = 0; t < duration; t += dt) {
            const tau = duration - t;
            if (tau < 0.001) break;

            const f = this.getGWFrequency(t, duration);
            const h = this.getStrainAmplitude(f);

            phase += 2 * Math.PI * f * dt;
            samples.push({
                time: t,
                strain: h * Math.cos(phase),
                frequency: f
            });
        }

        return samples;
    }
}

// === Visualization Data Generator ===
class VisualizationDataGenerator {
    // Generate (naive) Page curve data; phase switches where S_rad,thermal crosses S_BH
    static generatePageCurve(bh, points = 100) {
        const data = [];
        const tTotal = bh.evaporationTime;

        for (let i = 0; i <= points; i++) {
            const t = (i / points) * tTotal;
            const { sBH, sRadThermal, sPage } = bh.getPageCurveComponents(t, tTotal);
            data.push({
                time: t / CONSTANTS.year,
                entropy: sPage,
                sBH,
                sRadThermal,
                phase: sRadThermal < sBH ? "thermal" : "information"
            });
        }
        return data;
    }

    // Generate OTOC decay data
    static generateOTOCData(bh, points = 100) {
        const data = [];
        const tMax = 5 * bh.scramblingTime;

        for (let i = 0; i <= points; i++) {
            const t = (i / points) * tMax;
            const F = bh.getOTOCDecay(t);
            data.push({
                time: t,
                otoc: F,
                scrambled: t > bh.scramblingTime
            });
        }
        return data;
    }

    // Generate QNM spectrum (eikonal approximation, l ≥ 2 for gravitational modes)
    static generateQNMSpectrum(bh, lMax = 4, nMax = 5) {
        const modes = [];

        for (let l = 2; l <= lMax; l++) {
            for (let n = 0; n <= nMax; n++) {
                const omega = bh.getQNMFrequency(l, n);
                modes.push({
                    l,
                    n,
                    approximation: "eikonal",
                    omega_real: omega.real,
                    omega_imag: omega.imag,
                    frequency: omega.real / (2 * Math.PI),
                    damping_time: -1 / omega.imag
                });
            }
        }
        return modes;
    }

    // Kerr structure in the meridional plane using oblate (Kerr–Schild) coordinates:
    // x = √(r² + aLen²) sinθ, z = r cosθ, so the ring singularity (r = 0, θ = π/2)
    // sits at cylindrical radius aLen, inside the inner horizon.
    static generateKerrStructure(bh, thetaPoints = 50) {
        const structure = {
            eventHorizon: [],
            ergosphere: [],
            cauchyHorizon: [],
            singularity: []
        };

        const aLen = bh.aLen;
        const toOblate = (r, theta) => ({
            r,
            theta,
            x: Math.sqrt(r ** 2 + aLen ** 2) * Math.sin(theta),
            z: r * Math.cos(theta)
        });

        for (let i = 0; i <= thetaPoints; i++) {
            const theta = (i / thetaPoints) * Math.PI;
            structure.eventHorizon.push(toOblate(bh.rPlus, theta));
            structure.ergosphere.push(toOblate(bh.ergosphereRadius(theta), theta));
            structure.cauchyHorizon.push(toOblate(bh.rMinus, theta));
        }

        // Ring singularity of radius aLen in the equatorial plane
        for (let i = 0; i <= 50; i++) {
            const phi = (i / 50) * 2 * Math.PI;
            structure.singularity.push({
                x: aLen * Math.cos(phi),
                y: aLen * Math.sin(phi),
                z: 0
            });
        }

        return structure;
    }
}

// === Formatting Utilities (DOM-free) ===
function formatMass(mass) {
    if (mass >= 1e9) return (mass / 1e9).toFixed(1) + '×10⁹';
    if (mass >= 1e6) return (mass / 1e6).toFixed(1) + '×10⁶';
    if (mass >= 1000) return (mass / 1000).toFixed(1) + '×10³';
    return mass.toFixed(1);
}

function formatDistance(meters) {
    if (meters >= CONSTANTS.ly) return (meters / CONSTANTS.ly).toFixed(2) + ' ly';
    if (meters >= 1e12) return (meters / 1e12).toFixed(2) + ' Tm';
    if (meters >= 1e9) return (meters / 1e9).toFixed(2) + ' Gm';
    if (meters >= 1e6) return (meters / 1e6).toFixed(2) + ' Mm';
    if (meters >= 1e3) return (meters / 1e3).toFixed(2) + ' km';
    return meters.toFixed(2) + ' m';
}

function formatTemperature(kelvin) {
    if (kelvin < 1e-12) return formatScientific(kelvin) + ' K';
    if (kelvin < 1e-9) return (kelvin * 1e12).toFixed(2) + ' pK';
    if (kelvin < 1e-6) return (kelvin * 1e9).toFixed(2) + ' nK';
    if (kelvin < 1e-3) return (kelvin * 1e6).toFixed(2) + ' μK';
    return kelvin.toFixed(6) + ' K';
}

function formatTime(seconds) {
    const years = seconds / CONSTANTS.year;
    if (years >= 1e15) return formatScientific(years) + ' years';
    if (years >= 1e12) return (years / 1e12).toFixed(1) + ' trillion years';
    if (years >= 1e9) return (years / 1e9).toFixed(1) + ' billion years';
    if (years >= 1) return years.toFixed(1) + ' years';
    if (seconds >= 86400) return (seconds / 86400).toFixed(1) + ' days';
    if (seconds >= 3600) return (seconds / 3600).toFixed(1) + ' hours';
    if (seconds >= 1) return seconds.toFixed(2) + ' s';
    if (seconds >= 1e-3) return (seconds * 1e3).toFixed(2) + ' ms';
    if (seconds >= 1e-6) return (seconds * 1e6).toFixed(2) + ' μs';
    if (seconds >= 1e-9) return (seconds * 1e9).toFixed(2) + ' ns';
    return formatScientific(seconds) + ' s';
}

function formatScientific(num) {
    if (num === 0) return '0';
    const exp = Math.floor(Math.log10(Math.abs(num)));
    const mantissa = num / Math.pow(10, exp);

    if (exp >= -3 && exp <= 3) return num.toPrecision(4);
    return mantissa.toFixed(2) + '×10' + superscript(exp);
}

function superscript(n) {
    const chars = { '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', '-': '⁻' };
    return String(n).split('').map(c => chars[c] || c).join('');
}

const FORMATTERS = { formatMass, formatDistance, formatTemperature, formatTime, formatScientific, superscript };

if (typeof module !== 'undefined' && module.exports) module.exports = { CONSTANTS, BlackHolePhysics, GravitationalWaveSimulator, QuantumChaosAnalyzer, VisualizationDataGenerator, KNOWN_BLACK_HOLES, ...FORMATTERS };
else Object.assign(window, { CONSTANTS, BlackHolePhysics, GravitationalWaveSimulator, QuantumChaosAnalyzer, VisualizationDataGenerator, KNOWN_BLACK_HOLES, ...FORMATTERS });
