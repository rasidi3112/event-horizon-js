// === Pure physics engine (no DOM) ===

// === Physical Constants (SI Units) ===
const CONSTANTS = {
    G: 6.67430e-11,          // Gravitational constant (m³/kg/s²)
    c: 299792458,            // Speed of light (m/s)
    hbar: 1.054571817e-34,   // Reduced Planck constant (J·s)
    kB: 1.380649e-23,        // Boltzmann constant (J/K)
    lP: 1.616255e-35,        // Planck length (m)
    tP: 5.391247e-44,        // Planck time (s)
    mP: 2.176434e-8,         // Planck mass (kg)
    M_sun: 1.98892e30,       // Solar mass (kg)
    year: 3.15576e7,         // Seconds per year
    ly: 9.461e15,            // Light year in meters
};

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
        this.a = Math.min(Math.max(spin, 0), 0.998); // Kerr spin parameter (0 to 0.998)
    }

    // Schwarzschild radius
    get rs() {
        return (2 * CONSTANTS.G * this.M) / (CONSTANTS.c ** 2);
    }

    // Kerr horizons
    get rPlus() {
        const rg = this.rs / 2;
        return rg * (1 + Math.sqrt(1 - this.a ** 2));
    }

    get rMinus() {
        const rg = this.rs / 2;
        return rg * (1 - Math.sqrt(1 - this.a ** 2));
    }

    // Ergosphere outer boundary (at equator)
    get rErgo() {
        const rg = this.rs / 2;
        return rg * (1 + Math.sqrt(1 - this.a ** 2));
    }

    // Hawking temperature
    get hawkingTemperature() {
        if (this.a === 0) {
            return (CONSTANTS.hbar * CONSTANTS.c ** 3) / (8 * Math.PI * CONSTANTS.G * this.M * CONSTANTS.kB);
        }
        // Kerr black hole temperature
        const rPlus = this.rPlus;
        const rMinus = this.rMinus;
        return (CONSTANTS.hbar * CONSTANTS.c * (rPlus - rMinus)) /
            (4 * Math.PI * CONSTANTS.kB * (rPlus ** 2 + (this.a * this.rs / 2) ** 2));
    }

    // Bekenstein-Hawking entropy (in units of kB)
    get entropy() {
        const A = 4 * Math.PI * this.rPlus ** 2;
        return A / (4 * CONSTANTS.lP ** 2);
    }

    // Surface gravity
    get surfaceGravity() {
        const rPlus = this.rPlus;
        const rg = this.rs / 2;
        return (rPlus - rg) / (2 * (rPlus ** 2 + (this.a * rg) ** 2));
    }

    // MSS Lyapunov bound (maximum chaos rate)
    get lyapunovBound() {
        return (2 * Math.PI * CONSTANTS.kB * this.hawkingTemperature) / CONSTANTS.hbar;
    }

    // Scrambling time
    get scramblingTime() {
        const S = this.entropy;
        const T = this.hawkingTemperature;
        return (CONSTANTS.hbar / (2 * Math.PI * CONSTANTS.kB * T)) * Math.log(S);
    }

    // Page time (when entropy starts decreasing)
    get pageTime() {
        return this.evaporationTime / 2;
    }

    // Evaporation time
    get evaporationTime() {
        return (5120 * Math.PI * CONSTANTS.G ** 2 * this.M ** 3) /
            (CONSTANTS.hbar * CONSTANTS.c ** 4);
    }

    // Innermost Stable Circular Orbit (ISCO) for Kerr
    get isco() {
        if (this.a === 0) return 3 * this.rs;

        const Z1 = 1 + Math.cbrt(1 - this.a ** 2) * (Math.cbrt(1 + this.a) + Math.cbrt(1 - this.a));
        const Z2 = Math.sqrt(3 * this.a ** 2 + Z1 ** 2);
        const rg = this.rs / 2;

        // Prograde orbit
        return rg * (3 + Z2 - Math.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)));
    }

    // Photon sphere radius
    get photonSphere() {
        if (this.a === 0) return 1.5 * this.rs;
        // Simplified for equatorial prograde orbits
        return this.rs * (1 + Math.cos((2 / 3) * Math.acos(-this.a)));
    }

    // Penrose process maximum efficiency
    get penroseEfficiency() {
        return 1 - Math.sqrt(0.5 * (1 + Math.sqrt(1 - this.a ** 2)));
    }

    // Angular velocity at horizon (frame dragging)
    get horizonAngularVelocity() {
        return (this.a * CONSTANTS.c) / (2 * this.rPlus);
    }

    // Calculate OTOC decay rate
    getOTOCDecay(time) {
        const lambda = this.lyapunovBound;
        const tScramble = this.scramblingTime;
        const epsilon = 1e-4; // Initial perturbation

        if (time < tScramble) {
            return 1 - epsilon * Math.exp(lambda * time);
        }
        return 0.5 * (1 + Math.random() * 0.1); // Saturated regime with fluctuations
    }

    // Page curve entropy evolution
    getPageCurveEntropy(t, tTotal) {
        const ratio = t / tTotal;
        const Smax = this.entropy;
        const pageRatio = 0.5;

        if (ratio < pageRatio) {
            // Rising phase (thermal)
            return Smax * (ratio / pageRatio);
        } else {
            // Falling phase (information recovery)
            return Smax * (1 - (ratio - pageRatio) / (1 - pageRatio));
        }
    }

    // Quasinormal mode frequencies (Schwarzschild approximation)
    getQNMFrequency(l, n) {
        const omega_R = ((l + 0.5) * CONSTANTS.c) / (3 * Math.sqrt(3) * this.rs / 2);
        const omega_I = -(n + 0.5) * CONSTANTS.c / (3 * Math.sqrt(3) * this.rs / 2);
        return { real: omega_R, imag: omega_I };
    }

    // Complexity growth rate (Lloyd's bound)
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
                penroseEfficiency: this.penroseEfficiency * 100
            }
        };
    }
}

// === Quantum Chaos Analyzer ===
class QuantumChaosAnalyzer {
    constructor() {
        this.SYK_N = 20; // Number of Majorana fermions
        this.SYK_J = 1;  // Coupling strength
    }

    // Random Matrix Theory level spacing distribution (GUE)
    getGUESpacing(s) {
        return (32 / (Math.PI ** 2)) * s ** 2 * Math.exp(-(4 / Math.PI) * s ** 2);
    }

    // OTOC for SYK-like systems
    getSYKOTOC(t, beta) {
        const lambda = 2 * Math.PI / beta; // MSS saturating
        const tScramble = (beta / (2 * Math.PI)) * Math.log(this.SYK_N);

        if (t < tScramble) {
            return 1 - (1 / this.SYK_N) * Math.exp(lambda * t);
        }
        return 0.5 * (1 + 0.1 * Math.sin(t / beta));
    }
}

// === Gravitational Wave Signal Generator ===
class GravitationalWaveSimulator {
    constructor(m1, m2, distance) {
        this.m1 = m1 * CONSTANTS.M_sun;
        this.m2 = m2 * CONSTANTS.M_sun;
        this.M_total = this.m1 + this.m2;
        this.mu = (this.m1 * this.m2) / this.M_total;
        this.M_chirp = Math.pow(this.m1 * this.m2, 3 / 5) / Math.pow(this.M_total, 1 / 5);
        this.distance = distance * CONSTANTS.ly;
    }

    // Orbital frequency evolution
    getOrbitalFrequency(t, t_merge) {
        const tau = t_merge - t;
        if (tau <= 0) return Infinity;

        const f0 = Math.pow(
            (256 * Math.PI ** (8 / 3) * CONSTANTS.G ** (5 / 3) * this.M_chirp ** (5 / 3)) /
            (5 * CONSTANTS.c ** 5 * tau),
            3 / 8
        ) / (2 * Math.PI);

        return f0;
    }

    // Strain amplitude
    getStrainAmplitude(f) {
        const h = (4 * CONSTANTS.G ** (5 / 3) * this.M_chirp ** (5 / 3) *
            (Math.PI * f) ** (2 / 3)) / (CONSTANTS.c ** 4 * this.distance);
        return h;
    }

    // Generate waveform data
    generateWaveform(duration, sampleRate) {
        const samples = [];
        const dt = 1 / sampleRate;
        let phase = 0;

        for (let t = 0; t < duration; t += dt) {
            const tau = duration - t;
            if (tau < 0.001) break;

            const f = this.getOrbitalFrequency(t, duration);
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
    // Generate Page curve data
    static generatePageCurve(bh, points = 100) {
        const data = [];
        const tTotal = bh.evaporationTime;

        for (let i = 0; i <= points; i++) {
            const t = (i / points) * tTotal;
            const S = bh.getPageCurveEntropy(t, tTotal);
            data.push({
                time: t / CONSTANTS.year,
                entropy: S,
                phase: t < tTotal / 2 ? "thermal" : "information"
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

    // Generate QNM spectrum
    static generateQNMSpectrum(bh, lMax = 4, nMax = 5) {
        const modes = [];

        for (let l = 2; l <= lMax; l++) {
            for (let n = 0; n <= nMax; n++) {
                const omega = bh.getQNMFrequency(l, n);
                modes.push({
                    l,
                    n,
                    omega_real: omega.real,
                    omega_imag: omega.imag,
                    frequency: omega.real / (2 * Math.PI),
                    damping_time: -1 / omega.imag
                });
            }
        }
        return modes;
    }

    // Generate Kerr structure visualization data
    static generateKerrStructure(bh, thetaPoints = 50) {
        const structure = {
            eventHorizon: [],
            ergosphere: [],
            cauchyHorizon: [],
            singularity: []
        };

        const rg = bh.rs / 2;
        const a = bh.a;

        for (let i = 0; i <= thetaPoints; i++) {
            const theta = (i / thetaPoints) * Math.PI;

            // Event horizon (r+)
            const rPlus = rg * (1 + Math.sqrt(1 - a ** 2));
            structure.eventHorizon.push({
                r: rPlus,
                theta,
                x: rPlus * Math.sin(theta),
                z: rPlus * Math.cos(theta)
            });

            // Ergosphere
            const rErgo = rg * (1 + Math.sqrt(1 - a ** 2 * Math.cos(theta) ** 2));
            structure.ergosphere.push({
                r: rErgo,
                theta,
                x: rErgo * Math.sin(theta),
                z: rErgo * Math.cos(theta)
            });

            // Cauchy horizon (r-)
            const rMinus = rg * (1 - Math.sqrt(1 - a ** 2));
            structure.cauchyHorizon.push({
                r: rMinus,
                theta,
                x: rMinus * Math.sin(theta),
                z: rMinus * Math.cos(theta)
            });
        }

        // Ring singularity
        for (let i = 0; i <= 50; i++) {
            const phi = (i / 50) * 2 * Math.PI;
            structure.singularity.push({
                x: a * rg * Math.cos(phi),
                y: a * rg * Math.sin(phi),
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
    if (kelvin < 1e-9) return (kelvin * 1e12).toFixed(2) + ' pK';
    if (kelvin < 1e-6) return (kelvin * 1e9).toFixed(2) + ' nK';
    if (kelvin < 1e-3) return (kelvin * 1e6).toFixed(2) + ' μK';
    return kelvin.toFixed(6) + ' K';
}

function formatTime(seconds) {
    if (seconds >= 1e67) return formatScientific(seconds / CONSTANTS.year) + ' years';
    if (seconds >= CONSTANTS.year * 1e12) return (seconds / (CONSTANTS.year * 1e12)).toFixed(1) + ' trillion years';
    if (seconds >= CONSTANTS.year * 1e9) return (seconds / (CONSTANTS.year * 1e9)).toFixed(1) + ' billion years';
    if (seconds >= CONSTANTS.year) return (seconds / CONSTANTS.year).toFixed(1) + ' years';
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
