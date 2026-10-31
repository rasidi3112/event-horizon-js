// === Real-time ray-traced Schwarzschild black hole (WebGL 1, no dependencies) ===
// Units: Schwarzschild radius r_s = 2GM/c² = 1, so the photon sphere is at r = 1.5,
// the shadow (critical impact parameter) at b = √27/2 ≈ 2.598 and the ISCO at r = 3.
// Each pixel's photon is integrated backwards with
//     d²x/dλ² = −(3/2) h² x / r⁵,   h = |x × dx/dλ| (conserved),
// which reproduces Schwarzschild null geodesics in these Cartesian-like coordinates
// (orbit equation u'' + u = (3/2) r_s u², Binet form). The thin disk lies in y = 0
// between the ISCO and an outer radius; its colour uses an approximate blackbody ramp
// with Doppler beaming and gravitational redshift g = D·√(1 − r_s/r).
(function (global) {
    'use strict';

    const VERTEX_SHADER = `
attribute vec2 aPos;
void main() {
    gl_Position = vec4(aPos, 0.0, 1.0);
}`;

    const FRAGMENT_SHADER = `
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif

#define MAX_STEPS 320

uniform vec2 uRes;
uniform vec2 uCenter;
uniform float uTime;
uniform float uFov;
uniform vec3 uCamPos;
uniform mat3 uCamRot;
uniform int uSteps;
uniform float uDiskOuter;
uniform float uDiskIntensity;
uniform float uExposure;

const float R_IN = 3.0;      // ISCO for a Schwarzschild black hole (r_s = 1)
const float R_ESCAPE = 45.0;

float hash13(vec3 p) {
    p = fract(p * 0.1031);
    p += dot(p, p.zyx + 31.32);
    return fract((p.x + p.y) * p.z);
}

float noise3(vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
    f = f * f * (3.0 - 2.0 * f);
    return mix(
        mix(mix(hash13(i), hash13(i + vec3(1.0, 0.0, 0.0)), f.x),
            mix(hash13(i + vec3(0.0, 1.0, 0.0)), hash13(i + vec3(1.0, 1.0, 0.0)), f.x), f.y),
        mix(mix(hash13(i + vec3(0.0, 0.0, 1.0)), hash13(i + vec3(1.0, 0.0, 1.0)), f.x),
            mix(hash13(i + vec3(0.0, 1.0, 1.0)), hash13(i + vec3(1.0, 1.0, 1.0)), f.x), f.y),
        f.z);
}

float fbm(vec3 p) {
    float sum = 0.0;
    float amp = 0.5;
    for (int i = 0; i < 5; i++) {
        sum += amp * noise3(p);
        p = p * 2.03 + vec3(1.7, 9.2, 3.1);
        amp *= 0.5;
    }
    return sum;
}

vec3 starLayer(vec3 dir, float scale, float threshold) {
    vec3 q = dir * scale;
    vec3 id = floor(q);
    vec3 f = fract(q) - 0.5;
    float h = hash13(id);
    if (h < threshold) return vec3(0.0);
    vec3 jitter = vec3(hash13(id + 7.1), hash13(id + 3.7), hash13(id + 1.3)) - 0.5;
    float d = length(f - jitter * 0.6);
    float strength = (h - threshold) / (1.0 - threshold);
    float b = pow(max(0.0, 1.0 - d * 2.4), 14.0) * (0.4 + 1.6 * strength);
    vec3 tint = mix(vec3(0.70, 0.80, 1.00), vec3(1.00, 0.86, 0.72), hash13(id + 5.3));
    return tint * b;
}

vec3 background(vec3 dir) {
    vec3 col = vec3(0.0);
    float n1 = fbm(dir * 2.6 + vec3(3.0, 1.0, 7.0));
    float n2 = fbm(dir * 4.8 + vec3(11.0, 5.0, 2.0));
    col += vec3(0.020, 0.035, 0.080) * n1 * n1 * 2.2;
    col += vec3(0.090, 0.030, 0.050) * pow(n2, 3.0) * 1.6;
    col += starLayer(dir, 160.0, 0.965) * 1.4;
    col += starLayer(dir, 290.0, 0.975);
    col += starLayer(dir, 520.0, 0.982) * 0.7;
    return col;
}

// Approximate blackbody colour ramp, t ≈ 1 near the disk's peak temperature
vec3 temperatureRamp(float t) {
    t = clamp(t, 0.0, 1.6);
    vec3 deepRed = vec3(0.55, 0.10, 0.02);
    vec3 orange = vec3(1.00, 0.45, 0.10);
    vec3 warmWhite = vec3(1.00, 0.82, 0.55);
    vec3 blueWhite = vec3(0.72, 0.84, 1.00);
    if (t < 0.5) return mix(deepRed, orange, t / 0.5);
    if (t < 1.0) return mix(orange, warmWhite, (t - 0.5) / 0.5);
    return mix(warmWhite, blueWhite, (t - 1.0) / 0.6);
}

float diskPattern(float r, float angle, float seed) {
    vec3 q = vec3(cos(angle) * r, sin(angle) * r, seed);
    float clouds = fbm(q * 0.85);
    float streaks = noise3(vec3(r * 5.5, cos(angle) * 1.5 + seed, sin(angle) * 1.5));
    return mix(clouds, streaks, 0.35);
}

// Keplerian flow, Ω = √(M/r³) with M = 1/2; two phase-shifted samples avoid runaway shear
float animatedDisk(float r, float phi) {
    const float PERIOD = 16.0;
    float omega = sqrt(0.5 / (r * r * r)) * 4.0;
    float p1 = fract(uTime / PERIOD);
    float p2 = fract(uTime / PERIOD + 0.5);
    float w1 = 1.0 - abs(2.0 * p1 - 1.0);
    float a = diskPattern(r, phi - omega * p1 * PERIOD, 0.0);
    float b = diskPattern(r, phi - omega * p2 * PERIOD, 13.7);
    return mix(b, a, w1);
}

vec4 diskSample(vec3 hit, float r, vec3 rayVel) {
    float phi = atan(hit.z, hit.x);
    float pattern = animatedDisk(r, phi);

    // Novikov–Thorne-like temperature profile T ∝ r^(-3/4) (1 − √(r_in/r))^(1/4), peak ≈ 0.488
    float x = R_IN / r;
    float temp = pow(x, 0.75) * pow(max(1.0 - sqrt(x), 0.0), 0.25) / 0.488;

    // Orbital speed β = √(M/(r − 2M)) (static observer), prograde in +φ
    float beta = min(sqrt(0.5 / (r - 1.0)), 0.7);
    float gamma = 1.0 / sqrt(1.0 - beta * beta);
    vec3 orbitDir = normalize(vec3(-hit.z, 0.0, hit.x));
    float cosAngle = dot(orbitDir, -normalize(rayVel));
    float doppler = 1.0 / (gamma * (1.0 - beta * cosAngle));
    float g = doppler * sqrt(1.0 - 1.0 / r);

    float edge = smoothstep(R_IN, R_IN + 0.6, r) * (1.0 - smoothstep(uDiskOuter * 0.65, uDiskOuter, r));
    float density = clamp((0.30 + 0.95 * pattern) * edge, 0.0, 0.95);
    vec3 emission = temperatureRamp(temp * g * 1.05) * 2.4 * temp * temp * pow(g, 3.0) * (0.55 + pattern);
    return vec4(emission * uDiskIntensity, density);
}

vec3 aces(vec3 x) {
    return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
}

void main() {
    vec2 p = (gl_FragCoord.xy - 0.5 * uRes) / uRes.y - uCenter;
    vec3 vel = normalize(uCamRot * vec3(p * uFov, 1.0));
    vec3 pos = uCamPos;
    vec3 hv = cross(pos, vel);
    float h2 = dot(hv, hv);

    vec3 col = vec3(0.0);
    float alpha = 0.0;
    bool captured = false;

    for (int i = 0; i < MAX_STEPS; i++) {
        if (i >= uSteps) break;
        float r = length(pos);
        float dt = clamp(0.06 * r, 0.015, 2.5);

        // Velocity Verlet step of d²x/dλ² = −1.5 h² x / r⁵
        vec3 prev = pos;
        vec3 acc = -1.5 * h2 * pos / pow(r, 5.0);
        vel += acc * (0.5 * dt);
        pos += vel * dt;
        float rNew = length(pos);
        acc = -1.5 * h2 * pos / pow(rNew, 5.0);
        vel += acc * (0.5 * dt);

        if (prev.y * pos.y < 0.0) {
            float t = prev.y / (prev.y - pos.y);
            vec3 hit = mix(prev, pos, t);
            float rh = length(hit.xz);
            if (rh > R_IN && rh < uDiskOuter) {
                vec4 d = diskSample(hit, rh, vel);
                col += (1.0 - alpha) * d.rgb * d.a;
                alpha += (1.0 - alpha) * d.a;
            }
        }

        if (rNew < 1.0) { captured = true; break; }
        if (alpha > 0.985) break;
        if (rNew > R_ESCAPE && dot(pos, vel) > 0.0) break;
    }

    if (!captured) {
        col += (1.0 - alpha) * background(normalize(vel));
    }

    col *= uExposure;
    col = aces(col);
    vec2 q = gl_FragCoord.xy / uRes;
    col *= 0.85 + 0.15 * pow(16.0 * q.x * q.y * (1.0 - q.x) * (1.0 - q.y), 0.2);
    gl_FragColor = vec4(pow(col, vec3(1.0 / 2.2)), 1.0);
}`;

    const DEFAULTS = {
        distance: 22,          // camera radius (r_s)
        elevation: 8,          // degrees above the disk plane
        azimuth: 0,            // radians
        autoRotate: 0.035,     // rad/s orbit of the camera
        fov: 50,               // vertical field of view (degrees)
        offset: [0, 0],        // screen offset of the hole, units of viewport height (or function(w, h))
        diskOuter: 14,         // outer disk radius (r_s)
        diskIntensity: 1.0,
        exposure: 1.25,
        maxPixelRatio: 1.5,
        steps: null,           // ray-march steps (auto when null)
        interactive: true,
        onFirstFrame: null,
        onCameraChange: null
    };

    function compile(gl, type, source) {
        const shader = gl.createShader(type);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            const log = gl.getShaderInfoLog(shader);
            gl.deleteShader(shader);
            throw new Error(`Shader compile failed: ${log}`);
        }
        return shader;
    }

    function prefersReducedMotion() {
        return !!(global.matchMedia && global.matchMedia('(prefers-reduced-motion: reduce)').matches);
    }

    class BlackHole3D {
        static isSupported() {
            try {
                const canvas = document.createElement('canvas');
                return !!(global.WebGLRenderingContext &&
                    (canvas.getContext('webgl') || canvas.getContext('experimental-webgl')));
            } catch {
                return false;
            }
        }

        constructor(canvas, options = {}) {
            this.canvas = canvas;
            this.options = Object.assign({}, DEFAULTS, options);
            this.azimuth = this.options.azimuth;
            this.elevation = this.options.elevation * Math.PI / 180;
            this.time = 0;
            this.running = false;
            this.paused = false;
            this.visible = true;
            this.firstFrameDone = false;
            this.reducedMotion = prefersReducedMotion();
            this.isMobile = global.matchMedia ? global.matchMedia('(max-width: 768px)').matches : false;
            this.renderScale = this.isMobile ? 0.6 : 0.85;
            this.frameTimes = [];
            this.lastQualityChange = 0;
            this.rafId = 0;
            this.lastTimestamp = 0;

            this._frame = this._frame.bind(this);
            this._initGL();
            this._bindEvents();
            this.resize();
        }

        _initGL() {
            const attrs = { alpha: false, antialias: false, depth: false, powerPreference: 'high-performance' };
            const gl = this.canvas.getContext('webgl', attrs) || this.canvas.getContext('experimental-webgl', attrs);
            if (!gl) throw new Error('WebGL is not available');
            this.gl = gl;

            const program = gl.createProgram();
            gl.attachShader(program, compile(gl, gl.VERTEX_SHADER, VERTEX_SHADER));
            gl.attachShader(program, compile(gl, gl.FRAGMENT_SHADER, FRAGMENT_SHADER));
            gl.linkProgram(program);
            if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
                throw new Error(`Program link failed: ${gl.getProgramInfoLog(program)}`);
            }
            gl.useProgram(program);
            this.program = program;

            // Single full-screen triangle
            const buffer = gl.createBuffer();
            gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
            const loc = gl.getAttribLocation(program, 'aPos');
            gl.enableVertexAttribArray(loc);
            gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);

            this.uniforms = {};
            for (const name of ['uRes', 'uCenter', 'uTime', 'uFov', 'uCamPos', 'uCamRot', 'uSteps',
                'uDiskOuter', 'uDiskIntensity', 'uExposure']) {
                this.uniforms[name] = gl.getUniformLocation(program, name);
            }
        }

        _bindEvents() {
            if (global.ResizeObserver) {
                this.resizeObserver = new ResizeObserver(() => this.resize());
                this.resizeObserver.observe(this.canvas);
            } else {
                this._onResize = () => this.resize();
                global.addEventListener('resize', this._onResize);
            }

            if (global.IntersectionObserver) {
                this.intersectionObserver = new IntersectionObserver(entries => {
                    this.visible = entries.some(e => e.isIntersecting);
                    if (this.visible) this._requestFrame();
                });
                this.intersectionObserver.observe(this.canvas);
            }

            this._onVisibility = () => {
                if (!document.hidden) {
                    this.lastTimestamp = 0;
                    this._requestFrame();
                }
            };
            document.addEventListener('visibilitychange', this._onVisibility);

            this.canvas.addEventListener('webglcontextlost', e => {
                e.preventDefault();
                this.contextLost = true;
            });
            this.canvas.addEventListener('webglcontextrestored', () => {
                this.contextLost = false;
                this._initGL();
                this.resize();
                this._requestFrame();
            });

            if (this.options.interactive) {
                let dragging = false;
                let lastX = 0;
                let lastY = 0;
                this.canvas.addEventListener('pointerdown', e => {
                    dragging = true;
                    lastX = e.clientX;
                    lastY = e.clientY;
                    this.canvas.setPointerCapture?.(e.pointerId);
                });
                this.canvas.addEventListener('pointermove', e => {
                    if (!dragging) return;
                    const dx = e.clientX - lastX;
                    const dy = e.clientY - lastY;
                    lastX = e.clientX;
                    lastY = e.clientY;
                    this.azimuth -= dx * 0.006;
                    this.elevation = Math.max(-1.2, Math.min(1.2, this.elevation + dy * 0.004));
                    this.options.onCameraChange?.(this.elevation * 180 / Math.PI, this.azimuth);
                    this._requestFrame();
                });
                const stop = () => { dragging = false; };
                this.canvas.addEventListener('pointerup', stop);
                this.canvas.addEventListener('pointercancel', stop);
                this.canvas.addEventListener('pointerleave', stop);
            }
        }

        setOptions(options) {
            Object.assign(this.options, options);
            if ('elevation' in options) this.elevation = options.elevation * Math.PI / 180;
            if ('azimuth' in options) this.azimuth = options.azimuth;
            this._requestFrame();
        }

        setPaused(paused) {
            this.paused = paused;
            this.lastTimestamp = 0;
            this._requestFrame();
        }

        resize() {
            const rect = this.canvas.getBoundingClientRect();
            const dpr = Math.min(global.devicePixelRatio || 1, this.options.maxPixelRatio);
            const width = Math.max(1, Math.round(rect.width * dpr * this.renderScale));
            const height = Math.max(1, Math.round(rect.height * dpr * this.renderScale));
            if (this.canvas.width !== width || this.canvas.height !== height) {
                this.canvas.width = width;
                this.canvas.height = height;
            }
            this._requestFrame();
        }

        start() {
            this.running = true;
            this._requestFrame();
        }

        stop() {
            this.running = false;
            if (this.rafId) cancelAnimationFrame(this.rafId);
            this.rafId = 0;
        }

        destroy() {
            this.stop();
            this.resizeObserver?.disconnect();
            this.intersectionObserver?.disconnect();
            if (this._onResize) global.removeEventListener('resize', this._onResize);
            document.removeEventListener('visibilitychange', this._onVisibility);
        }

        _requestFrame() {
            if (!this.running || this.rafId) return;
            this.rafId = requestAnimationFrame(this._frame);
        }

        _animating() {
            return !this.paused && !this.reducedMotion && this.visible && !document.hidden;
        }

        _frame(timestamp) {
            this.rafId = 0;
            if (!this.running || this.contextLost) return;

            const dt = this.lastTimestamp ? Math.min((timestamp - this.lastTimestamp) / 1000, 0.1) : 0;
            this.lastTimestamp = timestamp;

            if (this._animating()) {
                this.time += dt;
                this.azimuth += dt * this.options.autoRotate;
                this._adaptQuality(dt, timestamp);
            }

            this._render();

            if (!this.firstFrameDone) {
                this.firstFrameDone = true;
                this.options.onFirstFrame?.();
            }
            if (this._animating()) this._requestFrame();
            else this.lastTimestamp = 0;
        }

        // Keep the frame time near 60 fps by scaling the internal resolution
        _adaptQuality(dt, timestamp) {
            if (dt <= 0) return;
            this.frameTimes.push(dt);
            if (this.frameTimes.length < 30) return;
            const avg = this.frameTimes.reduce((a, b) => a + b, 0) / this.frameTimes.length;
            this.frameTimes.length = 0;
            if (timestamp - this.lastQualityChange < 1500) return;

            let next = this.renderScale;
            if (avg > 0.024) next = Math.max(0.35, this.renderScale * 0.85);
            else if (avg < 0.0175) next = Math.min(1, this.renderScale * 1.08);
            if (Math.abs(next - this.renderScale) > 0.01) {
                this.renderScale = next;
                this.lastQualityChange = timestamp;
                this.resize();
            }
        }

        _render() {
            const gl = this.gl;
            const u = this.uniforms;
            const w = this.canvas.width;
            const h = this.canvas.height;
            const o = this.options;

            const ce = Math.cos(this.elevation);
            const cam = [
                o.distance * ce * Math.sin(this.azimuth),
                o.distance * Math.sin(this.elevation),
                o.distance * ce * Math.cos(this.azimuth)
            ];
            // Camera basis: forward toward the hole, right = forward × worldUp, up = right × forward
            const len = Math.hypot(cam[0], cam[1], cam[2]);
            const f = [-cam[0] / len, -cam[1] / len, -cam[2] / len];
            let r = [f[1] * 0 - f[2] * 1, f[2] * 0 - f[0] * 0, f[0] * 1 - f[1] * 0];
            const rl = Math.hypot(r[0], r[1], r[2]) || 1;
            r = [r[0] / rl, r[1] / rl, r[2] / rl];
            const up = [r[1] * f[2] - r[2] * f[1], r[2] * f[0] - r[0] * f[2], r[0] * f[1] - r[1] * f[0]];

            const offset = typeof o.offset === 'function' ? o.offset(w, h) : o.offset;
            const steps = o.steps || (this.isMobile ? 180 : 260);

            gl.viewport(0, 0, w, h);
            gl.uniform2f(u.uRes, w, h);
            gl.uniform2f(u.uCenter, offset[0], offset[1]);
            gl.uniform1f(u.uTime, this.time);
            gl.uniform1f(u.uFov, 2 * Math.tan((o.fov * Math.PI / 180) / 2));
            gl.uniform3f(u.uCamPos, cam[0], cam[1], cam[2]);
            gl.uniformMatrix3fv(u.uCamRot, false, new Float32Array([...r, ...up, ...f]));
            gl.uniform1i(u.uSteps, Math.min(steps, 320));
            gl.uniform1f(u.uDiskOuter, o.diskOuter);
            gl.uniform1f(u.uDiskIntensity, o.diskIntensity);
            gl.uniform1f(u.uExposure, o.exposure);
            gl.drawArrays(gl.TRIANGLES, 0, 3);
        }
    }

    global.BlackHole3D = BlackHole3D;
})(window);
