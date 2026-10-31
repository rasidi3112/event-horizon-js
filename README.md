<div align="center">
  <img src="hero-background.png" alt="Black Hole Explorer" width="100%" style="border-radius: 8px;">

  <h1>Black Hole Explorer & Quantum Chaos</h1>
  
  <p>
    <strong>A Computational Framework for Relativistic Astrophysics and Quantum Chaos Visualization</strong>
  </p>

  <p>
    <a href="https://api.nasa.gov/">
      <img src="https://img.shields.io/badge/Data-NASA_API-000000?style=flat-square" alt="NASA API">
    </a>
    <a href="LICENSE">
      <img src="https://img.shields.io/badge/License-MIT-000000?style=flat-square" alt="License">
    </a>
  </p>
</div>

<hr />

## Abstract

Black Hole Explorer is a web-based research interface designed to bridge the gap between theoretical physics and interactive data visualization. This platform integrates real-time astrophysical data from NASA's archives with client-side simulations of relativistic phenomena. The project focuses on visualizing Schwarzschild and Kerr horizon structure, Hawking radiation thermodynamics, and the diagnosis of quantum chaos through Out-of-Time-Order Correlators (OTOC) and Lyapunov exponents.

By leveraging modern web technologies (Vanilla JavaScript, CSS3), this repository demonstrates how complex mathematical models typically reserved for Python or MATLAB environments can be effectively rendered and explored in a browser environment with high performance.

## Core Architecture

### 1. Relativistic Physics Engine
The core application utilizes a custom JavaScript physics engine (`BlackHolePhysics` class in `physics.js`) to compute real-time properties of black holes based on observational parameters.
*   **Schwarzschild Radius Calculation:** $R_s = \frac{2GM}{c^2}$
*   **Hawking Temperature:** $T_H = \frac{\hbar \kappa}{2\pi c k_B}$, which reduces to $\frac{\hbar c^3}{8\pi G M k_B}$ for $a_* = 0$
*   **Bekenstein-Hawking Entropy:** $S_{BH} = \frac{k_B c^3 A}{4G\hbar}$ with the Kerr area $A = 4\pi(r_+^2 + a^2)$

### 2. Quantum Chaos Diagnostics
The simulation module implements numerical methods to visualize the butterfly effect in holographic systems.
*   **MSS Bound Visualization:** Demonstrates the saturation of the Lyapunov exponent bound $\lambda_L \le \frac{2\pi k_B T}{\hbar}$.
*   **OTOC Growth:** Simulator for $C(t) = -\langle [W(t), V(0)]^2 \rangle \sim e^{\lambda_L t}$.

### 3. Data Integration Layer
A robust fetch wrapper interfaces directly with the **NASA Image and Video Library API**, filtering for high-relevance astrophysical imagery (accretion disks, event horizons, active galactic nuclei) and presenting them in a curated, museum-grade gallery interface.

## Technical Stack

*   **Frontend:** HTML5, CSS3 (Custom Variables, Grid/Flexbox Layouts, Glassmorphism UI)
*   **Simulation Logic:** Vanilla JavaScript (ES6+ Classes, Async/Await)
*   **Visualization:** Chart.js 4.4.1 + chartjs-plugin-annotation 3.1.0 for thermodynamic curves and chaos bounds
*   **Math Rendering:** MathJax 3.2.2 for LaTeX equation support
*   **API:** NASA Open APIs

All CDN scripts are pinned to exact versions and loaded with Subresource Integrity hashes.

## Installation & Usage

This project is client-side only and requires no build pipeline or package manager.

1.  **Clone the repository**
    ```bash
    git clone https://github.com/rasidi3112/event-horizon-js.git
    ```

2.  **Navigate to directory**
    ```bash
    cd event-horizon-js
    ```

3.  **Launch**
    Open `index.html` in any modern web browser (Chrome, Firefox, Safari, Edge).
    
    *Note: For API functionality to work without strict CORS policies on some local environments, consider running a simple python server:*
    ```bash
    python3 -m http.server 8000
    ```

## Running the Tests

The physics engine has regression tests that use only the Node.js built-in test runner (Node 18+, no dependencies):

```bash
node --test tests/*.test.cjs
```

(`node --test` without arguments also discovers the file. On Node 21+ a bare directory argument such as `tests/` is treated as a glob pattern and does not match the test file.)

## Project Structure

```bash
.
├── index.html            # Main dashboard & core UI
├── advanced-math.html    # Mathematical derivations & formalisms
├── simulations.html      # Interactive chaos & physics simulators
├── styles.css            # Global design system & responsive layout
├── physics.js            # Pure physics engine, black hole catalog & formatters (no DOM)
├── app.js                # UI logic & NASA API integration
├── hero-background.png   # Hero image
├── promo-dashboard.png   # Promotional screenshots
├── promo-image.png
├── tests/
│   └── physics.test.cjs  # Physics regression tests (node --test)
└── LICENSE               # MIT License
```

## Known Approximations

*   **Quasinormal modes** are computed in the eikonal (geometric-optics) limit, $\omega = \Omega_c(\ell + \tfrac12) - i(n + \tfrac12)|\lambda_{\rm geo}|$. For the Schwarzschild $\ell = 2$ fundamental mode the real part is off by about 29% compared with Leaver's values, which are shown alongside for reference.
*   **Page curve** uses a naive model: $M(t) = M_0(1 - t/t_{\rm evap})^{1/3}$, $S_{\rm rad} = S_0 - S_{BH}$, no greybody factors and no entropy production in the emission. It peaks at $t = (1 - 2^{-3/2})\,t_{\rm evap} \approx 0.65\,t_{\rm evap}$; Page (2013) obtains $\approx 0.54\,t_{\rm evap}$ once the radiation entropy production is included.
*   **Evaporation time** $t_{\rm evap} = 5120\pi G^2 M^3/(\hbar c^4)$ ignores greybody factors, particle species and spin.
*   **Spins** flagged as *assumed* in the catalog (Sgr A*, M87*, TON 618, NGC 1277, Holmberg 15A*, Phoenix A*, GW150914 components) are placeholders, not measurements. The engine caps spin at the Thorne limit $a_* = 0.998$.
*   **OTOC** curves use $F(t) = \max(0, 1 - e^{\lambda_L t}/S)$, i.e. the early-time exponential growth with $\epsilon = 1/S$; the saturation regime is not modelled.
*   **Chaos-bound** and **complexity** charts are schematic/illustrative.

## References

*   J. Maldacena, S. Shenker, D. Stanford, *A bound on chaos*, JHEP 08 (2016) 106, [arXiv:1503.01409](https://arxiv.org/abs/1503.01409)
*   J. Maldacena, D. Stanford, *Comments on the Sachdev-Ye-Kitaev model*, Phys. Rev. D 94 (2016) 106002, [arXiv:1604.07818](https://arxiv.org/abs/1604.07818)
*   D. N. Page, *Time Dependence of Hawking Radiation Entropy*, JCAP 09 (2013) 028, [arXiv:1301.4995](https://arxiv.org/abs/1301.4995)
*   V. Cardoso et al., *Geodesic stability, Lyapunov exponents and quasinormal modes*, Phys. Rev. D 79 (2009) 064016, [arXiv:0812.1806](https://arxiv.org/abs/0812.1806)
*   K. D. Kokkotas, B. G. Schmidt, *Quasi-Normal Modes of Stars and Black Holes*, Living Rev. Rel. 2 (1999) 2, [arXiv:gr-qc/9909058](https://arxiv.org/abs/gr-qc/9909058) (Leaver's Schwarzschild QNM table)
*   D. Carmi et al., *On the Time Dependence of Holographic Complexity*, [arXiv:1709.10184](https://arxiv.org/abs/1709.10184)
*   C. Kehle, R. Unger, *Gravitational collapse to extremal black holes and the third law of black hole thermodynamics*, [arXiv:2211.15742](https://arxiv.org/abs/2211.15742)
*   J. C. A. Miller-Jones et al., *Cygnus X-1 contains a 21-solar mass black hole*, Science 371 (2021) 1046, [arXiv:2102.09091](https://arxiv.org/abs/2102.09091)
*   GRAVITY Collaboration, *A geometric distance measurement to the Galactic Center black hole with 0.3% uncertainty*, A&A 625 (2019) L10, [arXiv:1904.05721](https://arxiv.org/abs/1904.05721)
*   Event Horizon Telescope Collaboration, *First M87 Event Horizon Telescope Results. VI*, ApJL 875 (2019) L6, [arXiv:1906.11243](https://arxiv.org/abs/1906.11243)
*   LIGO Scientific & Virgo Collaborations, *GWTC-1*, Phys. Rev. X 9 (2019) 031040, [arXiv:1811.12907](https://arxiv.org/abs/1811.12907)
*   LIGO-Virgo-KAGRA Collaboration, *GWTC-4.0*, [arXiv:2508.18082](https://arxiv.org/abs/2508.18082)
*   A. W. Graham et al., *A normal supermassive black hole in NGC 1277*, [arXiv:1601.05151](https://arxiv.org/abs/1601.05151)
*   K. Mehrgan et al., *A 40-billion solar mass black hole in the extreme core of Holm 15A*, [arXiv:1907.10608](https://arxiv.org/abs/1907.10608)
*   M. Brockamp et al., *Unveiling Gargantua*, [arXiv:1509.04782](https://arxiv.org/abs/1509.04782)
*   M. J. Reid et al., *A Parallax Distance to the Microquasar GRS 1915+105*, [arXiv:1409.2453](https://arxiv.org/abs/1409.2453); J. E. McClintock et al., *The Spin of the Near-Extreme Kerr Black Hole GRS 1915+105*, [arXiv:astro-ph/0606076](https://arxiv.org/abs/astro-ph/0606076)

## Contributing

Contributions to the physics engine or UI refinements are welcome. Please ensure all mathematical implementations cite valid sources (arXiv, standard textbooks) in pull request descriptions, and that `node --test tests/*.test.cjs` passes.

1.  Fork the Project
2.  Create your Feature Branch (`git checkout -b feature/AmazingPhysics`)
3.  Commit your Changes (`git commit -m 'Add Kerr geodesic integrator'`)
4.  Push to the Branch (`git push origin feature/AmazingPhysics`)
5.  Open a Pull Request

## License

Distributed under the MIT License. See `LICENSE` for more information.

---

<div align="center">
  <small>&copy; 2026 Black Hole Explorer Research Platform. Data courtesy of NASA.</small>
</div>
