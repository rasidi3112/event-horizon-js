// Physics engine, catalog and formatters are loaded from physics.js (must be included first).

// === NASA API Configuration ===
const NASA_IMAGES = 'https://images-api.nasa.gov';

// === NASA Data Fetcher ===
class NASADataFetcher {
    constructor() {
        this.cache = new Map();
    }

    // Throws on network failure or non-2xx responses; callers report the error.
    async searchBlackHoleImages(query = "black hole", limit = 24) {
        const cacheKey = `images_${query}_${limit}`;
        if (this.cache.has(cacheKey)) {
            return this.cache.get(cacheKey);
        }

        const response = await fetch(
            `${NASA_IMAGES}/search?q=${encodeURIComponent(query)}&media_type=image&page_size=${limit}`
        );
        if (!response.ok) {
            throw new Error(`NASA API responded with HTTP ${response.status}`);
        }
        const data = await response.json();
        const results = data.collection.items.filter(item => item.links && item.links[0]);
        this.cache.set(cacheKey, results);
        return results;
    }

    // Get specific black hole related images
    async getBlackHoleObservations() {
        const queries = [
            "black hole EHT",
            "sagittarius A",
            "M87 black hole",
            "gravitational waves",
            "quasar",
            "active galactic nuclei"
        ];

        const allResults = [];
        for (const query of queries) {
            const results = await this.searchBlackHoleImages(query, 6);
            allResults.push(...results);
        }

        return allResults;
    }
}

// === Application State ===
const AppState = {
    selectedBlackHole: null,
    currentPhysics: null,
    chaosAnalyzer: new QuantumChaosAnalyzer(),
    nasaFetcher: new NASADataFetcher(),
    charts: {}
};

// === Initialize Application ===
document.addEventListener('DOMContentLoaded', () => {
    initNavigation();
    initEventListeners();
    initBlackHoleDatabase();
    updateAPIStatus('connecting');
    searchGallery('black hole');

    // Initialize with Sagittarius A*
    selectBlackHole(KNOWN_BLACK_HOLES[0]);

    console.log(`
╔═══════════════════════════════════════════════════════════════════════╗
║       BLACK HOLE EXPLORER - QUANTUM CHAOS RESEARCH ENGINE v2.0        ║
╠═══════════════════════════════════════════════════════════════════════╣
║                                                                       ║
║  PHYSICS ENGINE LOADED:                                               ║
║  • Kerr horizons, ergosphere, ISCO, surface gravity                   ║
║  • Hawking radiation & thermodynamics                                 ║
║  • Quantum chaos (OTOC, Lyapunov, scrambling)                        ║
║  • SYK model correlation functions                                    ║
║  • Quasinormal mode spectrum (eikonal)                                ║
║  • Gravitational wave templates                                       ║
║  • Page curve evolution (naive model)                                 ║
║                                                                       ║
║  KNOWN BLACK HOLES: ${KNOWN_BLACK_HOLES.length} objects loaded                                    ║
║  DATA SOURCE: NASA APIs, EHT, GRAVITY, LIGO/Virgo                    ║
║                                                                       ║
║  KEY FORMULAS:                                                        ║
║  • r_s = 2GM/c²              (Schwarzschild radius)                  ║
║  • T_H = ℏc³/(8πGMk_B)       (Hawking temperature)                   ║
║  • S_BH = kA/(4l_P²)         (Bekenstein-Hawking entropy)            ║
║  • λ_L ≤ 2πk_BT/ℏ           (MSS chaos bound)                       ║
║  • t_* ~ (β/2π)ln(S)         (Scrambling time)                       ║
║                                                                       ║
╚═══════════════════════════════════════════════════════════════════════╝
    `);
});

// === Black Hole Database UI ===
function initBlackHoleDatabase() {
    const container = document.getElementById('bhDatabase');
    if (!container) return;

    container.replaceChildren(...KNOWN_BLACK_HOLES.map((bh, index) => {
        const entry = document.createElement('div');
        entry.className = 'bh-entry';
        entry.dataset.index = String(index);
        entry.addEventListener('click', () => selectBlackHole(bh));

        const name = document.createElement('div');
        name.className = 'bh-name';
        name.textContent = bh.name;

        const meta = document.createElement('div');
        meta.className = 'bh-meta';
        const mass = document.createElement('span');
        mass.className = 'bh-mass';
        mass.textContent = `${formatMass(bh.mass)} M☉`;
        const type = document.createElement('span');
        type.className = 'bh-type';
        type.textContent = bh.type;
        meta.append(mass, type);

        const source = document.createElement('div');
        source.className = 'bh-source';
        source.textContent = `Source: ${bh.source}`;
        if (bh.ref) source.title = bh.ref;

        entry.append(name, meta, source);
        return entry;
    }));
}

function selectBlackHole(bh) {
    AppState.selectedBlackHole = bh;
    AppState.currentPhysics = new BlackHolePhysics(bh.mass, bh.spin);

    // Update UI
    document.querySelectorAll('.bh-entry').forEach(el => el.classList.remove('active'));
    const activeEntry = document.querySelector(`.bh-entry[data-index="${KNOWN_BLACK_HOLES.indexOf(bh)}"]`);
    if (activeEntry) activeEntry.classList.add('active');

    // Update physics display
    updatePhysicsDisplay();

    // Dispatch event for simulations page
    window.dispatchEvent(new CustomEvent('blackHoleSelected', { detail: bh }));
}

function updatePhysicsDisplay() {
    const physics = AppState.currentPhysics;
    if (!physics) return;

    const display = physics.toDisplayObject();
    const bh = AppState.selectedBlackHole;

    // Update various display elements if they exist
    const elements = {
        'displayMass': `${bh.name}: ${formatMass(bh.mass)} M☉`,
        'displayRadius': formatDistance(display.radii.schwarzschild),
        'displayTemperature': formatTemperature(display.thermodynamics.temperature),
        'displayEntropy': formatScientific(display.thermodynamics.entropy) + ' k_B',
        'displayEvapTime': formatTime(display.thermodynamics.evaporationTime),
        'displayLyapunov': formatScientific(display.quantumChaos.lyapunovBound) + ' s⁻¹',
        'displayScrambling': formatTime(display.quantumChaos.scramblingTime),
        'displayISCO': formatDistance(display.radii.isco),
        'displaySpin': `a* = ${display.spin.toFixed(3)}` + (bh.spinAssumed ? ' (assumed)' : '')
    };

    for (const [id, value] of Object.entries(elements)) {
        const el = document.getElementById(id);
        if (el) el.textContent = value;
    }
}

// === Navigation ===
function initNavigation() {
    const navbar = document.getElementById('navbar');
    const navToggle = document.getElementById('navToggle');
    const navMenu = document.getElementById('navMenu');

    if (navbar) {
        // Throttle scroll work to one update per animation frame
        let scrollTicking = false;
        window.addEventListener('scroll', () => {
            if (scrollTicking) return;
            scrollTicking = true;
            requestAnimationFrame(() => {
                navbar.classList.toggle('scrolled', window.scrollY > 50);
                updateActiveNav();
                scrollTicking = false;
            });
        }, { passive: true });
    }

    if (navToggle && navMenu) {
        navToggle.addEventListener('click', () => {
            navMenu.classList.toggle('active');
        });
    }

    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        if (anchor.id === 'modalDownload') return;

        anchor.addEventListener('click', (e) => {
            // Read href at click time; only in-page fragments are handled here
            const href = anchor.getAttribute('href') || '';
            if (!href.startsWith('#')) return;

            e.preventDefault();
            if (navMenu) navMenu.classList.remove('active');

            if (href.length === 1) {
                window.scrollTo({ top: 0, behavior: 'smooth' });
                return;
            }

            let target = null;
            try {
                target = document.querySelector(href);
            } catch {
                return; // Not a valid selector
            }
            if (target) {
                target.scrollIntoView({ behavior: 'smooth', block: 'start' });
            }
        });
    });
}

function updateActiveNav() {
    const sections = document.querySelectorAll('section[id]');
    const navLinks = document.querySelectorAll('.nav-link');

    let current = '';
    sections.forEach(section => {
        if (window.scrollY >= section.offsetTop - 200) {
            current = section.getAttribute('id');
        }
    });

    navLinks.forEach(link => {
        link.classList.toggle('active', link.getAttribute('href') === `#${current}`);
    });
}

// === Event Listeners ===
function initEventListeners() {
    // Gallery search
    const searchBtn = document.getElementById('searchBtn');
    const searchInput = document.getElementById('gallerySearch');

    if (searchBtn) {
        searchBtn.addEventListener('click', () => {
            const query = searchInput?.value.trim();
            if (query) {
                document.querySelectorAll('.tag').forEach(t => t.classList.remove('active'));
                searchGallery(query);
            }
        });
    }

    if (searchInput) {
        searchInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') searchBtn?.click();
        });
    }

    // Gallery tags
    document.querySelectorAll('.tag').forEach(tag => {
        tag.addEventListener('click', () => {
            document.querySelectorAll('.tag').forEach(t => t.classList.remove('active'));
            tag.classList.add('active');
            if (searchInput) searchInput.value = '';
            searchGallery(tag.dataset.query);
        });
    });

    // Modal
    const modalClose = document.getElementById('modalClose');
    const modalBackdrop = document.querySelector('.modal-backdrop');

    if (modalClose) modalClose.addEventListener('click', closeModal);
    if (modalBackdrop) modalBackdrop.addEventListener('click', closeModal);
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape') closeModal();
    });
}

// === API Status ===
// state: 'connecting' | 'connected' | 'error'
function updateAPIStatus(state) {
    const indicator = document.getElementById('apiIndicator');
    if (!indicator) return;

    const dot = indicator.querySelector('.indicator-dot');
    const text = indicator.querySelector('.indicator-text');
    const labels = {
        connecting: 'Connecting…',
        connected: 'NASA API Connected',
        error: 'API Error'
    };

    if (dot) {
        dot.classList.remove('connected', 'error');
        if (state === 'connected') dot.classList.add('connected');
        if (state === 'error') dot.classList.add('error');
    }
    if (text) {
        text.textContent = labels[state] || labels.connecting;
    }
}

// === NASA Image Gallery ===
function safeHttpsUrl(url) {
    try {
        const parsed = new URL(url);
        return parsed.protocol === 'https:' ? parsed.href : null;
    } catch {
        return null;
    }
}

function setGalleryMessage(text) {
    const message = document.getElementById('galleryMessage');
    if (!message) return;
    message.textContent = text;
    message.hidden = !text;
}

function createGalleryItem(item) {
    const info = item.data?.[0] || {};
    const imageUrl = safeHttpsUrl(item.links?.[0]?.href);
    if (!imageUrl) return null;

    const title = info.title || 'NASA image';
    const description = (info.description || '').substring(0, 500);
    const date = info.date_created?.split('T')[0] || '';

    const card = document.createElement('div');
    card.className = 'gallery-item';
    card.addEventListener('click', () => openModal(imageUrl, title, description, date));

    const img = document.createElement('img');
    img.src = imageUrl;
    img.alt = title;
    img.loading = 'lazy';

    const meta = document.createElement('div');
    meta.className = 'gallery-item-info';
    const titleEl = document.createElement('div');
    titleEl.className = 'gallery-item-title';
    titleEl.textContent = title;
    const dateEl = document.createElement('div');
    dateEl.className = 'gallery-item-date';
    dateEl.textContent = date || 'NASA';
    meta.append(titleEl, dateEl);

    card.append(img, meta);
    return card;
}

async function searchGallery(query) {
    const loader = document.getElementById('galleryLoader');
    const grid = document.getElementById('galleryGrid');
    if (!loader || !grid) return;

    loader.style.display = 'block';
    setGalleryMessage('');
    grid.replaceChildren();

    try {
        const results = await AppState.nasaFetcher.searchBlackHoleImages(query, 24);
        const cards = results.map(createGalleryItem).filter(Boolean);

        if (cards.length === 0) {
            setGalleryMessage('No images found.');
        } else {
            grid.replaceChildren(...cards);
        }
        updateAPIStatus('connected');
    } catch (error) {
        console.error('Gallery Error:', error);
        setGalleryMessage('Failed to load gallery. Check your connection.');
        updateAPIStatus('error');
    } finally {
        loader.style.display = 'none';
    }
}

// === Modal ===
function openModal(imageUrl, title, description, date) {
    const modal = document.getElementById('imageModal');
    if (!modal) return;

    const safeUrl = safeHttpsUrl(imageUrl);
    if (!safeUrl) return;

    const img = document.getElementById('modalImage');
    const titleEl = document.getElementById('modalTitle');
    const descEl = document.getElementById('modalDesc');
    const dateEl = document.getElementById('modalDate');
    const downloadEl = document.getElementById('modalDownload');

    if (img) {
        img.src = safeUrl;
        img.alt = title;
    }
    if (titleEl) titleEl.textContent = title;
    if (descEl) descEl.textContent = description;
    if (dateEl) dateEl.textContent = date ? formatDate(date) : '';
    if (downloadEl) downloadEl.href = safeUrl;

    modal.classList.add('active');
    document.body.style.overflow = 'hidden';
}

function closeModal() {
    const modal = document.getElementById('imageModal');
    if (modal) modal.classList.remove('active');
    document.body.style.overflow = '';
}

window.openModal = openModal;

// === Formatting Utilities ===
function formatDate(dateString) {
    try {
        const options = { year: 'numeric', month: 'long', day: 'numeric' };
        return new Date(dateString).toLocaleDateString('en-US', options);
    } catch {
        return dateString;
    }
}

// === Export for other pages ===
window.AppState = AppState;
window.selectBlackHole = selectBlackHole;
