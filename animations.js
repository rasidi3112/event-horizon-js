// === Scroll reveal animations (progressive enhancement) ===
// Content stays fully visible without JavaScript or with prefers-reduced-motion.
(function () {
    'use strict';

    const SELECTOR = [
        '.hero-content > *',
        '.section-header',
        '.bh-database-panel', '.bh-properties-panel',
        '.anatomy-card', '.type-card',
        '.concept-card', '.holo-card', '.app-item',
        '.law-card', '.formula-item', '.solution-card',
        '.obs-card', '.ref-card',
        '.sim-card', '.math-block'
    ].join(',');

    const reduceMotion = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    if (reduceMotion || !('IntersectionObserver' in window)) return;

    function init() {
        const observer = new IntersectionObserver(entries => {
            for (const entry of entries) {
                if (!entry.isIntersecting) continue;
                entry.target.classList.add('is-visible');
                observer.unobserve(entry.target);
            }
        }, { rootMargin: '0px 0px -8% 0px', threshold: 0.08 });

        document.querySelectorAll(SELECTOR).forEach(el => {
            // Stagger siblings of the same kind, capped so long grids do not lag
            const siblings = el.parentElement
                ? Array.from(el.parentElement.children).filter(child => child.matches(SELECTOR))
                : [];
            const index = Math.max(0, siblings.indexOf(el));
            el.style.setProperty('--reveal-delay', `${Math.min(index, 6) * 70}ms`);
            el.classList.add('reveal');
            observer.observe(el);
        });
    }

    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
    else init();
})();
