# ion-track-mc — web page (static HTML/JavaScript)

Showcase page for **ion-track-mc**, an educational event-by-event Monte Carlo code
for proton tracks in liquid water (https://github.com/juanmontiz/ion-track-mc).
It shows the capabilities of the code (animated proton tracks, cross sections,
stopping power vs NIST PSTAR and measurements, Bragg curve, radial dose, LET) and
points to the GitHub repository and the Colab notebooks.

The page **does not run simulations**: every number and track was computed beforehand
with the real code (`dev_tools/export_web_data.py`) and is stored in `data.js`.
No server-side software, no Python, no database, no external CDN.

## Deployment (for the webmaster)

Copy this folder to the web server, e.g. as `https://www.ifir-conicet.gov.ar/ion-track-mc/`,
**except the `dev_tools/` subdirectory — do NOT upload `dev_tools/`** (developer scripts,
not needed by the website). That's all. Any static web server works.

To test locally: double-click `index.html`.

### Language
Bilingual (Spanish / English) with an ES/EN switch (with flags) in the top bar. The default
language is Spanish; the visitor's last choice is remembered (localStorage); links can force it with `?lang=en` or `?lang=es`
(e.g. `https://www.ifir-conicet.gov.ar/ion-track-mc/?lang=es`).

## Files

| File | Purpose |
|------|---------|
| `index.html` | Page layout, styles and texts (both languages, `lang="en"` / `lang="es"` spans) |
| `app.js` | Language switch, track animation (canvas), interactive plots |
| `data.js` | Precomputed tracks and curves (generated, do not edit) |
| `plotly.min.js` | Plotly.js 2.35.2 plotting library (MIT license) |
| `docs/tutorial*.pdf`, `docs/tutorial*_cover.png` | Tutorial slides (copied from the repository) and their cover thumbnails |
| `dev_tools/` | **Developer tools — do NOT upload** |

## Updating the data (developers)

```bash
cd dev_tools
../../ion-track-mc/.venv/bin/python export_web_data.py               # everything, ~6 min (50-proton beam)
../../ion-track-mc/.venv/bin/python export_web_data.py --tracks-only # keep the beam already in data.js
```

The script expects the repository at `../ion-track-mc` (sibling folder) and also
copies the tutorial PDFs into `docs/`.

## How to cite

> Monti, J. M., & Quinto, M. A. (2026). *ion-track-mc: an educational event-by-event Monte Carlo code for proton tracks in liquid water* (Version 0.1.0) [Computer software]. https://github.com/juanmontiz/ion-track-mc

together with M. A. Quinto *et al.*, *Eur. Phys. J. D* **71**, 130 (2017), doi:10.1140/epjd/e2017-70709-6.

Contact: monti@ifir-conicet.gov.ar
