# Figure Orbital Order

Interactive Dash app to visualize real atomic orbital shapes on configurable lattice sites.

## Run locally

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python FigureOrbitalOrder.py
```

Then open http://127.0.0.1:8050.

## Deploy

This repository includes deployment files for Render:

- `Procfile`
- `wsgi.py`
- `render.yaml`

On Render, create a new Web Service from this GitHub repo and deploy.
