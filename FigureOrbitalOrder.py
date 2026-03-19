import numpy as np
from itertools import product
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, dcc, html, no_update

ORBITAL_ALIASES = {
    # Shell shortcuts.
    "s": "s",
    "p": "pz",
    "d": "dz2",
    "f": "fz3",
    "g": "gz4",
    # p orbitals.
    "pz": "pz",
    "px": "px",
    "py": "py",
    # d orbitals.
    "dz2": "dz2",
    "dxz": "dxz",
    "dyz": "dyz",
    "dx2-y2": "dx2-y2",
    "dxy": "dxy",
    # f orbitals (real tesseral forms).
    "fz3": "fz3",
    "fxz2": "fxz2",
    "fyz2": "fyz2",
    "fzx2-zy2": "fzx2-zy2",
    "fxyz": "fxyz",
    "fx3-3xy2": "fx3-3xy2",
    "f3x2y-y3": "f3x2y-y3",
    # g orbitals (real tesseral forms).
    "gz4": "gz4",
    "gxz3": "gxz3",
    "gyz3": "gyz3",
    "gz2x2-y2": "gz2x2-y2",
    "gxyz2": "gxyz2",
    "gx4-6x2y2+y4": "gx4-6x2y2+y4",
    "gy4-6x2y2+x4": "gy4-6x2y2+x4",
    "gzx3-zy3": "gzx3-zy3",
    "gxyx2-y2": "gxyx2-y2",
}

# Reflection parity exponents (mod 2) for each canonical real orbital.
# Reflection factor is sx^ex * sy^ey * sz^ez where each s in {-1, +1}.
ORBITAL_REFLECTION_EXPONENTS = {
    "s": (0, 0, 0),
    "px": (1, 0, 0),
    "py": (0, 1, 0),
    "pz": (0, 0, 1),
    "dz2": (0, 0, 0),
    "dxz": (1, 0, 1),
    "dyz": (0, 1, 1),
    "dx2-y2": (0, 0, 0),
    "dxy": (1, 1, 0),
    "fz3": (0, 0, 1),
    "fxz2": (1, 0, 0),
    "fyz2": (0, 1, 0),
    "fzx2-zy2": (0, 0, 1),
    "fxyz": (1, 1, 1),
    "fx3-3xy2": (1, 0, 0),
    "f3x2y-y3": (0, 1, 0),
    "gz4": (0, 0, 0),
    "gxz3": (1, 0, 1),
    "gyz3": (0, 1, 1),
    "gz2x2-y2": (0, 0, 0),
    "gxyz2": (1, 1, 0),
    "gx4-6x2y2+y4": (0, 0, 0),
    "gy4-6x2y2+x4": (1, 1, 0),
    "gzx3-zy3": (1, 0, 1),
    "gxyx2-y2": (0, 1, 1),
}

def orbital_angular(theta, phi, orbital_type="dz2"):
    """Return angular part from explicit real-orbital polynomial formulas."""
    key = orbital_type.lower()
    canonical = ORBITAL_ALIASES.get(key)
    if canonical is None:
        supported = ", ".join(sorted(ORBITAL_ALIASES.keys()))
        raise ValueError(f"Unsupported orbital '{orbital_type}'. Choose one of: {supported}")

    # Unit-sphere Cartesian coordinates from angular grid.
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    # Real (tesseral) orbital angular forms up to g shell.
    formulas = {
        # s
        "s": np.ones_like(theta),
        # p
        "px": x,
        "py": y,
        "pz": z,
        # d
        "dz2": 0.5 * (3.0 * z**2 - 1.0),
        "dxz": x * z,
        "dyz": y * z,
        "dx2-y2": x**2 - y**2,
        "dxy": x * y,
        # f
        "fz3": 0.5 * (5.0 * z**3 - 3.0 * z),
        "fxz2": x * (5.0 * z**2 - 1.0),
        "fyz2": y * (5.0 * z**2 - 1.0),
        "fzx2-zy2": z * (x**2 - y**2),
        "fxyz": x * y * z,
        "fx3-3xy2": x * (x**2 - 3.0 * y**2),
        "f3x2y-y3": y * (3.0 * x**2 - y**2),
        # g
        "gz4": 35.0 * z**4 - 30.0 * z**2 + 3.0,
        "gxz3": x * (7.0 * z**3 - 3.0 * z),
        "gyz3": y * (7.0 * z**3 - 3.0 * z),
        "gz2x2-y2": (x**2 - y**2) * (7.0 * z**2 - 1.0),
        "gxyz2": x * y * (7.0 * z**2 - 1.0),
        "gx4-6x2y2+y4": x**4 - 6.0 * x**2 * y**2 + y**4,
        "gy4-6x2y2+x4": 4.0 * x * y * (x**2 - y**2),
        "gzx3-zy3": z * (x**3 - 3.0 * x * y**2),
        "gxyx2-y2": z * (3.0 * x**2 * y - y**3),
    }

    return formulas[canonical]


def orbital_reflection_parity(orbital_type, sx, sy, sz):
    """Return +/-1 for orbital sign change under coordinate reflection."""
    key = orbital_type.lower()
    canonical = ORBITAL_ALIASES.get(key)
    if canonical is None:
        return 1.0

    exponents = ORBITAL_REFLECTION_EXPONENTS.get(canonical, (0, 0, 0))
    ex, ey, ez = exponents

    parity = 1.0
    if ex % 2 == 1:
        parity *= float(sx)
    if ey % 2 == 1:
        parity *= float(sy)
    if ez % 2 == 1:
        parity *= float(sz)
    return parity


def compute_scene_ranges(sites, orbital_extent=1.0, padding=0.4, min_half_span=2.0):
    """Compute symmetric 3D axis ranges that include all orbitals and preview."""
    if not sites:
        return {
            "x": [-min_half_span, min_half_span],
            "y": [-min_half_span, min_half_span],
            "z": [-min_half_span, min_half_span],
        }

    centers = np.array([site["pos"] for site in sites], dtype=float)
    mins = centers.min(axis=0) - orbital_extent - padding
    maxs = centers.max(axis=0) + orbital_extent + padding

    center_point = 0.5 * (mins + maxs)
    half_span = max(0.5 * np.max(maxs - mins), min_half_span)

    return {
        "x": [center_point[0] - half_span, center_point[0] + half_span],
        "y": [center_point[1] - half_span, center_point[1] + half_span],
        "z": [center_point[2] - half_span, center_point[2] + half_span],
    }


def rotate_coordinates(x, y, z, rx_deg=0.0, ry_deg=0.0, rz_deg=0.0):
    """Rotate coordinate grids around x, y, z axes by degrees."""
    rx = np.deg2rad(float(rx_deg or 0.0))
    ry = np.deg2rad(float(ry_deg or 0.0))
    rz = np.deg2rad(float(rz_deg or 0.0))

    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)

    # Rotation order: X, then Y, then Z.
    y1 = y * cx - z * sx
    z1 = y * sx + z * cx
    x1 = x

    x2 = x1 * cy + z1 * sy
    z2 = -x1 * sy + z1 * cy
    y2 = y1

    x3 = x2 * cz - y2 * sz
    y3 = x2 * sz + y2 * cz
    z3 = z2

    return x3, y3, z3


def _reflect_coordinate(value, sign, plane_offset):
    """Reflect scalar across an offset coordinate plane when sign is -1."""
    if sign >= 0:
        return value
    return 2.0 * plane_offset - value


def expand_sites_with_mirror_planes(sites, mirror_planes, mirror_offsets=None):
    """Duplicate sites by reflecting them across selected coordinate planes."""
    sites = sites or []
    mirror_planes = set(mirror_planes or [])
    mirror_offsets = mirror_offsets or {}
    x_plane = float(mirror_offsets.get("yz", 0.0) or 0.0)
    y_plane = float(mirror_offsets.get("xz", 0.0) or 0.0)
    z_plane = float(mirror_offsets.get("xy", 0.0) or 0.0)

    x_signs = [-1.0, 1.0] if "yz" in mirror_planes else [1.0]
    y_signs = [-1.0, 1.0] if "xz" in mirror_planes else [1.0]
    z_signs = [-1.0, 1.0] if "xy" in mirror_planes else [1.0]

    sign_combos = list(product(x_signs, y_signs, z_signs))
    expanded_sites = []

    for site in sites:
        pos = site["pos"]
        for sx, sy, sz in sign_combos:
            det_sign = sx * sy * sz
            new_site = dict(site)
            new_site["pos"] = [
                _reflect_coordinate(pos[0], sx, x_plane),
                _reflect_coordinate(pos[1], sy, y_plane),
                _reflect_coordinate(pos[2], sz, z_plane),
            ]
            new_site["mirror_sign"] = [sx, sy, sz]
            new_site["mirror_det"] = det_sign
            expanded_sites.append(new_site)

    return expanded_sites


def apply_superposed_opposite_opacity(sites, superposed_opacity=0.05, tol=1e-9):
    """Dim exact overlapping orbitals when opposite effective phases coexist."""
    sites = sites or []
    grouped_indices = {}

    def _round_tuple(values, ndigits=9):
        return tuple(round(float(v), ndigits) for v in values)

    for idx, site in enumerate(sites):
        canonical_type = ORBITAL_ALIASES.get(str(site.get("type", "")).lower(), str(site.get("type", "")).lower())
        key = (
            _round_tuple(site.get("pos", [0.0, 0.0, 0.0])),
            canonical_type,
            _round_tuple([
                site.get("rot_x", 0.0),
                site.get("rot_y", 0.0),
                site.get("rot_z", 0.0),
            ]),
            str(site.get("multipole_kind", "electric")).lower(),
        )
        grouped_indices.setdefault(key, []).append(idx)

    updated_sites = [dict(site) for site in sites]

    for indices in grouped_indices.values():
        if len(indices) < 2:
            continue

        has_positive = False
        has_negative = False
        for i in indices:
            site = updated_sites[i]
            sx, sy, sz = site.get("mirror_sign", [1.0, 1.0, 1.0])
            orbital_parity = orbital_reflection_parity(site.get("type", ""), sx, sy, sz)
            base_phase = float(site.get("phase", 1.0) or 1.0)
            multipole_kind = str(site.get("multipole_kind", "electric")).lower()
            mirror_det = float(site.get("mirror_det", 1.0) or 1.0)
            magnetic_axial = mirror_det if multipole_kind == "magnetic" else 1.0
            effective_phase = base_phase * orbital_parity * magnetic_axial

            if effective_phase > tol:
                has_positive = True
            elif effective_phase < -tol:
                has_negative = True

            if has_positive and has_negative:
                break

        if has_positive and has_negative:
            for i in indices:
                updated_sites[i]["opacity"] = float(superposed_opacity)

    return updated_sites


def add_mirror_plane_surfaces(fig, scene_ranges, mirror_planes, mirror_offsets=None, opacity=0.18):
    """Draw translucent mirror planes in the scene."""
    mirror_planes = set(mirror_planes or [])
    mirror_offsets = mirror_offsets or {}

    x_min, x_max = scene_ranges["x"]
    y_min, y_max = scene_ranges["y"]
    z_min, z_max = scene_ranges["z"]

    if "xy" in mirror_planes:
        z0 = float(mirror_offsets.get("xy", 0.0) or 0.0)
        fig.add_trace(
            go.Surface(
                x=np.array([[x_min, x_max], [x_min, x_max]]),
                y=np.array([[y_min, y_min], [y_max, y_max]]),
                z=np.array([[z0, z0], [z0, z0]]),
                surfacecolor=np.zeros((2, 2)),
                colorscale=[[0.0, "#2ca02c"], [1.0, "#2ca02c"]],
                showscale=False,
                opacity=opacity,
                hoverinfo="skip",
            )
        )

    if "xz" in mirror_planes:
        y0 = float(mirror_offsets.get("xz", 0.0) or 0.0)
        fig.add_trace(
            go.Surface(
                x=np.array([[x_min, x_max], [x_min, x_max]]),
                y=np.array([[y0, y0], [y0, y0]]),
                z=np.array([[z_min, z_min], [z_max, z_max]]),
                surfacecolor=np.zeros((2, 2)),
                colorscale=[[0.0, "#ff7f0e"], [1.0, "#ff7f0e"]],
                showscale=False,
                opacity=opacity,
                hoverinfo="skip",
            )
        )

    if "yz" in mirror_planes:
        x0 = float(mirror_offsets.get("yz", 0.0) or 0.0)
        fig.add_trace(
            go.Surface(
                x=np.array([[x0, x0], [x0, x0]]),
                y=np.array([[y_min, y_max], [y_min, y_max]]),
                z=np.array([[z_min, z_min], [z_max, z_max]]),
                surfacecolor=np.zeros((2, 2)),
                colorscale=[[0.0, "#1f77b4"], [1.0, "#1f77b4"]],
                showscale=False,
                opacity=opacity,
                hoverinfo="skip",
            )
        )


def plot_lattice_orbitals(sites, resolution=80, mirror_planes=None, mirror_offsets=None):
    fig = go.Figure()
    scene_ranges = compute_scene_ranges(sites)
    
    # Grid for the orbital surface
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    U, V = np.meshgrid(u, v)

    for site in sites:
        pos = site['pos']
        orb_type = site['type']
        phase = site.get('phase', 1.0)
        multipole_kind = site.get('multipole_kind', 'electric')
        rot_x = site.get('rot_x', 0.0)
        rot_y = site.get('rot_y', 0.0)
        rot_z = site.get('rot_z', 0.0)
        mirror_sign = site.get('mirror_sign', [1.0, 1.0, 1.0])
        mirror_det = float(site.get('mirror_det', 1.0))
        color_pos = site.get('color_pos', 'red')
        color_neg = site.get('color_neg', 'blue')
        opacity = site.get('opacity', 0.85)

        sx, sy, sz = mirror_sign

        # Geometry is already mirrored by applying mirror_sign to coordinates
        # below. Applying orbital reflection parity again here would double-count
        # the mirror effect and cancel expected color/phase inversion.
        # Magnetic multipoles still require the axial det factor.
        magnetic_axial = mirror_det if multipole_kind == 'magnetic' else 1.0
        mirror_phase = magnetic_axial
        
        # Calculate angular shape for the selected orbital.
        R_raw = phase * mirror_phase * orbital_angular(V, U, orb_type)
        R = R_raw / np.max(np.abs(R_raw))
        
        # Convert to Cartesian for plotting, centered at 'pos'
        X = np.abs(R) * np.sin(V) * np.cos(U)
        Y = np.abs(R) * np.sin(V) * np.sin(U)
        Z = np.abs(R) * np.cos(V)

        X, Y, Z = rotate_coordinates(X, Y, Z, rot_x, rot_y, rot_z)
        X = mirror_sign[0] * X
        Y = mirror_sign[1] * Y
        Z = mirror_sign[2] * Z

        # Split the surface by phase: positive and negative lobes.
        pos_mask = R >= 0
        neg_mask = R < 0

        X_pos = np.where(pos_mask, X, np.nan) + pos[0]
        Y_pos = np.where(pos_mask, Y, np.nan) + pos[1]
        Z_pos = np.where(pos_mask, Z, np.nan) + pos[2]

        X_neg = np.where(neg_mask, X, np.nan) + pos[0]
        Y_neg = np.where(neg_mask, Y, np.nan) + pos[1]
        Z_neg = np.where(neg_mask, Z, np.nan) + pos[2]

        fig.add_trace(
            go.Surface(
                x=X_pos,
                y=Y_pos,
                z=Z_pos,
                surfacecolor=np.zeros_like(X_pos),
                colorscale=[[0.0, color_pos], [1.0, color_pos]],
                showscale=False,
                opacity=opacity,
                hoverinfo='skip',
            )
        )

        fig.add_trace(
            go.Surface(
                x=X_neg,
                y=Y_neg,
                z=Z_neg,
                surfacecolor=np.zeros_like(X_neg),
                colorscale=[[0.0, color_neg], [1.0, color_neg]],
                showscale=False,
                opacity=opacity,
                hoverinfo='skip',
            )
        )

    fig.update_layout(
        scene=dict(
            xaxis=dict(title='X (A)', range=scene_ranges["x"]),
            yaxis=dict(title='Y (A)', range=scene_ranges["y"]),
            zaxis=dict(title='Z (A)', range=scene_ranges["z"]),
            aspectmode='cube',
        ),
        margin=dict(l=0, r=0, b=0, t=30),
    )
    add_mirror_plane_surfaces(fig, scene_ranges, mirror_planes, mirror_offsets)
    return fig

def build_orbital_app(resolution=80):
    app = Dash(__name__)

    orbital_options = sorted(set(ORBITAL_ALIASES.keys()))

    compact_control_style = {"fontSize": "13px", "marginBottom": "2px", "display": "block"}
    compact_input_style = {"width": "100%", "height": "32px", "fontSize": "13px", "marginBottom": "8px"}
    compact_button_style = {"fontSize": "12px", "padding": "6px 8px", "height": "32px"}

    app.layout = html.Div(
        style={"fontFamily": "Arial, sans-serif", "padding": "10px"},
        children=[
            dcc.Store(id="sites-store", data=[]),
            html.Div(
                style={
                    "display": "flex",
                    "gap": "10px",
                    "alignItems": "stretch",
                    "flexWrap": "wrap",
                },
                children=[
                    html.Div(
                        style={
                            "flex": "0 0 280px",
                            "maxWidth": "320px",
                            "minWidth": "240px",
                            "border": "1px solid #d5d5d5",
                            "borderRadius": "8px",
                            "padding": "10px",
                            "background": "#fafafa",
                        },
                        children=[
                            html.H4("Orbital Controls", style={"margin": "0 0 8px 0", "fontSize": "16px"}),
                            html.Label("Orbital", style=compact_control_style),
                            dcc.Dropdown(
                                id="orbital-type",
                                options=[{"label": o, "value": o} for o in orbital_options],
                                value="dz2",
                                clearable=False,
                                style={"marginBottom": "8px", "fontSize": "13px"},
                            ),
                            html.Label("X", style=compact_control_style),
                            dcc.Input(id="x-pos", type="number", value=0.0, step=0.1, style=compact_input_style),
                            html.Label("Y", style=compact_control_style),
                            dcc.Input(id="y-pos", type="number", value=0.0, step=0.1, style=compact_input_style),
                            html.Label("Z", style=compact_control_style),
                            dcc.Input(id="z-pos", type="number", value=0.0, step=0.1, style=compact_input_style),
                            html.Label("Rotate X (deg)", style=compact_control_style),
                            dcc.Input(id="rot-x", type="number", value=0.0, step=1.0, style=compact_input_style),
                            html.Label("Rotate Y (deg)", style=compact_control_style),
                            dcc.Input(id="rot-y", type="number", value=0.0, step=1.0, style=compact_input_style),
                            html.Label("Rotate Z (deg)", style=compact_control_style),
                            dcc.Input(id="rot-z", type="number", value=0.0, step=1.0, style=compact_input_style),
                            html.Label("+ Color", style=compact_control_style),
                            dcc.Input(id="color-pos", type="text", value="red", style=compact_input_style),
                            html.Label("- Color", style=compact_control_style),
                            dcc.Input(id="color-neg", type="text", value="blue", style=compact_input_style),
                            html.Label("Phase", style=compact_control_style),
                            dcc.Dropdown(
                                id="phase",
                                options=[{"label": "+1", "value": 1.0}, {"label": "-1", "value": -1.0}],
                                value=1.0,
                                clearable=False,
                                style={"marginBottom": "8px", "fontSize": "13px"},
                            ),
                            html.Label("Multipole", style=compact_control_style),
                            dcc.Dropdown(
                                id="multipole-kind",
                                options=[
                                    {"label": "Electric", "value": "electric"},
                                    {"label": "Magnetic", "value": "magnetic"},
                                ],
                                value="electric",
                                clearable=False,
                                style={"marginBottom": "8px", "fontSize": "13px"},
                            ),
                            html.Div(
                                style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "6px", "marginBottom": "8px"},
                                children=[
                                    html.Button("Add", id="add-btn", n_clicks=0, style=compact_button_style),
                                    html.Button("Swap", id="swap-colors-btn", n_clicks=0, style=compact_button_style),
                                    html.Button("Undo", id="undo-btn", n_clicks=0, style=compact_button_style),
                                    html.Button("Trash", id="trash-btn", n_clicks=0, style=compact_button_style),
                                    html.Button(
                                        "Clear All",
                                        id="clear-btn",
                                        n_clicks=0,
                                        style={**compact_button_style, "gridColumn": "1 / span 2"},
                                    ),
                                ],
                            ),
                            dcc.Checklist(
                                id="preview-toggle",
                                options=[{"label": "Preview Orbital", "value": "show"}],
                                value=["show"],
                                style={"fontSize": "13px", "marginBottom": "8px"},
                            ),
                            dcc.Checklist(
                                id="mirror-planes",
                                options=[
                                    {"label": "Mirror XY (z -> -z)", "value": "xy"},
                                    {"label": "Mirror XZ (y -> -y)", "value": "xz"},
                                    {"label": "Mirror YZ (x -> -x)", "value": "yz"},
                                ],
                                value=[],
                                style={"fontSize": "13px", "marginBottom": "8px"},
                            ),
                            html.Label("XY plane z", style=compact_control_style),
                            dcc.Input(id="mirror-xy-z", type="number", value=0.0, step=0.1, style=compact_input_style),
                            html.Label("XZ plane y", style=compact_control_style),
                            dcc.Input(id="mirror-xz-y", type="number", value=0.0, step=0.1, style=compact_input_style),
                            html.Label("YZ plane x", style=compact_control_style),
                            dcc.Input(id="mirror-yz-x", type="number", value=0.0, step=0.1, style=compact_input_style),
                            dcc.Dropdown(
                                id="site-select",
                                options=[],
                                placeholder="Select orbital to trash",
                                style={"marginBottom": "8px", "fontSize": "13px"},
                            ),
                            html.Div(id="status", style={"fontSize": "12px", "lineHeight": "1.35"}),
                        ],
                    ),
                    html.Div(
                        style={"flex": "1 1 640px", "minWidth": "320px"},
                        children=[
                            dcc.Graph(id="orbital-graph", style={"height": "92vh"}),
                        ],
                    ),
                ],
            ),
        ],
    )

    @app.callback(
        Output("sites-store", "data"),
        Output("status", "children"),
        Input("add-btn", "n_clicks"),
        Input("swap-colors-btn", "n_clicks"),
        Input("undo-btn", "n_clicks"),
        Input("trash-btn", "n_clicks"),
        Input("clear-btn", "n_clicks"),
        State("orbital-type", "value"),
        State("x-pos", "value"),
        State("y-pos", "value"),
        State("z-pos", "value"),
        State("rot-x", "value"),
        State("rot-y", "value"),
        State("rot-z", "value"),
        State("color-pos", "value"),
        State("color-neg", "value"),
        State("phase", "value"),
        State("multipole-kind", "value"),
        State("site-select", "value"),
        State("sites-store", "data"),
        prevent_initial_call=True,
    )
    def update_sites(add_clicks, swap_clicks, undo_clicks, trash_clicks, clear_clicks, orbital_type, x_pos, y_pos, z_pos,
                     rot_x, rot_y, rot_z, color_pos, color_neg, phase, multipole_kind, site_select, sites):
        from dash import ctx
        trigger = ctx.triggered_id
        if trigger is None:
            return sites, "No action performed."

        sites = sites or []

        if trigger == "add-btn":
            if x_pos is None or y_pos is None or z_pos is None:
                return sites, "Please provide numeric X, Y, Z values."
            if orbital_type not in ORBITAL_ALIASES:
                return sites, f"Unsupported orbital: {orbital_type}"
            new_site = {
                "pos": [float(x_pos), float(y_pos), float(z_pos)],
                "type": orbital_type,
                "phase": float(phase),
                "rot_x": float(rot_x or 0.0),
                "rot_y": float(rot_y or 0.0),
                "rot_z": float(rot_z or 0.0),
                "multipole_kind": multipole_kind or "electric",
                "color_pos": color_pos or "red",
                "color_neg": color_neg or "blue",
            }
            sites = [*sites, new_site]
            return sites, (
                f"Added {orbital_type} at ({x_pos}, {y_pos}, {z_pos}) with "
                f"rotation ({rot_x or 0.0}, {rot_y or 0.0}, {rot_z or 0.0}) deg as "
                f"{multipole_kind or 'electric'}. Total: {len(sites)}"
            )

        if trigger == "swap-colors-btn":
            # Input colors are swapped by a separate callback; here we also swap
            # the currently selected placed orbital to reflect "current" object.
            if site_select is None or not sites:
                return sites, "Swapped input colors."
            idx = int(site_select)
            if idx < 0 or idx >= len(sites):
                return sites, "Swapped input colors."

            updated_sites = [*sites]
            current = dict(updated_sites[idx])
            current_pos_color = current.get("color_pos", "red")
            current_neg_color = current.get("color_neg", "blue")
            current["color_pos"] = current_neg_color
            current["color_neg"] = current_pos_color
            updated_sites[idx] = current
            return updated_sites, f"Swapped colors for selected orbital {idx}."

        if trigger == "undo-btn":
            if not sites:
                return sites, "Nothing to undo."
            sites = sites[:-1]
            return sites, f"Removed last orbital. Total: {len(sites)}"

        if trigger == "trash-btn":
            if not sites:
                return sites, "Nothing to delete."
            if site_select is None:
                return sites, "Select an orbital to trash first."
            idx = int(site_select)
            if idx < 0 or idx >= len(sites):
                return sites, "Selected orbital no longer exists."
            removed = sites[idx]
            sites = [site for i, site in enumerate(sites) if i != idx]
            return sites, (
                f"Trashed orbital {idx}: {removed['type']} at ({removed['pos'][0]}, "
                f"{removed['pos'][1]}, {removed['pos'][2]}). Total: {len(sites)}"
            )

        if trigger == "clear-btn":
            return [], "Cleared all orbitals."

        return sites, "No action performed."

    @app.callback(
        Output("color-pos", "value"),
        Output("color-neg", "value"),
        Input("swap-colors-btn", "n_clicks"),
        State("color-pos", "value"),
        State("color-neg", "value"),
        prevent_initial_call=True,
    )
    def swap_colors(_, color_pos, color_neg):
        return (color_neg or "blue"), (color_pos or "red")

    @app.callback(
        Output("site-select", "options"),
        Output("site-select", "value"),
        Input("sites-store", "data"),
        State("site-select", "value"),
    )
    def refresh_site_selector(sites, selected_value):
        sites = sites or []
        options = [
            {
                "label": (
                    f"{i}: {site.get('multipole_kind', 'electric')[0].upper()}-{site['type']} @ "
                    f"({site['pos'][0]:.2f}, {site['pos'][1]:.2f}, {site['pos'][2]:.2f}) "
                    f"r=({site.get('rot_x', 0.0):.0f}, {site.get('rot_y', 0.0):.0f}, {site.get('rot_z', 0.0):.0f})"
                ),
                "value": i,
            }
            for i, site in enumerate(sites)
        ]

        if not options:
            return [], None

        if selected_value is None or int(selected_value) >= len(options):
            return options, len(options) - 1

        return options, selected_value

    @app.callback(
        Output("orbital-graph", "figure"),
        Input("sites-store", "data"),
        Input("preview-toggle", "value"),
        Input("mirror-planes", "value"),
        Input("mirror-xy-z", "value"),
        Input("mirror-xz-y", "value"),
        Input("mirror-yz-x", "value"),
        Input("orbital-type", "value"),
        Input("x-pos", "value"),
        Input("y-pos", "value"),
        Input("z-pos", "value"),
        Input("rot-x", "value"),
        Input("rot-y", "value"),
        Input("rot-z", "value"),
        Input("color-pos", "value"),
        Input("color-neg", "value"),
        Input("phase", "value"),
        Input("multipole-kind", "value"),
    )
    def update_graph(sites, preview_toggle, mirror_planes, mirror_xy_z, mirror_xz_y, mirror_yz_x, orbital_type, x_pos, y_pos, z_pos, rot_x, rot_y, rot_z, color_pos, color_neg, phase, multipole_kind):
        sites = sites or []
        show_preview = preview_toggle is not None and "show" in preview_toggle
        mirror_offsets = {
            "xy": float(mirror_xy_z or 0.0),
            "xz": float(mirror_xz_y or 0.0),
            "yz": float(mirror_yz_x or 0.0),
        }

        if show_preview and orbital_type in ORBITAL_ALIASES and None not in (x_pos, y_pos, z_pos):
            preview_site = {
                "pos": [float(x_pos), float(y_pos), float(z_pos)],
                "type": orbital_type,
                "phase": float(phase),
                "rot_x": float(rot_x or 0.0),
                "rot_y": float(rot_y or 0.0),
                "rot_z": float(rot_z or 0.0),
                "multipole_kind": multipole_kind or "electric",
                "color_pos": color_pos or "red",
                "color_neg": color_neg or "blue",
                "opacity": 0.45,
            }
            sites_to_draw = [*sites, preview_site]
        else:
            sites_to_draw = sites

        sites_to_draw = expand_sites_with_mirror_planes(sites_to_draw, mirror_planes, mirror_offsets)
        sites_to_draw = apply_superposed_opposite_opacity(sites_to_draw, superposed_opacity=0.1)
        scene_ranges = compute_scene_ranges(sites_to_draw)

        if not sites_to_draw:
            fig = go.Figure()
            fig.update_layout(
                scene=dict(
                    xaxis=dict(title="X (A)", range=scene_ranges["x"]),
                    yaxis=dict(title="Y (A)", range=scene_ranges["y"]),
                    zaxis=dict(title="Z (A)", range=scene_ranges["z"]),
                    aspectmode="cube",
                ),
                margin=dict(l=0, r=0, b=0, t=30),
                annotations=[
                    dict(
                        text="Add an orbital from the controls above.",
                        xref="paper",
                        yref="paper",
                        x=0.01,
                        y=0.99,
                        showarrow=False,
                    )
                ],
            )
            add_mirror_plane_surfaces(fig, scene_ranges, mirror_planes, mirror_offsets)
            return fig
        return plot_lattice_orbitals(
            sites_to_draw,
            resolution=resolution,
            mirror_planes=mirror_planes,
            mirror_offsets=mirror_offsets,
        )

    return app


if __name__ == "__main__":
    app = build_orbital_app(resolution=100)
    app.run(debug=False)