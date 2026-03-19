const ORBITAL_ALIASES = {
  s: "s",
  p: "pz",
  d: "dz2",
  f: "fz3",
  g: "gz4",
  pz: "pz",
  px: "px",
  py: "py",
  dz2: "dz2",
  dxz: "dxz",
  dyz: "dyz",
  "dx2-y2": "dx2-y2",
  dxy: "dxy",
  fz3: "fz3",
  fxz2: "fxz2",
  fyz2: "fyz2",
  "fzx2-zy2": "fzx2-zy2",
  fxyz: "fxyz",
  "fx3-3xy2": "fx3-3xy2",
  "f3x2y-y3": "f3x2y-y3",
  gz4: "gz4",
  gxz3: "gxz3",
  gyz3: "gyz3",
  "gz2x2-y2": "gz2x2-y2",
  gxyz2: "gxyz2",
  "gx4-6x2y2+y4": "gx4-6x2y2+y4",
  "gy4-6x2y2+x4": "gy4-6x2y2+x4",
  "gzx3-zy3": "gzx3-zy3",
  "gxyx2-y2": "gxyx2-y2"
};

const ORBITAL_REFLECTION_EXPONENTS = {
  s: [0, 0, 0],
  px: [1, 0, 0],
  py: [0, 1, 0],
  pz: [0, 0, 1],
  dz2: [0, 0, 0],
  dxz: [1, 0, 1],
  dyz: [0, 1, 1],
  "dx2-y2": [0, 0, 0],
  dxy: [1, 1, 0],
  fz3: [0, 0, 1],
  fxz2: [1, 0, 0],
  fyz2: [0, 1, 0],
  "fzx2-zy2": [0, 0, 1],
  fxyz: [1, 1, 1],
  "fx3-3xy2": [1, 0, 0],
  "f3x2y-y3": [0, 1, 0],
  gz4: [0, 0, 0],
  gxz3: [1, 0, 1],
  gyz3: [0, 1, 1],
  "gz2x2-y2": [0, 0, 0],
  gxyz2: [1, 1, 0],
  "gx4-6x2y2+y4": [0, 0, 0],
  "gy4-6x2y2+x4": [1, 1, 0],
  "gzx3-zy3": [1, 0, 1],
  "gxyx2-y2": [0, 1, 1]
};

const RESOLUTION = 58;
const EPSILON = 1e-9;
let sites = [];

const dom = {
  orbitalType: document.getElementById("orbitalType"),
  xPos: document.getElementById("xPos"),
  yPos: document.getElementById("yPos"),
  zPos: document.getElementById("zPos"),
  rotX: document.getElementById("rotX"),
  rotY: document.getElementById("rotY"),
  rotZ: document.getElementById("rotZ"),
  colorPos: document.getElementById("colorPos"),
  colorNeg: document.getElementById("colorNeg"),
  phase: document.getElementById("phase"),
  multipoleKind: document.getElementById("multipoleKind"),
  previewToggle: document.getElementById("previewToggle"),
  mirrorXY: document.getElementById("mirrorXY"),
  mirrorXZ: document.getElementById("mirrorXZ"),
  mirrorYZ: document.getElementById("mirrorYZ"),
  mirrorXyZ: document.getElementById("mirrorXyZ"),
  mirrorXzY: document.getElementById("mirrorXzY"),
  mirrorYzX: document.getElementById("mirrorYzX"),
  siteSelect: document.getElementById("siteSelect"),
  addBtn: document.getElementById("addBtn"),
  swapBtn: document.getElementById("swapBtn"),
  undoBtn: document.getElementById("undoBtn"),
  trashBtn: document.getElementById("trashBtn"),
  clearBtn: document.getElementById("clearBtn"),
  status: document.getElementById("status"),
  graph: document.getElementById("orbitalGraph")
};

function toNumber(value, fallback = 0) {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : fallback;
}

function canonicalOrbital(key) {
  if (!key) {
    return null;
  }
  return ORBITAL_ALIASES[String(key).toLowerCase()] || null;
}

function orbitalValue(canonical, x, y, z) {
  switch (canonical) {
    case "s": return 1;
    case "px": return x;
    case "py": return y;
    case "pz": return z;
    case "dz2": return 0.5 * (3 * z * z - 1);
    case "dxz": return x * z;
    case "dyz": return y * z;
    case "dx2-y2": return x * x - y * y;
    case "dxy": return x * y;
    case "fz3": return 0.5 * (5 * z * z * z - 3 * z);
    case "fxz2": return x * (5 * z * z - 1);
    case "fyz2": return y * (5 * z * z - 1);
    case "fzx2-zy2": return z * (x * x - y * y);
    case "fxyz": return x * y * z;
    case "fx3-3xy2": return x * (x * x - 3 * y * y);
    case "f3x2y-y3": return y * (3 * x * x - y * y);
    case "gz4": return 35 * z ** 4 - 30 * z * z + 3;
    case "gxz3": return x * (7 * z ** 3 - 3 * z);
    case "gyz3": return y * (7 * z ** 3 - 3 * z);
    case "gz2x2-y2": return (x * x - y * y) * (7 * z * z - 1);
    case "gxyz2": return x * y * (7 * z * z - 1);
    case "gx4-6x2y2+y4": return x ** 4 - 6 * x * x * y * y + y ** 4;
    case "gy4-6x2y2+x4": return 4 * x * y * (x * x - y * y);
    case "gzx3-zy3": return z * (x ** 3 - 3 * x * y * y);
    case "gxyx2-y2": return z * (3 * x * x * y - y ** 3);
    default: return 0;
  }
}

function rotatePoint(x, y, z, rxDeg, ryDeg, rzDeg) {
  const rx = (rxDeg * Math.PI) / 180;
  const ry = (ryDeg * Math.PI) / 180;
  const rz = (rzDeg * Math.PI) / 180;

  const cx = Math.cos(rx);
  const sx = Math.sin(rx);
  const cy = Math.cos(ry);
  const sy = Math.sin(ry);
  const cz = Math.cos(rz);
  const sz = Math.sin(rz);

  const y1 = y * cx - z * sx;
  const z1 = y * sx + z * cx;
  const x1 = x;

  const x2 = x1 * cy + z1 * sy;
  const z2 = -x1 * sy + z1 * cy;
  const y2 = y1;

  const x3 = x2 * cz - y2 * sz;
  const y3 = x2 * sz + y2 * cz;
  const z3 = z2;

  return [x3, y3, z3];
}

function reflectCoordinate(value, sign, planeOffset) {
  if (sign >= 0) {
    return value;
  }
  return 2 * planeOffset - value;
}

function orbitalReflectionParity(orbitalType, sx, sy, sz) {
  const canonical = canonicalOrbital(orbitalType);
  if (!canonical) {
    return 1;
  }

  const exponents = ORBITAL_REFLECTION_EXPONENTS[canonical] || [0, 0, 0];
  let parity = 1;

  if (exponents[0] % 2 === 1) {
    parity *= sx;
  }
  if (exponents[1] % 2 === 1) {
    parity *= sy;
  }
  if (exponents[2] % 2 === 1) {
    parity *= sz;
  }

  return parity;
}

function selectedMirrorPlanes() {
  const planes = [];
  if (dom.mirrorXY.checked) {
    planes.push("xy");
  }
  if (dom.mirrorXZ.checked) {
    planes.push("xz");
  }
  if (dom.mirrorYZ.checked) {
    planes.push("yz");
  }
  return planes;
}

function mirrorOffsets() {
  return {
    xy: toNumber(dom.mirrorXyZ.value, 0),
    xz: toNumber(dom.mirrorXzY.value, 0),
    yz: toNumber(dom.mirrorYzX.value, 0)
  };
}

function expandSitesWithMirrorPlanes(baseSites, mirrorPlanes, offsets) {
  const xSigns = mirrorPlanes.includes("yz") ? [-1, 1] : [1];
  const ySigns = mirrorPlanes.includes("xz") ? [-1, 1] : [1];
  const zSigns = mirrorPlanes.includes("xy") ? [-1, 1] : [1];

  const expanded = [];

  for (const site of baseSites) {
    for (const sx of xSigns) {
      for (const sy of ySigns) {
        for (const sz of zSigns) {
          const det = sx * sy * sz;
          expanded.push({
            ...site,
            pos: [
              reflectCoordinate(site.pos[0], sx, offsets.yz),
              reflectCoordinate(site.pos[1], sy, offsets.xz),
              reflectCoordinate(site.pos[2], sz, offsets.xy)
            ],
            mirrorSign: [sx, sy, sz],
            mirrorDet: det
          });
        }
      }
    }
  }

  return expanded;
}

function applySuperposedOppositeOpacity(inputSites, dimOpacity = 0.1) {
  const groups = new Map();
  const output = inputSites.map((site) => ({ ...site }));

  function rounded(values) {
    return values.map((value) => Number(value).toFixed(9)).join("|");
  }

  output.forEach((site, index) => {
    const canonical = canonicalOrbital(site.type) || site.type;
    const key = [
      rounded(site.pos),
      canonical,
      rounded([site.rotX, site.rotY, site.rotZ]),
      String(site.multipoleKind || "electric").toLowerCase()
    ].join("::");

    if (!groups.has(key)) {
      groups.set(key, []);
    }
    groups.get(key).push(index);
  });

  for (const indices of groups.values()) {
    if (indices.length < 2) {
      continue;
    }

    let hasPositive = false;
    let hasNegative = false;

    for (const idx of indices) {
      const site = output[idx];
      const [sx, sy, sz] = site.mirrorSign || [1, 1, 1];
      const parity = orbitalReflectionParity(site.type, sx, sy, sz);
      const basePhase = toNumber(site.phase, 1);
      const det = toNumber(site.mirrorDet, 1);
      const magneticAxial = String(site.multipoleKind).toLowerCase() === "magnetic" ? det : 1;
      const effectivePhase = basePhase * parity * magneticAxial;

      if (effectivePhase > EPSILON) {
        hasPositive = true;
      } else if (effectivePhase < -EPSILON) {
        hasNegative = true;
      }

      if (hasPositive && hasNegative) {
        break;
      }
    }

    if (hasPositive && hasNegative) {
      for (const idx of indices) {
        output[idx].opacity = dimOpacity;
      }
    }
  }

  return output;
}

function computeSceneRanges(drawSites, orbitalExtent = 1, padding = 0.4, minHalfSpan = 2) {
  if (drawSites.length === 0) {
    return {
      x: [-minHalfSpan, minHalfSpan],
      y: [-minHalfSpan, minHalfSpan],
      z: [-minHalfSpan, minHalfSpan]
    };
  }

  const mins = [Infinity, Infinity, Infinity];
  const maxs = [-Infinity, -Infinity, -Infinity];

  for (const site of drawSites) {
    for (let k = 0; k < 3; k += 1) {
      mins[k] = Math.min(mins[k], site.pos[k] - orbitalExtent - padding);
      maxs[k] = Math.max(maxs[k], site.pos[k] + orbitalExtent + padding);
    }
  }

  const center = [
    0.5 * (mins[0] + maxs[0]),
    0.5 * (mins[1] + maxs[1]),
    0.5 * (mins[2] + maxs[2])
  ];

  const halfSpan = Math.max(
    0.5 * Math.max(maxs[0] - mins[0], maxs[1] - mins[1], maxs[2] - mins[2]),
    minHalfSpan
  );

  return {
    x: [center[0] - halfSpan, center[0] + halfSpan],
    y: [center[1] - halfSpan, center[1] + halfSpan],
    z: [center[2] - halfSpan, center[2] + halfSpan]
  };
}

function addMirrorPlaneSurfaces(data, ranges, mirrorPlanes, offsets) {
  const [xMin, xMax] = ranges.x;
  const [yMin, yMax] = ranges.y;
  const [zMin, zMax] = ranges.z;

  if (mirrorPlanes.includes("xy")) {
    const z0 = offsets.xy;
    data.push({
      type: "surface",
      x: [[xMin, xMax], [xMin, xMax]],
      y: [[yMin, yMin], [yMax, yMax]],
      z: [[z0, z0], [z0, z0]],
      surfacecolor: [[0, 0], [0, 0]],
      colorscale: [[0, "#2ca02c"], [1, "#2ca02c"]],
      showscale: false,
      opacity: 0.18,
      hoverinfo: "skip"
    });
  }

  if (mirrorPlanes.includes("xz")) {
    const y0 = offsets.xz;
    data.push({
      type: "surface",
      x: [[xMin, xMax], [xMin, xMax]],
      y: [[y0, y0], [y0, y0]],
      z: [[zMin, zMin], [zMax, zMax]],
      surfacecolor: [[0, 0], [0, 0]],
      colorscale: [[0, "#ff7f0e"], [1, "#ff7f0e"]],
      showscale: false,
      opacity: 0.18,
      hoverinfo: "skip"
    });
  }

  if (mirrorPlanes.includes("yz")) {
    const x0 = offsets.yz;
    data.push({
      type: "surface",
      x: [[x0, x0], [x0, x0]],
      y: [[yMin, yMax], [yMin, yMax]],
      z: [[zMin, zMin], [zMax, zMax]],
      surfacecolor: [[0, 0], [0, 0]],
      colorscale: [[0, "#1f77b4"], [1, "#1f77b4"]],
      showscale: false,
      opacity: 0.18,
      hoverinfo: "skip"
    });
  }
}

function makeGrid(resolution) {
  const uVals = [];
  const vVals = [];

  for (let i = 0; i < resolution; i += 1) {
    uVals.push((2 * Math.PI * i) / (resolution - 1));
    vVals.push((Math.PI * i) / (resolution - 1));
  }

  return { uVals, vVals };
}

const GRID = makeGrid(RESOLUTION);

function buildSiteTraces(site) {
  const canonical = canonicalOrbital(site.type);
  if (!canonical) {
    return [];
  }

  const [sx, sy, sz] = site.mirrorSign || [1, 1, 1];
  const mirrorDet = toNumber(site.mirrorDet, 1);
  const mirrorPhase = String(site.multipoleKind || "electric").toLowerCase() === "magnetic" ? mirrorDet : 1;
  const phase = toNumber(site.phase, 1);

  const X = [];
  const Y = [];
  const Z = [];
  const R = [];

  let maxAbsR = 0;

  for (let i = 0; i < GRID.vVals.length; i += 1) {
    const v = GRID.vVals[i];
    const rowX = [];
    const rowY = [];
    const rowZ = [];
    const rowR = [];

    for (let j = 0; j < GRID.uVals.length; j += 1) {
      const u = GRID.uVals[j];
      const ux = Math.sin(v) * Math.cos(u);
      const uy = Math.sin(v) * Math.sin(u);
      const uz = Math.cos(v);

      const raw = phase * mirrorPhase * orbitalValue(canonical, ux, uy, uz);
      rowR.push(raw);
      maxAbsR = Math.max(maxAbsR, Math.abs(raw));

      rowX.push(ux);
      rowY.push(uy);
      rowZ.push(uz);
    }

    X.push(rowX);
    Y.push(rowY);
    Z.push(rowZ);
    R.push(rowR);
  }

  if (maxAbsR < EPSILON) {
    return [];
  }

  const xPos = [];
  const yPos = [];
  const zPos = [];
  const xNeg = [];
  const yNeg = [];
  const zNeg = [];

  for (let i = 0; i < GRID.vVals.length; i += 1) {
    const rowXPos = [];
    const rowYPos = [];
    const rowZPos = [];
    const rowXNeg = [];
    const rowYNeg = [];
    const rowZNeg = [];

    for (let j = 0; j < GRID.uVals.length; j += 1) {
      const normR = R[i][j] / maxAbsR;
      const absR = Math.abs(normR);

      const bx = absR * X[i][j];
      const by = absR * Y[i][j];
      const bz = absR * Z[i][j];

      const rotated = rotatePoint(bx, by, bz, toNumber(site.rotX, 0), toNumber(site.rotY, 0), toNumber(site.rotZ, 0));
      const finalX = sx * rotated[0] + site.pos[0];
      const finalY = sy * rotated[1] + site.pos[1];
      const finalZ = sz * rotated[2] + site.pos[2];

      if (normR >= 0) {
        rowXPos.push(finalX);
        rowYPos.push(finalY);
        rowZPos.push(finalZ);
        rowXNeg.push(NaN);
        rowYNeg.push(NaN);
        rowZNeg.push(NaN);
      } else {
        rowXPos.push(NaN);
        rowYPos.push(NaN);
        rowZPos.push(NaN);
        rowXNeg.push(finalX);
        rowYNeg.push(finalY);
        rowZNeg.push(finalZ);
      }
    }

    xPos.push(rowXPos);
    yPos.push(rowYPos);
    zPos.push(rowZPos);
    xNeg.push(rowXNeg);
    yNeg.push(rowYNeg);
    zNeg.push(rowZNeg);
  }

  const opacity = toNumber(site.opacity, 0.85);

  return [
    {
      type: "surface",
      x: xPos,
      y: yPos,
      z: zPos,
      surfacecolor: xPos.map((row) => row.map(() => 0)),
      colorscale: [[0, site.colorPos || "red"], [1, site.colorPos || "red"]],
      showscale: false,
      opacity,
      hoverinfo: "skip"
    },
    {
      type: "surface",
      x: xNeg,
      y: yNeg,
      z: zNeg,
      surfacecolor: xNeg.map((row) => row.map(() => 0)),
      colorscale: [[0, site.colorNeg || "blue"], [1, site.colorNeg || "blue"]],
      showscale: false,
      opacity,
      hoverinfo: "skip"
    }
  ];
}

function currentInputSite() {
  return {
    pos: [toNumber(dom.xPos.value, 0), toNumber(dom.yPos.value, 0), toNumber(dom.zPos.value, 0)],
    type: dom.orbitalType.value,
    phase: toNumber(dom.phase.value, 1),
    rotX: toNumber(dom.rotX.value, 0),
    rotY: toNumber(dom.rotY.value, 0),
    rotZ: toNumber(dom.rotZ.value, 0),
    multipoleKind: dom.multipoleKind.value || "electric",
    colorPos: dom.colorPos.value || "red",
    colorNeg: dom.colorNeg.value || "blue"
  };
}

function refreshSiteSelector() {
  const selected = Number(dom.siteSelect.value);
  dom.siteSelect.innerHTML = "";

  for (let i = 0; i < sites.length; i += 1) {
    const site = sites[i];
    const option = document.createElement("option");
    const kind = String(site.multipoleKind || "electric");
    option.value = String(i);
    option.textContent = `${i}: ${kind[0].toUpperCase()}-${site.type} @ (${site.pos[0].toFixed(2)}, ${site.pos[1].toFixed(2)}, ${site.pos[2].toFixed(2)}) r=(${site.rotX.toFixed(0)}, ${site.rotY.toFixed(0)}, ${site.rotZ.toFixed(0)})`;
    dom.siteSelect.appendChild(option);
  }

  if (sites.length === 0) {
    dom.siteSelect.disabled = true;
    return;
  }

  dom.siteSelect.disabled = false;
  if (Number.isFinite(selected) && selected >= 0 && selected < sites.length) {
    dom.siteSelect.value = String(selected);
  } else {
    dom.siteSelect.value = String(sites.length - 1);
  }
}

function status(message) {
  dom.status.textContent = message;
}

function plot() {
  const mirrorPlanes = selectedMirrorPlanes();
  const offsets = mirrorOffsets();

  let drawSites = sites.map((site) => ({ ...site }));

  if (dom.previewToggle.checked) {
    drawSites.push({ ...currentInputSite(), opacity: 0.45 });
  }

  drawSites = expandSitesWithMirrorPlanes(drawSites, mirrorPlanes, offsets);
  drawSites = applySuperposedOppositeOpacity(drawSites, 0.1);

  const ranges = computeSceneRanges(drawSites);
  const data = [];

  if (drawSites.length === 0) {
    addMirrorPlaneSurfaces(data, ranges, mirrorPlanes, offsets);

    Plotly.react(dom.graph, data, {
      scene: {
        xaxis: { title: "X (A)", range: ranges.x },
        yaxis: { title: "Y (A)", range: ranges.y },
        zaxis: { title: "Z (A)", range: ranges.z },
        aspectmode: "cube"
      },
      margin: { l: 0, r: 0, b: 0, t: 30 },
      annotations: [
        {
          text: "Add an orbital from the controls to begin.",
          xref: "paper",
          yref: "paper",
          x: 0.01,
          y: 0.99,
          showarrow: false
        }
      ]
    }, { responsive: true });

    return;
  }

  for (const site of drawSites) {
    const traces = buildSiteTraces(site);
    data.push(...traces);
  }

  addMirrorPlaneSurfaces(data, ranges, mirrorPlanes, offsets);

  Plotly.react(dom.graph, data, {
    scene: {
      xaxis: { title: "X (A)", range: ranges.x },
      yaxis: { title: "Y (A)", range: ranges.y },
      zaxis: { title: "Z (A)", range: ranges.z },
      aspectmode: "cube"
    },
    margin: { l: 0, r: 0, b: 0, t: 30 }
  }, { responsive: true });
}

function wireControls() {
  dom.addBtn.addEventListener("click", () => {
    const candidate = currentInputSite();
    if (!canonicalOrbital(candidate.type)) {
      status(`Unsupported orbital: ${candidate.type}`);
      return;
    }
    sites.push(candidate);
    refreshSiteSelector();
    plot();
    status(`Added ${candidate.type} at (${candidate.pos[0]}, ${candidate.pos[1]}, ${candidate.pos[2]}). Total: ${sites.length}`);
  });

  dom.swapBtn.addEventListener("click", () => {
    const prevPos = dom.colorPos.value;
    dom.colorPos.value = dom.colorNeg.value || "blue";
    dom.colorNeg.value = prevPos || "red";

    const idx = Number(dom.siteSelect.value);
    if (Number.isFinite(idx) && idx >= 0 && idx < sites.length) {
      const temp = sites[idx].colorPos;
      sites[idx].colorPos = sites[idx].colorNeg;
      sites[idx].colorNeg = temp;
    }

    plot();
    status("Swapped +/- colors.");
  });

  dom.undoBtn.addEventListener("click", () => {
    if (sites.length === 0) {
      status("Nothing to undo.");
      return;
    }
    sites.pop();
    refreshSiteSelector();
    plot();
    status(`Removed last orbital. Total: ${sites.length}`);
  });

  dom.trashBtn.addEventListener("click", () => {
    if (sites.length === 0) {
      status("Nothing to trash.");
      return;
    }

    const idx = Number(dom.siteSelect.value);
    if (!Number.isFinite(idx) || idx < 0 || idx >= sites.length) {
      status("Select an orbital to trash.");
      return;
    }

    const removed = sites[idx];
    sites = sites.filter((_, i) => i !== idx);
    refreshSiteSelector();
    plot();
    status(`Trashed ${removed.type} at (${removed.pos[0]}, ${removed.pos[1]}, ${removed.pos[2]}). Total: ${sites.length}`);
  });

  dom.clearBtn.addEventListener("click", () => {
    sites = [];
    refreshSiteSelector();
    plot();
    status("Cleared all orbitals.");
  });

  const reactiveInputs = [
    dom.orbitalType,
    dom.xPos,
    dom.yPos,
    dom.zPos,
    dom.rotX,
    dom.rotY,
    dom.rotZ,
    dom.colorPos,
    dom.colorNeg,
    dom.phase,
    dom.multipoleKind,
    dom.previewToggle,
    dom.mirrorXY,
    dom.mirrorXZ,
    dom.mirrorYZ,
    dom.mirrorXyZ,
    dom.mirrorXzY,
    dom.mirrorYzX,
    dom.siteSelect
  ];

  for (const input of reactiveInputs) {
    input.addEventListener("input", plot);
    input.addEventListener("change", plot);
  }
}

function initOrbitalOptions() {
  const orbitalKeys = Array.from(new Set(Object.keys(ORBITAL_ALIASES))).sort();

  for (const key of orbitalKeys) {
    const option = document.createElement("option");
    option.value = key;
    option.textContent = key;
    dom.orbitalType.appendChild(option);
  }

  dom.orbitalType.value = "dz2";
}

function init() {
  initOrbitalOptions();
  refreshSiteSelector();
  wireControls();
  plot();
}

init();
