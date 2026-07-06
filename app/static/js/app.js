/* =========================================================
   VHELIBS Web – frontend logic
   ========================================================= */

"use strict";

// ── Tab navigation ────────────────────────────────────────

const tabLinks = document.querySelectorAll(".nav-link[data-tab]");
const tabPanels = document.querySelectorAll(".tab-panel");

function showTab(name) {
  tabPanels.forEach(p => p.classList.toggle("active", p.id === name));
  tabLinks.forEach(l => l.classList.toggle("active", l.dataset.tab === name));

  // Mol*'s canvas can end up with a stale/zero size if it was created (or
  // last resized) while this tab was hidden (display:none). Force a re-pin
  // + resize whenever the viewer tab becomes visible, not just on load.
  if (name === "viewer" && viewerInstance) {
    requestAnimationFrame(() => {
      pinContainerSize();
      viewerInstance.handleResize();
    });
    setTimeout(() => {
      pinContainerSize();
      viewerInstance.handleResize();
    }, 300);
  }
}

tabLinks.forEach(link => {
  link.addEventListener("click", e => {
    e.preventDefault();
    showTab(link.dataset.tab);
  });
});

// ── Analysis submission ───────────────────────────────────

const analyseBtn    = document.getElementById("analyseBtn");
const pdbInput      = document.getElementById("pdbInput");
const progressArea  = document.getElementById("progressArea");
const progressBar   = document.getElementById("progressBar");
const progressLabel = document.getElementById("progressLabel");

analyseBtn.addEventListener("click", startAnalysis);
pdbInput.addEventListener("keydown", e => {
  if (e.key === "Enter" && e.ctrlKey) startAnalysis();
});

// ── Load PDB IDs from a file ──────────────────────────────
// Accepts plain text/CSV/TSV files: any mix of commas, whitespace or
// newlines works, same as manual textarea input, since it's just merged
// into pdbInput.value and sent to /api/analyse unchanged — the backend
// (and gatherConfig above) already knows how to split that string.
const pdbFileInput = document.getElementById("pdbFileInput");
const pdbFileName  = document.getElementById("pdbFileName");

pdbFileInput.addEventListener("change", () => {
  const file = pdbFileInput.files && pdbFileInput.files[0];
  if (!file) return;

  const MAX_SIZE = 1024 * 1024; // 1 MB is more than enough for a list of PDB IDs
  if (file.size > MAX_SIZE) {
    alert("File is too large. Please provide a plain text file with PDB IDs.");
    pdbFileInput.value = "";
    return;
  }

  pdbFileName.textContent = file.name;
  const reader = new FileReader();
  reader.onload = () => {
    const text = String(reader.result || "").trim();
    if (!text) {
      alert(`"${file.name}" appears to be empty.`);
      return;
    }
    const existing = pdbInput.value.trim();
    pdbInput.value = existing ? existing + "\n" + text : text;
  };
  reader.onerror = () => alert(`Could not read file "${file.name}".`);
  reader.readAsText(file);

  // Reset so re-selecting the same file still fires 'change'.
  pdbFileInput.value = "";
});

function gatherConfig() {
  const v = id => document.getElementById(id).value;
  const b = id => document.getElementById(id).checked;
  return {
    pdbids:          pdbInput.value,
    rsr_upper:       parseFloat(v("th_rsr_upper")),
    rsr_lower:       parseFloat(v("th_rsr_lower")),
    rscc_min:        parseFloat(v("th_rscc_min")),
    rfree_max:       parseFloat(v("th_rfree_max")),
    occupancy_min:   parseFloat(v("th_occupancy_min")),
    tolerance:       parseInt(v("th_tolerance")),
    distance:        parseFloat(v("th_distance")),
    use_pdb_redo:    b("usePdbRedo"),
    check_owab:      b("chk_owab"),
    owab_max:        parseFloat(v("th_owab_max")),
    check_resolution: b("chk_resolution"),
    resolution_max:  parseFloat(v("th_resolution_max")),
    use_rdiff:       b("chk_rdiff"),
    rdiff_max:       parseFloat(v("th_rdiff_max")),
    use_dpi:         b("chk_dpi"),
    dpi_max:         parseFloat(v("th_dpi_max")),
  };
}

function startAnalysis() {
  const cfg = gatherConfig();
  const ids = cfg.pdbids.trim();
  if (!ids) { alert("Please enter at least one PDB ID."); return; }

  analyseBtn.disabled = true;
  analyseBtn.textContent = "Analysing…";
  progressArea.classList.remove("hidden");
  progressBar.style.width = "0%";
  progressLabel.textContent = "Submitting job…";

  clearResults();

  fetch("/api/analyse", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(cfg),
  })
  .then(r => r.json())
  .then(data => {
    if (data.error) {
      alert("Error: " + data.error);
      resetAnalyseBtn();
      return;
    }
    pollJob(data.job_id, data.total);
  })
  .catch(err => {
    alert("Network error: " + err.message);
    resetAnalyseBtn();
  });
}

function pollJob(jobId, total) {
  const POLL_MS = 1500;
  const poll = () => {
    fetch("/api/status/" + jobId)
      .then(r => r.json())
      .then(data => {
        const done = data.progress || 0;
        const pct  = total > 0 ? Math.round((done / total) * 100) : 0;
        progressBar.style.width = pct + "%";
        progressLabel.textContent = `Processed ${done} / ${total} structure(s)…`;

        if (data.status === "done") {
          progressBar.style.width = "100%";
          progressLabel.textContent = "Analysis complete.";
          resetAnalyseBtn();
          renderResults(data.results);
          showTab("results");
        } else {
          setTimeout(poll, POLL_MS);
        }
      })
      .catch(() => setTimeout(poll, POLL_MS * 2));
  };
  setTimeout(poll, POLL_MS);
}

function resetAnalyseBtn() {
  analyseBtn.disabled = false;
  analyseBtn.textContent = "Analyse";
}

// ── Results rendering ─────────────────────────────────────

const resultsEmpty     = document.getElementById("resultsEmpty");
const resultsSummary   = document.getElementById("resultsSummary");
const resultsContainer = document.getElementById("resultsContainer");

function clearResults() {
  resultsEmpty.classList.remove("hidden");
  resultsSummary.classList.add("hidden");
  resultsContainer.classList.add("hidden");
  resultsSummary.innerHTML = "";
  resultsContainer.innerHTML = "";
}

function qualityClass(q) {
  if (!q || q === "Good")    return "good";
  if (q === "Dubious")       return "dubious";
  if (q === "Bad")           return "bad";
  return "error";
}

function renderResults(results) {
  if (!results || !results.length) {
    resultsEmpty.classList.remove("hidden");
    return;
  }
  resultsEmpty.classList.add("hidden");

  // Tally summary
  let counts = { good: 0, dubious: 0, bad: 0, error: 0 };
  results.forEach(r => {
    if (r.error) { counts.error++; return; }
    (r.ligands || []).forEach(l => {
      counts[qualityClass(l.ligand_quality)]++;
    });
  });

  resultsSummary.innerHTML = `
    <span class="summary-pill pill-good">✓ Good ${counts.good}</span>
    <span class="summary-pill pill-dubious">~ Dubious ${counts.dubious}</span>
    <span class="summary-pill pill-bad">✕ Bad ${counts.bad}</span>
    ${counts.error ? `<span class="summary-pill pill-error">⚠ Errors ${counts.error}</span>` : ""}
  `;
  resultsSummary.classList.remove("hidden");

  resultsContainer.innerHTML = "";
  results.forEach(r => {
    resultsContainer.appendChild(buildResultCard(r));
  });
  resultsContainer.classList.remove("hidden");
}

function buildResultCard(r) {
  if (r.error) {
    const div = document.createElement("div");
    div.className = "error-card";
    div.innerHTML = `
      <span class="error-pdbid">${esc(r.pdbid.toUpperCase())}</span>
      <span class="error-reason">${esc(r.error)}</span>
    `;
    return div;
  }

  const card = document.createElement("div");
  card.className = "result-card";

  // Overall quality: worst ligand or "No ligands"
  const ligands = r.ligands || [];
  let overallQ = ligands.length ? "Good" : null;
  ligands.forEach(l => {
    if (l.ligand_quality === "Bad") overallQ = "Bad";
    else if (l.ligand_quality === "Dubious" && overallQ !== "Bad") overallQ = "Dubious";
  });
  const qc = qualityClass(overallQ);

  const header = document.createElement("div");
  header.className = "result-card-header";
  header.innerHTML = `
    <span class="result-pdbid">${esc(r.pdbid.toUpperCase())}</span>
    <span class="result-badges">
      ${ligands.length} ligand(s)
      ${overallQ ? `<span class="badge badge-${qc}">${overallQ}</span>` : ""}
    </span>
  `;
  card.appendChild(header);

  const body = document.createElement("div");
  body.className = "result-card-body";

  // Structure stats
  const sd = r.struc_dict || {};
  const fmt = v => (v != null ? parseFloat(v).toFixed(4) : "N/A");
  body.innerHTML = `
    <div class="struc-stats">
      <span class="stat-chip"><strong>R-free</strong> ${fmt(sd.rFree)}</span>
      <span class="stat-chip"><strong>R-work</strong> ${fmt(sd.rWork)}</span>
      ${sd.Resolution != null ? `<span class="stat-chip"><strong>Resolution</strong> ${fmt(sd.Resolution)} Å</span>` : ""}
    </div>
  `;

  // Ligand entries
  const ligDiv = document.createElement("div");
  ligDiv.className = "ligand-entries";
  ligands.forEach((l, i) => {
    const qc2 = qualityClass(l.ligand_quality);
    const bsQc = qualityClass(l.binding_site_quality);
    const entry = document.createElement("div");
    entry.className = "ligand-entry";
    entry.innerHTML = `
      <div class="ligand-entry-header">
        <div>
          <span class="badge badge-${qc2}">Ligand: ${l.ligand_quality}</span>
          <span class="badge badge-${bsQc}" style="margin-left:6px">BS: ${l.binding_site_quality}</span>
        </div>
        <button class="view-btn" data-pdbid="${esc(r.pdbid)}"
          data-ligands='${JSON.stringify(l.ligand_residues)}'
          data-bs='${JSON.stringify(l.binding_site_residues)}'
          data-rte='${JSON.stringify(l.residues_to_examine || [])}'
          data-boxes='${JSON.stringify(l.density_boxes || {})}'
          data-atoms='${JSON.stringify(l.density_atoms || {})}'>
          View 3D
        </button>
      </div>
      <div class="ligand-residues">
        <strong>Ligand:</strong>
        ${l.ligand_residues.map(s => `<span class="residue-tag">${esc(s)}</span>`).join("")}
      </div>
      ${l.binding_site_residues.length ? `
      <div class="ligand-residues" style="margin-top:6px">
        <strong>Binding site:</strong>
        ${l.binding_site_residues.slice(0, 8).map(s => `<span class="residue-tag">${esc(s)}</span>`).join("")}
        ${l.binding_site_residues.length > 8 ? `<span class="residue-tag">+${l.binding_site_residues.length-8} more</span>` : ""}
      </div>` : ""}
      ${l.low_occupancy && l.low_occupancy.length ? `
      <div class="ligand-residues" style="margin-top:6px;color:var(--clr-dubious)">
        ⚠ Low occupancy: ${l.low_occupancy.map(s => `<span class="residue-tag">${esc(s)}</span>`).join("")}
      </div>` : ""}
    `;
    ligDiv.appendChild(entry);
  });

  // Rejected residues
  const rejected = r.rejected || {};
  const rejKeys = Object.keys(rejected);
  if (rejKeys.length) {
    const rejDiv = document.createElement("details");
    rejDiv.style.marginTop = "12px";
    rejDiv.innerHTML = `<summary style="cursor:pointer;font-size:.8125rem;color:var(--clr-muted)">${rejKeys.length} rejected residue(s)</summary>
      <div style="margin-top:6px;font-size:.8rem;font-family:monospace;color:var(--clr-muted)">
        ${rejKeys.map(k => `<div><span class="residue-tag">${esc(k)}</span> ${esc(rejected[k])}</div>`).join("")}
      </div>`;
    ligDiv.appendChild(rejDiv);
  }

  body.appendChild(ligDiv);
  card.appendChild(body);

  header.addEventListener("click", () => {
    const open = body.classList.toggle("open");
    header.classList.toggle("open", open);
  });

  // Wire "View 3D" buttons
  card.querySelectorAll(".view-btn").forEach(btn => {
    btn.addEventListener("click", e => {
      e.stopPropagation();
      const pdbid   = btn.dataset.pdbid;
      const ligands = JSON.parse(btn.dataset.ligands || "[]");
      const bs      = JSON.parse(btn.dataset.bs || "[]");
      const rte     = JSON.parse(btn.dataset.rte || "[]");
      const boxes   = JSON.parse(btn.dataset.boxes || "{}");
      const atoms   = JSON.parse(btn.dataset.atoms || "{}");
      openViewer(pdbid, ligands, bs, rte, boxes, atoms);
    });
  });

  return card;
}

// ── 3D Viewer (Mol*) ──────────────────────────────────────
//
// Migrated from NGL to Mol*. Instead of imperatively adding/removing
// representations on a live component (the NGL approach), the scene is
// described declaratively as a MolViewSpec (MVS) tree and (re)loaded via
// the MolViewSpec extension. This has a useful side effect for us: Mol*
// diffs the new state tree against the current one and only recomputes
// the parts that actually changed, so re-loading the MVS tree on every
// checkbox toggle does NOT re-download/re-parse the structure file — only
// the affected representations are touched. Camera focus/jump-to-residue
// is handled separately via viewer.structureInteractivity(), which doesn't
// require touching the state tree at all.

let viewerInitPromise = null;   // Promise<Viewer>, created once
let viewerInstance    = null;   // resolved Mol* Viewer (has .plugin)
let isLoadingStructure = false;
let currentPdbId      = null;
let currentLigandRes  = [];
let currentBsRes      = [];
let currentRteRes     = [];
let currentDensityBoxes = null; // {ligand, binding_site, residues_to_examine} bboxes, or null if unknown
let currentDensityAtoms = null; // {ligand, binding_site, residues_to_examine} per-atom centers, or null
let currentAtomRadius = 1.6; // Å, radius of the per-atom density clip sphere (user-adjustable)
let currentFocusRes = null; // residue string clicked in the "residues to examine" list, or null
const layerState = { protein: true, ligand: true, bs: true }; // structure checkbox state
const densityLayerState = { ligand: true, bs: true, rte: true }; // density checkbox state
let currentIsovalue = 1.0; // relative sigma units

const molContainer        = document.getElementById("nglContainer"); // id kept for CSS compat
const viewerPdbInput      = document.getElementById("viewerPdbInput");
const loadStructureBtn    = document.getElementById("loadStructureBtn");
const viewerLigandList    = document.getElementById("viewerLigandList");
const viewerResiduePicker = document.getElementById("viewerResiduePicker");
const viewerDensityControls = document.getElementById("viewerDensityControls");
const isovalueSlider      = document.getElementById("isovalueSlider");
const isovalueReadout     = document.getElementById("isovalueReadout");
const atomRadiusSlider    = document.getElementById("atomRadiusSlider");
const atomRadiusReadout   = document.getElementById("atomRadiusReadout");

function pinContainerSize() {
  // Mol*'s WebGL canvas needs the container to report a real, non-zero
  // pixel size. Relying purely on CSS percentages inside a Grid cell has
  // proven unreliable here — the canvas ends up clamped to a few px. So we
  // explicitly read the parent grid cell's box and set hard pixel
  // dimensions on the container itself. Same approach as before, just
  // pointed at Mol* instead of NGL.
  const rect = molContainer.getBoundingClientRect();
  const w = Math.max(Math.round(rect.width), 500);
  const h = Math.max(Math.round(rect.height), 400);
  molContainer.style.width  = w + "px";
  molContainer.style.height = h + "px";
  return { w, h };
}

function initMolstar() {
  if (viewerInitPromise) return viewerInitPromise;
  molContainer.innerHTML = ""; // remove placeholder

  const { w, h } = pinContainerSize();
  console.log(`[VHELIBS] molContainer pinned to: ${w}x${h}`);

  viewerInitPromise = molstar.Viewer.create(molContainer, {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowLeftPanel: false,
    layoutShowRemoteState: false,
    viewportShowExpand: true,
    viewportShowControls: false,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false,
    viewportShowSettings: false,
    viewportBackgroundColor: "#1a1d27",
    pdbProvider: "rcsb",
  }).then(viewer => {
    viewerInstance = viewer;

    // Re-pin and resize on window resize, and once more shortly after
    // creation in case fonts/layout shifted the grid track size.
    window.addEventListener("resize", () => {
      pinContainerSize();
      viewerInstance.handleResize();
    });

    setTimeout(() => {
      const sz = pinContainerSize();
      viewerInstance.handleResize();
      console.log(`[VHELIBS] molContainer re-pinned after settle: ${sz.w}x${sz.h}`);
    }, 300);

    return viewer;
  });

  return viewerInitPromise;
}

loadStructureBtn.addEventListener("click", () => {
  if (isLoadingStructure) return;
  const pdbid = viewerPdbInput.value.trim().toUpperCase();
  if (!pdbid) return;
  loadMolstarStructure(pdbid, [], [], [], null, null);
});

function openViewer(pdbid, ligandRes, bsRes, rteRes, densityBoxes, densityAtoms) {
  showTab("viewer");
  viewerPdbInput.value = pdbid.toUpperCase();
  loadMolstarStructure(pdbid, ligandRes, bsRes, rteRes || [], densityBoxes || null, densityAtoms || null);
}

async function loadMolstarStructure(pdbid, ligandRes, bsRes, rteRes, densityBoxes, densityAtoms) {
  if (isLoadingStructure) {
    console.warn("loadMolstarStructure called while a load is already in progress — ignoring.");
    return;
  }
  isLoadingStructure = true;
  loadStructureBtn.disabled = true;
  loadStructureBtn.textContent = "Loading…";
  viewerResiduePicker.classList.add("hidden");
  viewerDensityControls.classList.add("hidden");

  currentPdbId = pdbid;
  currentLigandRes = ligandRes;
  currentBsRes = bsRes;
  currentRteRes = rteRes || [];
  currentDensityBoxes = densityBoxes || null;
  currentDensityAtoms = densityAtoms || null;
  currentAtomRadius = 1.6;
  currentFocusRes = null;
  layerState.protein = true;
  layerState.ligand = true;
  layerState.bs = true;
  densityLayerState.ligand = true;
  densityLayerState.bs = true;
  densityLayerState.rte = true;
  currentIsovalue = 1.0;
  document.getElementById("chkProtein").checked = true;
  document.getElementById("chkLigand").checked = true;
  document.getElementById("chkBS").checked = true;
  document.getElementById("chkDensityLigand").checked = true;
  document.getElementById("chkDensityBS").checked = true;
  document.getElementById("chkDensityRTE").checked = true;
  isovalueSlider.value = currentIsovalue;
  isovalueReadout.textContent = currentIsovalue.toFixed(1) + "σ";
  atomRadiusSlider.value = currentAtomRadius;
  atomRadiusReadout.textContent = currentAtomRadius.toFixed(1) + " Å";

  try {
    await initMolstar();
    pinContainerSize();
    await applyMolstarScene({ focus: true, keepCamera: false });
    populateResiduePicker(currentRteRes.length ? currentRteRes : ligandRes);
    // Only show density controls when we actually have region boxes to
    // query (i.e. came from an analysis result, not a bare PDB ID load).
    if (currentDensityBoxes && hasAnyBox(currentDensityBoxes)) {
      viewerDensityControls.classList.remove("hidden");
    }
  } catch (err) {
    console.error("[VHELIBS] Mol* load failed:", err);
    alert(`Could not load structure ${pdbid}.\nNone of the RCSB endpoints responded with a valid structure file. Check your internet connection or try a different PDB ID.`);
  } finally {
    finishLoading();
    // Re-fit shortly after, in case the canvas was 0x0 (hidden tab, layout
    // not yet settled) at the time this ran.
    setTimeout(() => {
      if (viewerInstance) {
        pinContainerSize();
        viewerInstance.handleResize();
      }
    }, 350);
  }
}

function finishLoading() {
  isLoadingStructure = false;
  loadStructureBtn.disabled = false;
  loadStructureBtn.textContent = "Load Structure";
}

// RCSB download endpoints to try, in order of preference.
function molstarSources(pdbidLower) {
  return [
    { url: `https://files.rcsb.org/download/${pdbidLower}.cif`, format: "mmcif" },
    { url: `https://files.rcsb.org/download/${pdbidLower}.pdb`, format: "pdb" },
  ];
}

/**
 * Build and load the MVS scene for the current PDB ID, honouring
 * layerState (which representations to include) and residue selections.
 * Tries each RCSB source in turn, the same fallback behaviour the NGL
 * version had.
 */
async function applyMolstarScene({ focus, keepCamera }) {
  const sources = molstarSources(currentPdbId.toLowerCase());
  let lastErr = null;

  for (const src of sources) {
    try {
      const mvsData = buildMvsData(
        src.url, src.format, currentLigandRes, currentBsRes, layerState, focus,
        currentDensityBoxes, currentDensityAtoms, currentAtomRadius, densityLayerState, currentIsovalue,
        currentFocusRes
      );
      await molstar.PluginExtensions.mvs.loadMVS(viewerInstance.plugin, mvsData, {
        replaceExisting: true,
        keepCamera: !!keepCamera,
        sanityChecks: true,
      });
      return;
    } catch (err) {
      console.warn(`[VHELIBS] Mol* load failed for ${src.url}:`, err && err.message ? err.message : err);
      lastErr = err;
    }
  }
    throw lastErr || new Error("Unknown Mol* load error");
  }

// Region -> {color, layerKey, boxKey} mapping shared by buildMvsData and
// the checkbox wiring below.
const DENSITY_REGIONS = [
  { boxKey: "ligand",               layerKey: "ligand", color: "#e879f9" }, // magenta
  { boxKey: "binding_site",         layerKey: "bs",      color: "#22d3ee" }, // cyan
  { boxKey: "residues_to_examine",  layerKey: "rte",     color: "#facc15" }, // yellow
];

function hasAnyBox(densityBoxes) {
  return DENSITY_REGIONS.some(r => densityBoxes && densityBoxes[r.boxKey]);
}

function buildMvsData(sourceUrl, sourceFormat, ligandRes, bsRes, layers, focus,
                       densityBoxes, densityAtoms, atomRadius, densityLayers, isovalue,
                       focusResidue) {
  const builder = molstar.PluginExtensions.mvs.MVSData.createBuilder();
  const structure = builder
    .download({ url: sourceUrl })
    .parse({ format: sourceFormat })
    .modelStructure({});

  // Protein cartoon (muted) — mirrors the old NGL "protein" representation.
  if (layers.protein) {
    structure.component({ selector: "protein" })
      .representation({ type: "cartoon" })
      .color({ color: "#64748b" })
      .opacity({ opacity: 0.6 });
  }

  // Ligand residues (amber ball+stick). Falls back to all non-water
  // hetero atoms when no explicit residue list is known, same as before.
  if (layers.ligand) {
    const ligandSelectors = ligandRes.map(residueSelector).filter(Boolean);
    const ligComp = ligandSelectors.length
      ? structure.component({ selector: ligandSelectors })
      : structure.component({ selector: "ligand" });
    ligComp.representation({ type: "ball_and_stick" }).color({ color: "#f59e0b" });
    // Auto-focus the camera on the ligand on first load, like the old
    // comp.autoView() call. MVS has no direct "licorice" type, so the
    // binding site below reuses ball_and_stick, which is visually close.
    // Skipped when the person has explicitly picked a residue to examine
    // below — that selection takes priority over the default ligand focus.
    if (focus && !focusResidue && ligandSelectors.length) {
      ligComp.focus({});
    }
  }

  // A residue clicked in the "Residues to examine" list: highlight it in
  // white and point the camera at it, overriding the default ligand focus.
  if (focusResidue) {
    const sel = residueSelector(focusResidue);
    if (sel) {
      const focusComp = structure.component({ selector: sel });
      focusComp.representation({ type: "ball_and_stick" })
        .color({ color: "#ffffff" })
        .opacity({ opacity: 1 });
      focusComp.focus({});
    }
  }

  // Binding site (blue, semi-transparent ball+stick).
  if (layers.bs && bsRes.length) {
    const bsSelectors = bsRes.map(residueSelector).filter(Boolean);
    if (bsSelectors.length) {
      structure.component({ selector: bsSelectors })
        .representation({ type: "ball_and_stick" })
        .color({ color: "#5b7cf6" })
        .opacity({ opacity: 0.85 });
    }
  }

  // ── Segmented electron density ────────────────────────────────────────
  // Each region (ligand / binding site / residues-to-examine) gets its
  // own density query, scoped to that region's padded bounding box (see
  // core.rsr_core.residues_bbox on the backend). This is the "masking"
  // step from the spec: rather than loading the whole map and clipping it
  // in the browser, we only ever ask EBI for the chunk of density that
  // covers the region we want to show, which is both cheaper and avoids
  // having to reimplement spatial masking client-side.
// Hard cap on clip representations per region — one Mol* representation
// node is created per atom, and a very large binding site (hundreds of
// atoms) would otherwise generate an excessive number of GPU draws. This
// keeps things responsive while still giving true per-atom masking for
// the vast majority of ligands/binding sites.
const MAX_CLIP_ATOMS_PER_REGION = 300;

  if (densityBoxes) {
    for (const region of DENSITY_REGIONS) {
      if (!densityLayers[region.layerKey]) continue;
      const box = densityBoxes[region.boxKey];
      if (!box || !box.min || !box.max) continue;

      // The box only sizes the *download window* — the EBI density server
      // has no per-atom masking, so the box for e.g. "binding site" also
      // contains the ligand and unrelated solvent. To actually show density
      // that differs per residue/atom, the same downloaded volume is
      // rendered as one representation per atom, each clipped to a small
      // sphere of `atomRadius` around that atom (see
      // rsr_core.residue_atom_centers). No extra network cost: it's one
      // download, just multiple clipped representations of it.
      const densityUrl = densityBoxUrl(currentPdbId, box);
      const volume = builder.download({ url: densityUrl }).parse({ format: "bcif" })
        .volume({ channel_id: "2FO-FC" });

      let atoms = (densityAtoms && densityAtoms[region.boxKey]) || [];
    if (atoms.length > MAX_CLIP_ATOMS_PER_REGION) {
      console.warn(`[VHELIBS] "${region.boxKey}" has ${atoms.length} atoms — capping per-atom density masking to ${MAX_CLIP_ATOMS_PER_REGION}.`);
      atoms = atoms.slice(0, MAX_CLIP_ATOMS_PER_REGION);
    }

    const firstRep = volume.representation({ type: "isosurface", relative_isovalue: isovalue, show_wireframe: false, show_faces: true });
    // volume-representation .clip() requires Mol* >= 5.0 (see index.html
    // comment). Detected once per representation rather than assumed, so
    // an older/unexpected Mol* build degrades gracefully instead of
    // throwing and aborting the whole structure load.
    //
    // IMPORTANT: a clip node's default behaviour is to cut its shape OUT
    // of the representation (hide what's inside, keep what's outside) —
    // like slicing through a structure to see inside it. We want the
    // opposite: keep only what's inside each atom's sphere and hide
    // everything else, so `invert: true` is required below, otherwise
    // each clip just punches an invisible pinhole out of the full box
    // and the whole box still renders.
    const supportsClip = atoms.length > 0 && typeof firstRep.clip === "function";

    if (supportsClip) {
      firstRep.color({ color: region.color }).opacity({ opacity: 0.35 })
        .clip({ type: "sphere", center: atoms[0].center, radius: atomRadius, invert: true });
      for (let i = 1; i < atoms.length; i++) {
        volume.representation({ type: "isosurface", relative_isovalue: isovalue, show_wireframe: false, show_faces: true })
          .color({ color: region.color })
          .opacity({ opacity: 0.35 })
          .clip({ type: "sphere", center: atoms[i].center, radius: atomRadius, invert: true });
      }
    } else {
      if (atoms.length) {
        console.warn(`[VHELIBS] This Mol* build has no volume clip() support (requires Mol* >= 5.0) — showing the unmasked region box for "${region.boxKey}" instead of per-atom density.`);
      }
      firstRep.color({ color: region.color }).opacity({ opacity: 0.35 });
    }
  }
  }

  return builder.getState();
}

// Builds the URL for a region's density chunk. Points at our own backend
// proxy (see edm_routes.py) by default, which caches repeated queries and
// keeps a single external-network boundary; EBI's endpoint could also be
// called directly from the browser if the proxy isn't deployed.
function densityBoxUrl(pdbid, box, detail = 3) {
  const min = box.min.join(",");
  const max = box.max.join(",");
  return `/api/density-box/${pdbid.toLowerCase()}?min=${min}&max=${max}&detail=${detail}`;
}

function populateResiduePicker(residues) {
  if (!residues.length) {
    viewerResiduePicker.classList.add("hidden");
    return;
  }
  viewerLigandList.innerHTML = "";
  residues.forEach(res => {
    const btn = document.createElement("button");
    btn.className = "residue-btn";
    btn.textContent = res;
    btn.dataset.residue = res;
    btn.addEventListener("click", () => {
      if (!residueSelector(res) || !currentPdbId || isLoadingStructure || !viewerInstance) return;
      currentFocusRes = res;
      viewerLigandList.querySelectorAll(".residue-btn").forEach(b => {
        const isActive = b.dataset.residue === res;
        b.classList.toggle("active", isActive);
        // Inline fallback in case style.css has no .residue-btn.active rule.
        b.style.outline = isActive ? "2px solid #fff" : "";
      });
      applyMolstarScene({ focus: false, keepCamera: false })
        .catch(err => console.error("[VHELIBS] Failed to focus residue:", err));
    });
    viewerLigandList.appendChild(btn);
  });
  viewerResiduePicker.classList.remove("hidden");
}

// Layer checkboxes: rebuild the MVS scene with the same structure/camera
// but only the checked representations included. Because Mol* diffs the
// state tree, this only touches the representations that actually change.
["chkProtein", "chkLigand", "chkBS"].forEach((id, i) => {
  const key = ["protein", "ligand", "bs"][i];
  document.getElementById(id).addEventListener("change", function () {
    layerState[key] = this.checked;
    if (!currentPdbId || isLoadingStructure || !viewerInstance) return;
    applyMolstarScene({ focus: false, keepCamera: true })
      .catch(err => console.error("[VHELIBS] Failed to update layer visibility:", err));
  });
});

// Density region checkboxes (ligand / binding site / residues to examine).
["chkDensityLigand", "chkDensityBS", "chkDensityRTE"].forEach((id, i) => {
  const key = ["ligand", "bs", "rte"][i];
  document.getElementById(id).addEventListener("change", function () {
    densityLayerState[key] = this.checked;
    if (!currentPdbId || isLoadingStructure || !viewerInstance) return;
    applyMolstarScene({ focus: false, keepCamera: true })
      .catch(err => console.error("[VHELIBS] Failed to update density visibility:", err));
  });
});

// Isovalue (contour level) slider — applies to all visible density
// volumes at once. Debounced slightly so dragging doesn't spam reloads.
let isovalueDebounce = null;
isovalueSlider.addEventListener("input", function () {
  currentIsovalue = parseFloat(this.value);
  isovalueReadout.textContent = currentIsovalue.toFixed(1) + "σ";
  if (!currentPdbId || isLoadingStructure || !viewerInstance) return;
  clearTimeout(isovalueDebounce);
  isovalueDebounce = setTimeout(() => {
    applyMolstarScene({ focus: false, keepCamera: true })
      .catch(err => console.error("[VHELIBS] Failed to update isovalue:", err));
  }, 150);
});

// Atom mask radius slider — controls the radius of the per-atom clip
// sphere used to mask density (see buildMvsData). Larger radius = more of
// the surrounding density bleeds in per atom; smaller = tighter to the
// atom itself. Debounced the same way as the isovalue slider.
let atomRadiusDebounce = null;
atomRadiusSlider.addEventListener("input", function () {
  currentAtomRadius = parseFloat(this.value);
  atomRadiusReadout.textContent = currentAtomRadius.toFixed(1) + " Å";
  if (!currentPdbId || isLoadingStructure || !viewerInstance) return;
  clearTimeout(atomRadiusDebounce);
  atomRadiusDebounce = setTimeout(() => {
    applyMolstarScene({ focus: false, keepCamera: true })
      .catch(err => console.error("[VHELIBS] Failed to update atom mask radius:", err));
  }, 150);
});

/**
 * Convert a VHELIBS residue string like "ATP A  42" to a Mol* MVS
 * ComponentExpression selector: { auth_asym_id, auth_seq_id }.
 * Format: "RESNAME CHAIN RESNUM" with fixed-width padding, same layout
 * used throughout the rest of the app (residue tags, binding-site lists).
 */
function residueSelector(res) {
  res = res.trim();
  const parts = res.split(/\s+/);
  if (parts.length >= 3) {
    const chain  = parts[1];
    const resnum = parseInt(parts[2], 10);
    if (Number.isNaN(resnum)) return null;
    return { auth_asym_id: chain, auth_seq_id: resnum };
  }
  if (parts.length === 2) {
    const chain  = parts[0];
    const resnum = parseInt(parts[1], 10);
    if (Number.isNaN(resnum)) return null;
    return { auth_asym_id: chain, auth_seq_id: resnum };
  }
  return null;
}

// ── Utilities ─────────────────────────────────────────────

function esc(s) {
  return String(s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
