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
          data-bs='${JSON.stringify(l.binding_site_residues)}'>
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
      openViewer(pdbid, ligands, bs);
    });
  });

  return card;
}

// ── 3D Viewer (NGL) ───────────────────────────────────────

let nglStage = null;
let currentComponent = null;
let isLoadingStructure = false;
let nglReprs = {};  // { protein: repr, ligand: repr, bs: repr }

const nglContainer      = document.getElementById("nglContainer");
const viewerPdbInput    = document.getElementById("viewerPdbInput");
const loadStructureBtn  = document.getElementById("loadStructureBtn");
const viewerLegend      = document.getElementById("viewerLegend");
const viewerLigandList  = document.getElementById("viewerLigandList");
const viewerResiduePicker = document.getElementById("viewerResiduePicker");

function pinContainerSize() {
  // NGL's WebGL canvas needs the container to report a real, non-zero pixel
  // size. Relying purely on CSS percentages inside a Grid cell has proven
  // unreliable here — the canvas ends up clamped to a few px. So we
  // explicitly read the parent grid cell's box and set hard pixel
  // dimensions on #nglContainer itself.
  const rect = nglContainer.getBoundingClientRect();
  const w = Math.max(Math.round(rect.width), 500);
  const h = Math.max(Math.round(rect.height), 400);
  nglContainer.style.width  = w + "px";
  nglContainer.style.height = h + "px";
  return { w, h };
}

function initNGL() {
  if (nglStage) return;
  // Remove placeholder
  nglContainer.innerHTML = "";

  const { w, h } = pinContainerSize();
  console.log(`[VHELIBS] nglContainer pinned to: ${w}x${h}`);

  nglStage = new NGL.Stage(nglContainer, {
    backgroundColor: "#1a1d27",
  });

  // Re-pin and resize on window resize, and once more shortly after
  // creation in case fonts/layout shifted the grid track size.
  window.addEventListener("resize", () => {
    pinContainerSize();
    nglStage.handleResize();
  });

  setTimeout(() => {
    const sz = pinContainerSize();
    nglStage.handleResize();
    console.log(`[VHELIBS] nglContainer re-pinned after settle: ${sz.w}x${sz.h}`);
  }, 300);
}

loadStructureBtn.addEventListener("click", () => {
  if (isLoadingStructure) return;
  const pdbid = viewerPdbInput.value.trim().toUpperCase();
  if (!pdbid) return;
  loadNGLStructure(pdbid, [], []);
});

function openViewer(pdbid, ligandRes, bsRes) {
  showTab("viewer");
  viewerPdbInput.value = pdbid.toUpperCase();
  loadNGLStructure(pdbid, ligandRes, bsRes);
}

function loadNGLStructure(pdbid, ligandRes, bsRes) {
  if (isLoadingStructure) {
    console.warn("loadNGLStructure called while a load is already in progress — ignoring.");
    return;
  }
  isLoadingStructure = true;
  loadStructureBtn.disabled = true;
  loadStructureBtn.textContent = "Loading…";

  initNGL();
  if (currentComponent) {
    nglStage.removeAllComponents();
    currentComponent = null;
  }
  viewerLegend.classList.add("hidden");
  viewerResiduePicker.classList.add("hidden");
  nglReprs = {};
  document.getElementById("chkProtein").checked = true;
  document.getElementById("chkLigand").checked = true;
  document.getElementById("chkBS").checked = true;

  const idLower = pdbid.toLowerCase();
  const idUpper = pdbid.toUpperCase();

  // Known-working RCSB download endpoints, in order of preference.
  // Each entry gets a unique "name" so NGL/the CIF parser never reuses an
  // internal cache key across different attempts.
  const sources = [
    { url: `https://files.rcsb.org/download/${idLower}.cif`,  ext: "cif", name: `${idLower}-cif` },
    { url: `https://files.rcsb.org/download/${idUpper}.pdb`,  ext: "pdb", name: `${idLower}-pdb` },
  ];

  tryLoadSource(sources, 0, pdbid, ligandRes, bsRes);
}

function finishLoading() {
  isLoadingStructure = false;
  loadStructureBtn.disabled = false;
  loadStructureBtn.textContent = "Load Structure";
}

function tryLoadSource(sources, index, pdbid, ligandRes, bsRes) {
  if (index >= sources.length) {
    finishLoading();
    alert(`Could not load structure ${pdbid}.\nNone of the RCSB endpoints responded with a valid structure file. Check your internet connection or try a different PDB ID.`);
    return;
  }

  const { url, ext, name } = sources[index];
  console.log(`[VHELIBS] Attempting load #${index}: ${url}`);
  nglStage.loadFile(url, { ext: ext, name: name, defaultRepresentation: false })
    .then(comp => {
      const n = comp && comp.structure ? comp.structure.atomCount : 0;
      console.log(`[VHELIBS] Load succeeded for ${url}, atomCount=${n}`);
      if (!comp || !comp.structure || !n) {
        if (comp) nglStage.removeComponent(comp);
        throw new Error("Empty structure returned");
      }
      currentComponent = comp;
      try {
        renderStructure(comp, ligandRes, bsRes);
      } catch (renderErr) {
        // A bug in our own representation/selection code should NOT trigger
        // a retry against a different file format — that would just load
        // the same (good) structure again. Surface it directly instead.
        console.error("[VHELIBS] renderStructure failed:", renderErr);
        alert("Structure loaded, but rendering failed: " + renderErr.message);
      }
      finishLoading();
    })
    .catch(err => {
      console.warn(`[VHELIBS] Load failed for ${url} (${ext}):`, err && err.message ? err.message : err);
      tryLoadSource(sources, index + 1, pdbid, ligandRes, bsRes);
    });
}

function renderStructure(comp, ligandRes, bsRes) {
  // Protein cartoon (muted)
  const proteinRepr = comp.addRepresentation("cartoon", {
    sele: "protein",
    color: "#64748b",
    opacity: 0.6,
    name: "protein"
  });

  // Ligand residues (amber ball+stick)
  let ligandRepr;
  if (ligandRes.length) {
    const ligSele = ligandRes.map(nglSele).join(" or ");
    ligandRepr = comp.addRepresentation("ball+stick", {
      sele: ligSele,
      color: "#f59e0b",
      multipleBond: "symmetric",
      name: "ligand"
    });
  } else {
    // Fallback: show all HETATMs
    ligandRepr = comp.addRepresentation("ball+stick", {
      sele: "hetero and not water",
      color: "#f59e0b",
      name: "ligand"
    });
  }

  // Binding site (blue licorice)
  let bsRepr = null;
  if (bsRes.length) {
    const bsSele = bsRes.map(nglSele).join(" or ");
    bsRepr = comp.addRepresentation("licorice", {
      sele: bsSele,
      color: "#5b7cf6",
      opacity: 0.85,
      name: "bs"
    });
  }

  // Guardar referencias
  nglReprs.protein = proteinRepr;
  nglReprs.ligand = ligandRepr;
  nglReprs.bs = bsRepr;

  // Configurar checkboxes
  const chkProtein = document.getElementById("chkProtein");
  const chkLigand = document.getElementById("chkLigand");
  const chkBS = document.getElementById("chkBS");

  chkProtein.onchange = function() {
    if (nglReprs.protein) nglReprs.protein.setVisibility(this.checked);
  };
  chkLigand.onchange = function() {
    if (nglReprs.ligand) nglReprs.ligand.setVisibility(this.checked);
  };
  chkBS.onchange = function() {
    if (nglReprs.bs) nglReprs.bs.setVisibility(this.checked);
  };
  // Forzar estado inicial
  chkProtein.onchange.call(chkProtein);
  chkLigand.onchange.call(chkLigand);
  chkBS.onchange.call(chkBS);

  comp.autoView();
  viewerLegend.classList.remove("hidden");

  // Re-fit the camera shortly after, in case the canvas was 0x0 (hidden tab,
  // layout not yet settled) at the time of the first autoView() call above.
  setTimeout(() => {
    if (nglStage) {
      pinContainerSize();
      nglStage.handleResize();
      comp.autoView();
    }
  }, 350);

  // Populate residue picker
  if (ligandRes.length) {
    viewerLigandList.innerHTML = "";
    ligandRes.forEach(res => {
      const btn = document.createElement("button");
      btn.className = "residue-btn";
      btn.textContent = res;
      btn.addEventListener("click", () => {
        const s = nglSele(res);
        comp.autoView(s, 1000);
      });
      viewerLigandList.appendChild(btn);
    });
    viewerResiduePicker.classList.remove("hidden");
  }
}

/**
 * Convert a VHELIBS residue string like "ATP A  42" to an NGL selection string.
 * Format: "RES CHAIN RESNUM" with fixed-width padding.
 */
function nglSele(res) {
  // res is like "ATP A  42" or "GLY A   5"
  res = res.trim();
  const parts = res.split(/\s+/);
  if (parts.length >= 3) {
    const resname = parts[0];
    const chain   = parts[1];
    const resnum  = parts[2];
    return `(${resnum} and :${chain} and [${resname}])`;
  }
  if (parts.length === 2) {
    const chain  = parts[0];
    const resnum = parts[1];
    return `(${resnum} and :${chain})`;
  }
  return res;
}

// ── Utilities ─────────────────────────────────────────────

function esc(s) {
  return String(s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
