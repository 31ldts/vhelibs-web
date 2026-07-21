/* =========================================================
   VHELIBS Web – frontend logic
   ========================================================= */

"use strict";

// ── Viewer layout sizing ──────────────────────────────────
// .viewer-layout (see style.css) sizes itself as
// calc(100vh - var(--header-h) - var(--footer-h) - 56px) so the mol*
// viewport gets a fixed, absolute height that doesn't depend on how much
// content is currently inside the sidebar (.viewer-controls). Both
// --header-h and --footer-h have CSS fallbacks, but header/footer height
// isn't hardcoded in HTML/CSS (the header can wrap to two lines on narrow
// viewports), so we measure them for real here and keep them in sync on
// resize.
const siteHeaderEl = document.querySelector(".site-header");
const siteFooterEl = document.querySelector(".site-footer");

function updateHeaderHeightVar() {
  if (siteHeaderEl) {
    document.documentElement.style.setProperty("--header-h", siteHeaderEl.getBoundingClientRect().height + "px");
  }
  if (siteFooterEl) {
    document.documentElement.style.setProperty("--footer-h", siteFooterEl.getBoundingClientRect().height + "px");
  }
}

updateHeaderHeightVar();
window.addEventListener("resize", updateHeaderHeightVar);

// ── Tab navigation ────────────────────────────────────────

const tabLinks = document.querySelectorAll(".nav-link[data-tab]");
const tabPanels = document.querySelectorAll(".tab-panel");

// Switching tabs hides the previously-active panel (display:none), which
// shrinks the document and makes the browser clamp/lose the window scroll
// position. To make each tab "remember" where the user was scrolled to
// (e.g. scrolled down in Results, hopping to the 3D Viewer and back), we
// save window.scrollY for the outgoing tab and restore it for the
// incoming one once the new panel's layout has settled.
let activeTabName = null;
const tabScrollPositions = Object.create(null);

function showTab(name) {
  if (activeTabName && activeTabName !== name) {
    tabScrollPositions[activeTabName] = window.scrollY;
  }

  tabPanels.forEach(p => p.classList.toggle("active", p.id === name));
  tabLinks.forEach(l => l.classList.toggle("active", l.dataset.tab === name));

  const restoreY = tabScrollPositions[name] || 0;
  requestAnimationFrame(() => window.scrollTo(0, restoreY));

  activeTabName = name;

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

// ── Ligand blacklist management ───────────────────────────
// Lets the user see the built-in metal/blacklist tables (previously only
// visible by reading core/cofactors.py) and toggle, extend, or fully
// replace them for their own analysis runs. Nothing here mutates any
// server-side state directly: the effective lists are recomputed by the
// server per-request from what we send in gatherConfig()'s `blacklist`
// field (see core.cofactors.build_effective_lists), so concurrent users
// never affect each other's blacklist. State is kept in localStorage only
// so a single browser remembers its own customization across reloads.

const BLACKLIST_STORAGE_KEY = "vhelibs_blacklist_v1";

// `activeEntries`: the current base list to show (server defaults, or a
// user-uploaded replacement) — [{code, name, category}].
// `disabledCodes`: Set of codes from `activeEntries` the user unchecked.
// `customEntries`: user-added entries layered on top, each
//   {code, name, category, enabled}.
// `replaceDicts`: null, or {metals, ligand_blacklist} when activeEntries
//   came from an uploaded file (sent to the server as blacklist.replace).
let blacklistActiveEntries = [];
let blacklistDisabledCodes = new Set();
let blacklistCustomEntries = [];
let blacklistReplaceDicts = null;
let blacklistPendingUpload = null; // {entries, metals, ligand_blacklist} awaiting Apply/Cancel

function loadBlacklistState() {
  try {
    const raw = localStorage.getItem(BLACKLIST_STORAGE_KEY);
    if (!raw) return;
    const saved = JSON.parse(raw);
    blacklistDisabledCodes = new Set(saved.disabled || []);
    blacklistCustomEntries = Array.isArray(saved.custom) ? saved.custom : [];
    blacklistReplaceDicts = saved.replace || null;
  } catch (e) {
    // Corrupt/old localStorage payload — ignore and start fresh rather than breaking the page.
    blacklistDisabledCodes = new Set();
    blacklistCustomEntries = [];
    blacklistReplaceDicts = null;
  }
}

function saveBlacklistState() {
  localStorage.setItem(BLACKLIST_STORAGE_KEY, JSON.stringify({
    disabled: Array.from(blacklistDisabledCodes),
    custom: blacklistCustomEntries,
    replace: blacklistReplaceDicts,
  }));
}

function entriesFromReplaceDicts(dicts) {
  const entries = [];
  Object.entries(dicts.ligand_blacklist || {}).forEach(([code, name]) => entries.push({ code, name, category: "blacklist" }));
  Object.entries(dicts.metals || {}).forEach(([code, name]) => entries.push({ code, name, category: "metal" }));
  return entries;
}

function fetchBlacklistDefaults() {
  // If the user previously replaced the list from a file, keep showing
  // that instead of re-fetching the server defaults over it.
  if (blacklistReplaceDicts) {
    blacklistActiveEntries = entriesFromReplaceDicts(blacklistReplaceDicts);
    renderBlacklist();
    return;
  }
  fetch("/api/blacklist")
    .then(r => r.json())
    .then(data => {
      blacklistActiveEntries = data.entries || [];
      renderBlacklist();
    })
    .catch(() => {
      const list = document.getElementById("blacklistListBlacklist");
      const metal = document.getElementById("blacklistListMetal");
      const msg = '<p class="blacklist-entry-empty">Could not load the blacklist from the server.</p>';
      if (list) list.innerHTML = msg;
      if (metal) metal.innerHTML = "";
    });
}

function makeBlacklistRow(entry, isCustom) {
  const row = document.createElement("label");
  row.className = "blacklist-entry" + (isCustom ? " is-custom" : "");

  const checkbox = document.createElement("input");
  checkbox.type = "checkbox";
  checkbox.checked = isCustom ? !!entry.enabled : !blacklistDisabledCodes.has(entry.code);
  checkbox.addEventListener("change", () => {
    if (isCustom) {
      entry.enabled = checkbox.checked;
    } else if (checkbox.checked) {
      blacklistDisabledCodes.delete(entry.code);
    } else {
      blacklistDisabledCodes.add(entry.code);
    }
    saveBlacklistState();
  });

  const code = document.createElement("span");
  code.className = "blacklist-entry-code";
  code.textContent = entry.code;

  const name = document.createElement("span");
  name.className = "blacklist-entry-name";
  name.title = entry.name;
  name.textContent = entry.name;

  row.appendChild(checkbox);
  row.appendChild(code);
  row.appendChild(name);

  if (isCustom) {
    const tag = document.createElement("span");
    tag.className = "blacklist-entry-custom-tag";
    tag.textContent = "custom";
    row.appendChild(tag);

    const removeBtn = document.createElement("button");
    removeBtn.type = "button";
    removeBtn.className = "blacklist-entry-remove";
    removeBtn.textContent = "✕";
    removeBtn.title = "Remove this custom entry";
    removeBtn.addEventListener("click", ev => {
      ev.preventDefault();
      ev.stopPropagation();
      blacklistCustomEntries = blacklistCustomEntries.filter(e => e !== entry);
      saveBlacklistState();
      renderBlacklist();
    });
    row.appendChild(removeBtn);
  }

  return row;
}

function renderBlacklist() {
  const search = (document.getElementById("blacklistSearch").value || "").trim().toLowerCase();
  const matches = e => !search || e.code.toLowerCase().includes(search) || (e.name || "").toLowerCase().includes(search);

  const columns = {
    blacklist: { container: document.getElementById("blacklistListBlacklist"), count: document.getElementById("blacklistCountBlacklist") },
    metal:     { container: document.getElementById("blacklistListMetal"),     count: document.getElementById("blacklistCountMetal") },
  };

  Object.entries(columns).forEach(([category, { container, count }]) => {
    container.innerHTML = "";
    const active = blacklistActiveEntries.filter(e => e.category === category && matches(e));
    const custom = blacklistCustomEntries.filter(e => e.category === category && matches(e));

    let enabledCount = 0;
    active.forEach(e => {
      if (!blacklistDisabledCodes.has(e.code)) enabledCount++;
      container.appendChild(makeBlacklistRow(e, false));
    });
    custom.forEach(e => {
      if (e.enabled) enabledCount++;
      container.appendChild(makeBlacklistRow(e, true));
    });

    if (!active.length && !custom.length) {
      container.innerHTML = '<p class="blacklist-entry-empty">No entries match.</p>';
    }
    const total = active.length + custom.length;
    count.textContent = total ? `(${enabledCount}/${total} active)` : "";
  });
}

document.getElementById("blacklistSearch").addEventListener("input", renderBlacklist);

document.getElementById("blacklistSelectAllBtn").addEventListener("click", () => {
  blacklistDisabledCodes.clear();
  blacklistCustomEntries.forEach(e => { e.enabled = true; });
  saveBlacklistState();
  renderBlacklist();
});

document.getElementById("blacklistSelectNoneBtn").addEventListener("click", () => {
  blacklistActiveEntries.forEach(e => blacklistDisabledCodes.add(e.code));
  blacklistCustomEntries.forEach(e => { e.enabled = false; });
  saveBlacklistState();
  renderBlacklist();
});

document.getElementById("blacklistResetBtn").addEventListener("click", () => {
  blacklistDisabledCodes = new Set();
  blacklistCustomEntries = [];
  blacklistReplaceDicts = null;
  saveBlacklistState();
  fetchBlacklistDefaults();
});

// Add a custom entry
document.getElementById("blacklistAddBtn").addEventListener("click", () => {
  const codeInput = document.getElementById("blacklistNewCode");
  const nameInput = document.getElementById("blacklistNewName");
  const categorySelect = document.getElementById("blacklistNewCategory");

  const code = codeInput.value.trim().toUpperCase();
  if (!code) { alert("Please enter a component code (e.g. the 3-letter PDB ligand ID)."); return; }

  const allCodes = new Set([
    ...blacklistActiveEntries.map(e => e.code),
    ...blacklistCustomEntries.map(e => e.code),
  ]);
  if (allCodes.has(code)) { alert(`"${code}" is already in the list.`); return; }

  blacklistCustomEntries.push({
    code,
    name: nameInput.value.trim() || code,
    category: categorySelect.value,
    enabled: true,
  });
  saveBlacklistState();
  codeInput.value = "";
  nameInput.value = "";
  renderBlacklist();
});

// Replace whole list from an uploaded file
const blacklistFileInput = document.getElementById("blacklistFileInput");
const blacklistFileName = document.getElementById("blacklistFileName");
const blacklistFilePreview = document.getElementById("blacklistFilePreview");
const blacklistFilePreviewText = document.getElementById("blacklistFilePreviewText");

blacklistFileInput.addEventListener("change", () => {
  const file = blacklistFileInput.files && blacklistFileInput.files[0];
  if (!file) return;

  const MAX_SIZE = 1024 * 1024;
  if (file.size > MAX_SIZE) {
    alert("File is too large for a blacklist file.");
    blacklistFileInput.value = "";
    return;
  }

  blacklistFileName.textContent = file.name;
  const reader = new FileReader();
  reader.onload = () => {
    const text = String(reader.result || "");
    fetch("/api/blacklist/parse", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ text }),
    })
      .then(r => r.json())
      .then(data => {
        if (data.error) { alert(data.error); return; }
        blacklistPendingUpload = data;
        const nBlacklist = Object.keys(data.ligand_blacklist || {}).length;
        const nMetal = Object.keys(data.metals || {}).length;
        blacklistFilePreviewText.textContent =
          `"${file.name}" defines ${nBlacklist} blacklist + ${nMetal} metal entries. ` +
          `Applying will replace your current list entirely (custom entries you added are kept).`;
        blacklistFilePreview.classList.remove("hidden");
      })
      .catch(() => alert("Could not parse this file."));
  };
  reader.onerror = () => alert(`Could not read file "${file.name}".`);
  reader.readAsText(file);
  blacklistFileInput.value = "";
});

document.getElementById("blacklistApplyFileBtn").addEventListener("click", () => {
  if (!blacklistPendingUpload) return;
  blacklistReplaceDicts = {
    metals: blacklistPendingUpload.metals || {},
    ligand_blacklist: blacklistPendingUpload.ligand_blacklist || {},
  };
  blacklistDisabledCodes = new Set(); // fresh start: everything in the new list is enabled
  blacklistActiveEntries = entriesFromReplaceDicts(blacklistReplaceDicts);
  blacklistPendingUpload = null;
  blacklistFilePreview.classList.add("hidden");
  blacklistFileName.textContent = "";
  saveBlacklistState();
  renderBlacklist();
});

document.getElementById("blacklistCancelFileBtn").addEventListener("click", () => {
  blacklistPendingUpload = null;
  blacklistFilePreview.classList.add("hidden");
  blacklistFileName.textContent = "";
});

loadBlacklistState();
fetchBlacklistDefaults();

// Human-readable labels for the numeric fields below, used to build a
// useful error message if one of them fails to parse (e.g. left empty or
// containing non-numeric text) instead of silently sending NaN -> null to
// the backend.
const THRESHOLD_LABELS = {
  rsr_upper:      "RSR upper threshold",
  rsr_lower:      "RSR lower threshold",
  rscc_min:       "RSCC minimum",
  rfree_max:      "R-free maximum",
  occupancy_min:  "Occupancy minimum",
  tolerance:      "Tolerance",
  distance:       "Distance",
  owab_max:       "OWAB maximum",
  resolution_max: "Resolution maximum",
  rdiff_max:      "R-diff maximum",
  dpi_max:        "DPI maximum",
};

/**
 * Collect the analysis form into a config object ready to POST.
 *
 * Throws a plain Error with a user-facing message if any numeric field
 * fails to parse, so the caller can surface it (via alert()) instead of
 * silently sending NaN -> null in the JSON payload, which would otherwise
 * make the backend fall back to its own defaults without the user knowing
 * their input was ignored.
 */
function gatherConfig() {
  const v = id => document.getElementById(id).value;
  const b = id => document.getElementById(id).checked;
  const cfg = {
    pdbids:          pdbInput.value,
    rsr_upper:       parseFloat(v("th_rsr_upper")),
    rsr_lower:       parseFloat(v("th_rsr_lower")),
    rscc_min:        parseFloat(v("th_rscc_min")),
    rfree_max:       parseFloat(v("th_rfree_max")),
    occupancy_min:   parseFloat(v("th_occupancy_min")),
    tolerance:       parseInt(v("th_tolerance"), 10),
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
    use_cache:       b("chk_use_cache"),
    blacklist: {
      disabled: Array.from(blacklistDisabledCodes),
      custom: blacklistCustomEntries
        .filter(e => e.enabled)
        .map(({ code, name, category }) => ({ code, name, category })),
      replace: blacklistReplaceDicts,
    },
  };

  // Only validate the thresholds that are actually "active" for this run:
  // e.g. rdiff_max/dpi_max are irrelevant (and may legitimately be blank)
  // when use_rdiff/use_dpi are unchecked.
  const mustValidate = Object.keys(THRESHOLD_LABELS).filter(key => {
    if (key === "owab_max")       return cfg.check_owab;
    if (key === "resolution_max") return cfg.check_resolution;
    if (key === "rdiff_max")      return cfg.use_rdiff;
    if (key === "dpi_max")        return cfg.use_dpi;
    return true;
  });
  const invalid = mustValidate.filter(key => Number.isNaN(cfg[key]));
  if (invalid.length) {
    const names = invalid.map(k => THRESHOLD_LABELS[k]).join(", ");
    throw new Error(`Please enter a valid number for: ${names}.`);
  }

  return cfg;
}

function startAnalysis() {
  let cfg;
  try {
    cfg = gatherConfig();
  } catch (err) {
    alert(err.message);
    return;
  }
  const ids = cfg.pdbids.trim();
  if (!ids) { alert("Please enter at least one PDB ID."); return; }

  // Kept so the Results-tab export (see "Results export" below) can
  // report the exact parameters this run was submitted with, even though
  // gatherConfig() itself isn't re-callable after the form may have
  // changed by the time the user clicks "Export".
  lastAnalysisConfig = cfg;

  stopPolling(); // cancel any poll loop still running for a previous job

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

// Tracks the pending poll timer so a fresh analysis (or a page-level
// cleanup) can cancel any poll loop still running for a previous job —
// without this, starting a second analysis while the first is still
// polling would leave two loops fighting over the same progress bar.
let activePollTimer = null;

function stopPolling() {
  if (activePollTimer) {
    clearTimeout(activePollTimer);
    activePollTimer = null;
  }
}

function pollJob(jobId, total) {
  stopPolling();

  const POLL_MS = 1500;
  const MAX_CONSECUTIVE_ERRORS = 5;   // ~ a few network hiccups is fine
  const MAX_TOTAL_MS = 30 * 60 * 1000; // give up after 30 min either way
  const startedAt = Date.now();
  let consecutiveErrors = 0;

  const giveUp = message => {
    stopPolling();
    alert(message);
    resetAnalyseBtn();
  };

  const poll = () => {
    if (Date.now() - startedAt > MAX_TOTAL_MS) {
      giveUp("The analysis is taking much longer than expected and has been abandoned. It may still finish server-side; check back later or re-submit.");
      return;
    }

    fetch("/api/status/" + jobId)
      .then(r => r.json())
      .then(data => {
        consecutiveErrors = 0;
        const done = data.progress || 0;
        const pct  = total > 0 ? Math.round((done / total) * 100) : 0;
        progressBar.style.width = pct + "%";
        progressLabel.textContent = `Processed ${done} / ${total} structure(s)…`;

        if (data.status === "done") {
          progressBar.style.width = "100%";
          progressLabel.textContent = "Analysis complete.";
          stopPolling();
          resetAnalyseBtn();
          lastAnalysisResults = data.results;
          renderResults(data.results);
          showTab("results");
        } else {
          activePollTimer = setTimeout(poll, POLL_MS);
        }
      })
      .catch(() => {
        consecutiveErrors++;
        if (consecutiveErrors >= MAX_CONSECUTIVE_ERRORS) {
          giveUp("Lost contact with the server while checking analysis progress. Please check your connection and try again.");
          return;
        }
        activePollTimer = setTimeout(poll, POLL_MS * 2);
      });
  };
  activePollTimer = setTimeout(poll, POLL_MS);
}

function resetAnalyseBtn() {
  analyseBtn.disabled = false;
  analyseBtn.textContent = "Analyse";
}

// ── Results rendering ─────────────────────────────────────

const resultsEmpty     = document.getElementById("resultsEmpty");
const resultsSummary   = document.getElementById("resultsSummary");
const resultsContainer = document.getElementById("resultsContainer");

// Quality filter toggle buttons (Ligand: Good/Dubious/Bad, Binding site:
// Good/Dubious/Bad, plus a separate toggle for structures that errored out
// entirely / weren't found). Markup is static (see index.html); every
// button starts "active" (full color) and toggles to a dimmed "inactive"
// state on click.
const filterToggles      = document.querySelectorAll(".filter-toggle");
const filterResetBtn     = document.getElementById("filterResetBtn");
const resultsFilterCount = document.getElementById("resultsFilterCount");

function clearResults() {
  resultsEmpty.classList.remove("hidden");
  resultsSummary.classList.add("hidden");
  resultsContainer.classList.add("hidden");
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

  // Fresh set of results → every filter toggle starts active again.
  filterToggles.forEach(btn => btn.classList.remove("inactive"));

  updateFilterCounts(results);

  resultsSummary.classList.remove("hidden");

  // Build every card in an off-DOM fragment first, then attach it in one
  // shot — avoids a reflow per card when the result set is large (many
  // PDB IDs submitted at once).
  resultsContainer.innerHTML = "";
  const frag = document.createDocumentFragment();
  results.forEach(r => {
    frag.appendChild(buildResultCard(r));
  });
  resultsContainer.appendChild(frag);
  resultsContainer.classList.remove("hidden");

  applyResultsFilter();
}

// Populates the "(N)" count next to each filter badge with how many
// ligands/binding sites/errored structures actually have that quality —
// always reflects the full result set, independent of which toggles are
// currently active.
function updateFilterCounts(results) {
  const tally = {
    ligand: { Good: 0, Dubious: 0, Bad: 0 },
    bs:     { Good: 0, Dubious: 0, Bad: 0 },
    errors: 0,
  };
  results.forEach(r => {
    if (r.error) { tally.errors++; return; }
    (r.ligands || []).forEach(l => {
      if (tally.ligand[l.ligand_quality] != null) tally.ligand[l.ligand_quality]++;
      if (tally.bs[l.binding_site_quality] != null) tally.bs[l.binding_site_quality]++;
    });
  });

  filterToggles.forEach(btn => {
    const span = btn.querySelector(".filter-count");
    if (!span) return;
    const axis = btn.dataset.axis;
    const n = axis === "errors" ? tally.errors : ((tally[axis] || {})[btn.dataset.value] || 0);
    span.textContent = `(${n})`;
  });
}

function toggleIsActive(btn) {
  return !btn.classList.contains("inactive");
}

// The set of quality values currently switched "on" for one axis (ligand
// or bs) — a ligand/binding-site entry is shown only if its own quality is
// in this set.
function activeValuesFor(axis) {
  const values = new Set();
  filterToggles.forEach(btn => {
    if (btn.dataset.axis === axis && toggleIsActive(btn)) values.add(btn.dataset.value);
  });
  return values;
}

// Shows only the ligand entries (and their parent structure card) whose
// ligand AND binding-site quality both have an active toggle, and shows/
// hides error cards ("Not found" structures) based on their own toggle.
function applyResultsFilter() {
  const ligandActive = activeValuesFor("ligand");
  const bsActive      = activeValuesFor("bs");
  const errorsBtn     = document.querySelector('.filter-toggle[data-axis="errors"]');
  const showErrors    = errorsBtn ? toggleIsActive(errorsBtn) : true;

  let visibleLigands = 0;
  let visibleCards = 0;

  resultsContainer.querySelectorAll(".result-card").forEach(card => {
    let cardHasMatch = false;
    card.querySelectorAll(".ligand-entry").forEach(entry => {
      const matches =
        ligandActive.has(entry.dataset.ligandQuality) &&
        bsActive.has(entry.dataset.bsQuality);
      entry.classList.toggle("hidden", !matches);
      if (matches) { cardHasMatch = true; visibleLigands++; }
    });
    card.classList.toggle("hidden", !cardHasMatch);
    if (cardHasMatch) visibleCards++;
  });

  let visibleErrors = 0;
  resultsContainer.querySelectorAll(".error-card").forEach(card => {
    card.classList.toggle("hidden", !showErrors);
    if (showErrors) visibleErrors++;
  });

  resultsFilterCount.textContent =
    `Showing ${visibleLigands} ligand(s) in ${visibleCards} structure(s)` +
    (visibleErrors ? `, plus ${visibleErrors} not-found structure(s).` : ".");
}

filterToggles.forEach(btn => {
  btn.addEventListener("click", () => {
    btn.classList.toggle("inactive");
    applyResultsFilter();
  });
});
filterResetBtn.addEventListener("click", () => {
  filterToggles.forEach(btn => btn.classList.remove("inactive"));
  applyResultsFilter();
});


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

  const ligands = r.ligands || [];

  const header = document.createElement("div");
  header.className = "result-card-header";
  header.innerHTML = `
    <span class="result-pdbid">${esc(r.pdbid.toUpperCase())}</span>
    ${r.title ? `<span class="result-title" style="margin-left:8px;font-weight:normal;color:var(--clr-muted,#666)">${esc(r.title)}</span>` : ""}
    ${r.uniprot ? `<span class="badge" style="margin-left:6px">UniProt: ${esc(r.uniprot)}</span>` : ""}
    <span class="result-badges">
      ${ligands.length} ligand(s)
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
    // Used by applyResultsFilter() to show/hide this entry based on the
    // Ligand-quality / Binding-site-quality filter selectors.
    entry.dataset.ligandQuality = l.ligand_quality || "";
    entry.dataset.bsQuality = l.binding_site_quality || "";
    entry.innerHTML = `
      <div class="ligand-entry-header">
        <div>
          <span class="badge badge-${qc2}">Ligand: ${l.ligand_quality}</span>
          <span class="badge badge-${bsQc}" style="margin-left:6px">BS: ${l.binding_site_quality}</span>
          ${l._manuallyEdited ? `<span class="review-edited-flag" title="Manually reviewed in the 3D Viewer">✎ edited</span>` : ""}
        </div>
        <button class="view-btn" data-pdbid="${esc(r.pdbid)}"
          data-source="${esc(l.source || 'PDB')}"
          data-ligand-index="${i}"
          data-ligand-quality="${esc(l.ligand_quality || '')}"
          data-bs-quality="${esc(l.binding_site_quality || '')}"
          data-residue-qualities='${escAttr(JSON.stringify(l.residue_qualities || {}))}'
          data-ligands='${escAttr(JSON.stringify(l.ligand_residues))}'
          data-bs='${escAttr(JSON.stringify(l.binding_site_residues))}'
          data-rte='${escAttr(JSON.stringify(l.residues_to_examine || []))}'
          data-boxes='${escAttr(JSON.stringify(l.density_boxes || {}))}'
          data-atoms='${escAttr(JSON.stringify(l.density_atoms || {}))}'>
          View 3D
        </button>
      </div>
      <div class="ligand-residues">
        <strong>Ligand:</strong>
        ${l.ligand_residues.map(s => `<span class="residue-tag">${esc(s)}</span>`).join("")}
        ${ligandNameLabel(l)}
      </div>
      ${l.binding_site_residues.length ? `
      <div class="ligand-residues" style="margin-top:6px">
        <strong>Binding site:</strong>
        ${l.binding_site_residues.slice(0, 8).map(s => `<span class="residue-tag">${esc(s)}</span>`).join("")}
        ${l.binding_site_residues.length > 8 ? `
        <span class="residue-extra-tags hidden">
          ${l.binding_site_residues.slice(8).map(s => `<span class="residue-tag">${esc(s)}</span>`).join("")}
        </span>
        <span class="residue-tag residue-toggle" role="button" tabindex="0"
          style="cursor:pointer" data-more-count="${l.binding_site_residues.length - 8}">+${l.binding_site_residues.length-8} more</span>` : ""}
      </div>` : ""}
      ${l.low_occupancy && l.low_occupancy.length ? `
      <div class="ligand-residues" style="margin-top:6px;color:var(--clr-dubious)">
        ⚠ Low occupancy: ${l.low_occupancy.map(s => `<span class="residue-tag">${esc(s)}</span>`).join("")}
      </div>` : ""}
      ${l.other_ligands && l.other_ligands.length ? `
      <div class="ligand-residues" style="margin-top:6px;color:var(--clr-muted);font-size:.75rem">
        ${l.other_ligands.length} other ligand residue(s) in this structure are not shown in the 3D view above.
      </div>` : ""}
    `;
    ligDiv.appendChild(entry);
  });

  // Rejected molecules
  const rejected = r.rejected || {};
  const rejKeys = Object.keys(rejected);
  if (rejKeys.length) {
    const rejDiv = document.createElement("details");
    rejDiv.style.marginTop = "12px";
    rejDiv.innerHTML = `<summary style="cursor:pointer;font-size:.8125rem;color:var(--clr-muted)">${rejKeys.length} rejected molecule(s)</summary>
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
      const pdbid  = btn.dataset.pdbid;
      const source = btn.dataset.source || "PDB";
      let ligands, bs, rte, boxes, atoms, residueQualities;
      try {
        ligands = JSON.parse(btn.dataset.ligands || "[]");
        bs      = JSON.parse(btn.dataset.bs || "[]");
        rte     = JSON.parse(btn.dataset.rte || "[]");
        boxes   = JSON.parse(btn.dataset.boxes || "{}");
        atoms   = JSON.parse(btn.dataset.atoms || "{}");
        residueQualities = JSON.parse(btn.dataset.residueQualities || "{}");
      } catch (err) {
        console.error("[VHELIBS] Corrupt residue/box data on View 3D button:", err);
        alert(`Could not open the 3D viewer for ${pdbid.toUpperCase()}: its residue data looks corrupted. Try re-running the analysis.`);
        return;
      }
      const reviewData = {
        ligandIndex: Number(btn.dataset.ligandIndex),
        ligandQuality: btn.dataset.ligandQuality || null,
        bsQuality: btn.dataset.bsQuality || null,
        residueQualities,
      };
      openViewer(pdbid, ligands, bs, rte, boxes, atoms, source, reviewData);
    });
  });

  // Wire binding-site "+N more" tags: clicking (or Enter/Space, since
  // these are <span role="button"> rather than real buttons) reveals the
  // rest of the binding-site residues and turns the tag into a "Show
  // less" toggle to collapse them again.
  card.querySelectorAll(".residue-toggle").forEach(toggle => {
    const expandToggle = e => {
      e.stopPropagation();
      const extra = toggle.previousElementSibling;
      if (!extra || !extra.classList.contains("residue-extra-tags")) return;
      const nowHidden = extra.classList.toggle("hidden");
      toggle.textContent = nowHidden ? `+${toggle.dataset.moreCount} more` : "Show less";
    };
    toggle.addEventListener("click", expandToggle);
    toggle.addEventListener("keydown", e => {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        expandToggle(e);
      }
    });
  });

  return card;
}

// ── Results export (.xlsx) ────────────────────────────────
// Builds a two-sheet workbook from the most recently completed analysis:
//   "Parameters" – the thresholds/options that run was submitted with.
//   "Ligands"    – one row per ligand entry, plus one (mostly-blank) row
//                  per structure that errored out / wasn't found.
// Built entirely client-side with SheetJS (loaded from CDN in
// index.html) — a real spreadsheet with two sheets, since plain .csv has
// no concept of multiple sheets. The only network calls made here are the
// lightweight electron-density-map existence checks below (never a full
// map download).

let lastAnalysisConfig  = null; // cfg object POSTed for the most recent /api/analyse call
let lastAnalysisResults = null; // data.results from the most recently completed job

const exportResultsBtn = document.getElementById("exportResultsBtn");
const exportStatus     = document.getElementById("exportStatus");

// [config key, English column label] — order here is the order rows
// appear in the "Parameters" sheet.
const PARAM_LABELS = [
  ["pdbids",           "PDB IDs / UniProt IDs (as entered)"],
  ["rsr_upper",        "RSR upper threshold (Bad above this)"],
  ["rsr_lower",        "RSR lower threshold (Good below this)"],
  ["rscc_min",         "RSCC minimum"],
  ["rfree_max",        "R-free maximum"],
  ["occupancy_min",    "Occupancy minimum"],
  ["tolerance",        "Tolerance"],
  ["distance",         "Binding site distance (Å)"],
  ["use_pdb_redo",     "Use PDB-REDO structures"],
  ["check_owab",       "Check OWAB"],
  ["owab_max",         "OWAB maximum"],
  ["check_resolution", "Check resolution"],
  ["resolution_max",   "Resolution maximum (Å)"],
  ["use_rdiff",        "Use R-diff filter"],
  ["rdiff_max",        "R-diff maximum"],
  ["use_dpi",          "Use DPI filter"],
  ["dpi_max",          "DPI maximum"],
  ["use_cache",        "Use cached downloads"],
];

function formatParamValue(v) {
  if (v === true) return "Yes";
  if (v === false) return "No";
  if (v === null || v === undefined) return "";
  return v;
}

function buildParametersSheetData(cfg) {
  const rows = [["Parameter", "Value"]];
  PARAM_LABELS.forEach(([key, label]) => rows.push([label, formatParamValue(cfg[key])]));

  const bl = cfg.blacklist || {};
  rows.push(["Blacklist: disabled default entries", (bl.disabled || []).join(", ")]);
  rows.push(["Blacklist: custom entries added", (bl.custom || []).map(e => `${e.code} (${e.name})`).join(", ")]);
  rows.push(["Blacklist: default list replaced by uploaded file", bl.replace ? "Yes" : "No"]);
  rows.push(["Export date", new Date().toISOString()]);
  return rows;
}

function fmtNum(v) {
  return (v != null && !Number.isNaN(v)) ? v : null;
}

function rejectedResiduesText(rejected) {
  const keys = Object.keys(rejected || {});
  if (!keys.length) return "";
  return keys.map(k => `${k}`).join("; ");
}

// Placeholder used throughout the xlsx export whenever a piece of data
// genuinely doesn't exist for a row.
const XLSX_NOT_AVAILABLE = "Not Available";

function orNA(value) {
  return (value === null || value === undefined || value === "") ? XLSX_NOT_AVAILABLE : value;
}

// One name per entry in l.ligand_residues (same order).
function ligandNameLabel(l) {
  const names = [...new Set((l.ligand_names || []).filter(Boolean))];
  if (!names.length) return "";
  return `<span class="ligand-name-label" style="margin-left:6px;color:var(--clr-muted,#666);white-space:normal;overflow-wrap:break-word;display:inline-block;max-width:100%">${esc(names.join(" / "))}</span>`;
}

function formatLigandNames(l) {
  const names = l.ligand_names || [];
  if (!names.length) return "";
  return names.map(n => n || XLSX_NOT_AVAILABLE).join("; ");
}

// Binding-site residues annotated with how each was classified. A residue
// absent from residue_qualities was never flagged as dubious/bad to begin
// with, so it's Good by construction.
const QUALITY_LETTER = { Bad: "B", Dubious: "D", Good: "G" };
const QUALITY_SEVERITY = { Bad: 0, Dubious: 1, Good: 2 };

function formatBindingSiteResidues(l) {
  const residues = l.binding_site_residues || [];
  if (!residues.length) return "";
  const qualities = l.residue_qualities || {};
  return residues
    .map(res => ({ res, quality: qualities[res] || "Good" }))
    .sort((a, b) => QUALITY_SEVERITY[a.quality] - QUALITY_SEVERITY[b.quality])
    .map(({ res, quality }) => `${res} (${QUALITY_LETTER[quality]})`)
    .join("; ");
}

/**
 * Resolve electron-density-map availability for every successfully
 * analysed structure in `results`, WITHOUT downloading any map.
 *
 * - PDB-REDO analyses: no request needed at all. A structure that
 *   analysed successfully already required fetching its PDB-REDO stats
 *   (see rsr_core._fetch_structure_data), which only succeeds if
 *   PDB-REDO has data for that entry — and therefore a map obtainable
 *   from its map-maker service — so this is inferred as "available".
 * - Standard PDB analyses: one lightweight HEAD-based existence check
 *   per structure via /api/edm-exists/<pdbid> (see routes.py), a few
 *   requests at a time.
 *
 * Returns a Map<lowercase pdbid, true|false|null>, null meaning
 * "could not be determined".
 */
async function fetchDensityMapAvailability(results) {
  const availability = new Map();
  const ok = results.filter(r => !r.error);
  if (!ok.length) return availability;

  const firstSource = (ok.find(r => r.ligands && r.ligands.length) || {}).ligands?.[0]?.source;
  const isPdbRedo = firstSource === "PDB_REDO";

  if (isPdbRedo) {
    ok.forEach(r => availability.set(r.pdbid.toLowerCase(), true));
    return availability;
  }

  const pdbids = ok.map(r => r.pdbid.toLowerCase());
  const CONCURRENCY = 6;
  let i = 0;
  async function worker() {
    while (i < pdbids.length) {
      const pdbid = pdbids[i++];
      try {
        const resp = await fetch(`/api/edm-exists/${pdbid}`);
        const data = await resp.json();
        availability.set(pdbid, !!data.exists);
      } catch (e) {
        availability.set(pdbid, null); // couldn't determine — left blank in the sheet
      }
    }
  }
  await Promise.all(Array.from({ length: Math.min(CONCURRENCY, pdbids.length) }, worker));
  return availability;
}

function buildLigandsSheetData(results, densityAvailability) {
  const rows = [[
    "UniProt ID", "PDB ID", "Ligand ID", "Ligand Name",
    "Ligand Class.", "BS Class.", "BS Residues",
    "R-free", "R-work", "Rejected Molecules",
    "EDM Available",
  ]];

  results.forEach(r => {
    const complex = (r.pdbid || "").toUpperCase();
    const uniprot = orNA(r.uniprot);

    if (r.error) {
      // The backend never got past the initial fetch/parse for this entry
      // (not found, download failure, no ligand atoms at all, …), so
      // nothing beyond the UniProt tag and complex ID was ever computed.
      rows.push([uniprot, complex,
        XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE,
        XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE]);
      return;
    }

    const avail = densityAvailability.get(r.pdbid.toLowerCase());
    const mapAvail = avail === true ? "Yes" : avail === false ? "No" : XLSX_NOT_AVAILABLE;
    const sd = r.struc_dict || {};
    const rFree = orNA(fmtNum(sd.rFree));
    const rWork = orNA(fmtNum(sd.rWork));
    const rejectedText = rejectedResiduesText(r.rejected);
    const ligands = r.ligands || [];

    if (!ligands.length) {
      // Analysed successfully, but nothing here qualified as a ligand
      // (e.g. everything present was blacklisted): fill in every
      // structure-level field we do have, and mark the per-ligand columns
      // — which genuinely don't apply to this row — as not available.
      rows.push([uniprot, complex,
        XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE, XLSX_NOT_AVAILABLE,
        rFree, rWork, rejectedText, mapAvail]);
      return;
    }

    ligands.forEach(l => {
      rows.push([
        uniprot,
        complex,
        orNA((l.ligand_residues || []).join("; ")),
        orNA(formatLigandNames(l)),
        orNA(l.ligand_quality),
        orNA(l.binding_site_quality),
        orNA(formatBindingSiteResidues(l)),
        rFree,
        rWork,
        rejectedText,
        mapAvail,
      ]);
    });
  });

  return rows;
}

const XLSX_MAX_COL_WIDTH = 30;

// Longest cell value per column, capped at XLSX_MAX_COL_WIDTH
function computeColumnWidths(rows) {
  const colCount = rows.reduce((max, row) => Math.max(max, row.length), 0);
  const widths = new Array(colCount).fill(8);
  rows.forEach(row => {
    row.forEach((cell, i) => {
      const len = cell == null ? 0 : String(cell).length;
      widths[i] = Math.min(XLSX_MAX_COL_WIDTH, Math.max(widths[i], len));
    });
  });
  return widths;
}

// Adds `rows` (array-of-arrays, first row = header) as a new worksheet on
// `workbook`, with column widths capped at XLSX_MAX_COL_WIDTH and real
// "Wrap Text" cell formatting applied throughout, so long content
// continues on further lines inside the same cell instead of widening the
// column indefinitely.
function addWrappedSheet(workbook, name, rows, noWrapColumns = []) {
  const ws = workbook.addWorksheet(name);
  ws.addRows(rows);
  computeColumnWidths(rows).forEach((w, i) => { ws.getColumn(i + 1).width = w; });

  const header = rows[0] || [];
  const noWrapIndices = new Set(
    noWrapColumns.map(colName => header.indexOf(colName)).filter(i => i >= 0)
  );

  ws.eachRow(row => {
    row.eachCell({ includeEmpty: true }, (cell, colNumber) => {
      cell.alignment = { wrapText: !noWrapIndices.has(colNumber - 1), vertical: "top" };
    });
  });
  return ws;
}

async function exportResults() {
  if (!lastAnalysisResults || !lastAnalysisResults.length) {
    alert("There are no results to export yet. Run an analysis first.");
    return;
  }
  if (typeof ExcelJS === "undefined") {
    alert("The export library failed to load. Check your internet connection and try again.");
    return;
  }

  exportResultsBtn.disabled = true;
  exportStatus.classList.remove("hidden");
  exportStatus.textContent = "Checking electron density map availability…";

  try {
    const densityAvailability = await fetchDensityMapAvailability(lastAnalysisResults);

    const wb = new ExcelJS.Workbook();

    // Ligands sheet first, since it's the one people actually want to
    // look at; Parameters is reference material for the run that produced
    // it, kept second.
    addWrappedSheet(wb, "Ligands", buildLigandsSheetData(lastAnalysisResults, densityAvailability),
      ["BS Residues"]);
    addWrappedSheet(wb, "Parameters", buildParametersSheetData(lastAnalysisConfig || {}));

    const buffer = await wb.xlsx.writeBuffer();
    const blob = new Blob([buffer], {
      type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    });
    const stamp = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `vhelibs_results_${stamp}.xlsx`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  } catch (err) {
    console.error("[VHELIBS] Export failed:", err);
    alert("Could not export the results: " + err.message);
  } finally {
    exportResultsBtn.disabled = false;
    exportStatus.classList.add("hidden");
  }
}

exportResultsBtn.addEventListener("click", exportResults);

// ── 3D Viewer (Mol*) ──────────────────────────────────────

let viewerInitPromise = null;   // Promise<Viewer>, created once
let viewerInstance    = null;   // resolved Mol* Viewer (has .plugin)
let isLoadingStructure = false;
let currentPdbId      = null;
let currentLigandRes  = [];
let currentBsRes      = [];
let currentRteRes     = [];
let currentDensityBoxes = null; // {ligand, binding_site, residues_to_examine} bboxes, or null if unknown
let currentDensityAtoms = null; // {ligand, binding_site, residues_to_examine} per-atom centers, or null
let currentSource = "PDB"; // "PDB" or "PDB_REDO" — decides where density comes from, see buildMvsData
let currentAtomRadius = 1.6; // Å, atom-mask radius sent to /api/density-mask (user-adjustable)
let currentFocusRes = null; // residue string clicked in the "residues to examine" list, or null
const layerState = { protein: true, ligand: true, bs: true }; // structure checkbox state
const densityLayerState = { ligand: true, bs: true, rte: true }; // density checkbox state
let currentIsovalue = 1.0; // relative sigma units

const molContainer        = document.getElementById("nglContainer"); // id kept for CSS compat
const viewerPdbInput      = document.getElementById("viewerPdbInput");
const viewerDensityControls = document.getElementById("viewerDensityControls");

// Quality review panel (right-hand column) — see "Quality review panel" section below.
const viewerReviewPanel      = document.getElementById("viewerReviewPanel");
const reviewPanelPlaceholder = document.getElementById("reviewPanelPlaceholder");
const reviewComponentsSection = document.getElementById("reviewComponentsSection");
const reviewOverallSection   = document.getElementById("reviewOverallSection");
const reviewActions          = document.getElementById("reviewActions");
const reviewComponentsList   = document.getElementById("reviewComponentsList");
const reviewComponentsEmpty  = document.getElementById("reviewComponentsEmpty");
const reviewLigandQualityBtns = document.getElementById("reviewLigandQualityBtns");
const reviewBsQualityBtns     = document.getElementById("reviewBsQualityBtns");
const reviewResetBtn         = document.getElementById("reviewResetBtn");
const reviewConfirmBtn       = document.getElementById("reviewConfirmBtn");
const reviewStatus           = document.getElementById("reviewStatus");

// "Reload Viewer" button — recovers from a hung/frozen Mol* instance
// (typically a stuck WebGL context after heavy density streaming) without
// requiring a full page reload, which would otherwise wipe out the
// Results tab and force re-running the whole analysis job. Declared in
// index.html; here we just grab the reference and wire it up below.
const reloadViewerBtn = document.getElementById("reloadViewerBtn");
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

function openViewer(pdbid, ligandRes, bsRes, rteRes, densityBoxes, densityAtoms, source, reviewData) {
  showTab("viewer");
  viewerPdbInput.value = pdbid.toUpperCase();
  loadMolstarStructure(pdbid, ligandRes, bsRes, rteRes || [], densityBoxes || null, densityAtoms || null, source || "PDB", reviewData || null);
}

async function loadMolstarStructure(pdbid, ligandRes, bsRes, rteRes, densityBoxes, densityAtoms, source, reviewData) {
  if (isLoadingStructure) {
    console.warn("loadMolstarStructure called while a load is already in progress — ignoring.");
    return;
  }
  isLoadingStructure = true;
  reloadViewerBtn.disabled = true;
  reloadViewerBtn.textContent = "Loading…";
  showReviewPlaceholder();
  viewerDensityControls.classList.add("hidden");

  currentPdbId = pdbid;
  currentLigandRes = ligandRes;
  currentBsRes = bsRes;
  currentRteRes = rteRes || [];
  currentDensityBoxes = densityBoxes || null;
  currentDensityAtoms = densityAtoms || null;
  currentSource = source || "PDB";
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
    setupReviewContext(pdbid, reviewData);
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
  reloadViewerBtn.disabled = false;
  reloadViewerBtn.textContent = "Reload Viewer";
}

// ── Reload Viewer ────────────────────────────────────────────────────────
//
// Mol* occasionally gets stuck — most often after heavy density-streaming
// work (segmented 2Fo-Fc boxes clipped per-atom) — leaving the canvas
// frozen/unresponsive. Previously the only fix was a full page reload,
// which wipes the Results tab and forces re-running the whole analysis
// job. This tears down and recreates just the Mol* plugin instance,
// reusing the region/PDB state already sitting in memory (currentPdbId,
// currentLigandRes, etc.) — nothing round-trips to the server, and
// nothing else in the page is touched.
//
// Note: this restores *which* structure/ligand/binding-site/RTE was
// loaded, but resets the viewer's own cosmetic state (layer/density
// checkboxes, isovalue, atom-mask radius, focused residue) to the same
// defaults a fresh "View 3D" click would use — same as loadMolstarStructure
// always does. If you'd rather have those survive a reload too, they're
// all cheap to snapshot/restore here; just say the word.
async function reloadViewer() {
  if (isLoadingStructure) return;
  if (!currentPdbId) {
    alert("No structure is currently loaded — nothing to reload.");
    return;
  }

  // Snapshot the region/PDB state to restore after recreating the plugin.
  const pdbid        = currentPdbId;
  const ligandRes     = currentLigandRes;
  const bsRes         = currentBsRes;
  const rteRes        = currentRteRes;
  const densityBoxes  = currentDensityBoxes;
  const densityAtoms  = currentDensityAtoms;
  const source        = currentSource;
  // Also snapshot any not-yet-confirmed review-panel edits — a reload is a
  // recovery action for a stuck viewer, not a reason to throw away
  // in-progress classification overrides the user hasn't confirmed yet.
  const reviewSnapshot = reviewContext ? JSON.parse(JSON.stringify(reviewContext)) : null;

  reloadViewerBtn.disabled = true;
  reloadViewerBtn.textContent = "Reloading…";

  // Tear down the (possibly hung) plugin. dispose() forces WebGL context
  // loss by default, which is exactly what's needed to actually free a
  // stuck GPU/WebGL context — re-rendering into the same context can't do
  // that, only a full dispose+recreate can. Wrapped in try/catch because a
  // frozen instance is precisely the case where dispose() itself might
  // throw or misbehave; either way we still want to proceed and rebuild.
  try {
    if (viewerInstance) {
      viewerInstance.dispose();
    }
  } catch (err) {
    console.warn("[VHELIBS] viewer.dispose() threw during reload (continuing anyway):", err);
  } finally {
    viewerInstance = null;
    viewerInitPromise = null;
    molContainer.innerHTML = "";
  }

  try {
    await loadMolstarStructure(pdbid, ligandRes, bsRes, rteRes, densityBoxes, densityAtoms, source, null);
    if (reviewSnapshot) {
      reviewContext = reviewSnapshot;
      reviewPanelPlaceholder.classList.add("hidden");
      reviewComponentsSection.classList.remove("hidden");
      reviewOverallSection.classList.remove("hidden");
      reviewActions.classList.remove("hidden");
      renderReviewPanel();
    }
  } finally {
    reloadViewerBtn.disabled = false;
    reloadViewerBtn.textContent = "Reload Viewer";
  }
}

reloadViewerBtn.addEventListener("click", () => {
  reloadViewer().catch(err => console.error("[VHELIBS] Reload viewer failed:", err));
});

// RCSB download endpoints to try, in order of preference. For PDB-REDO
// results we must NOT use the original RCSB coordinates: PDB-REDO
// re-refines/rebuilds the model (peptide flips, rotamer changes, moved
// waters, etc.), and its density map (see pdbRedoMapUrl) is calculated
// against ITS OWN "final" coordinates, not the original deposited ones.
// The per-atom centers used to mask the density (density_atoms) are
// already computed from that same PDB-REDO model (see
// core.pdb_utils.get_pdb_file(pdbid, pdb_redo=True) on the backend), so the
// displayed structure must match it too, or the density can end up
// masked around positions that don't correspond to the shown atoms.
function molstarSources(pdbidLower) {
  if (currentSource === "PDB_REDO") {
    return [
      { url: `https://pdb-redo.eu/db/${pdbidLower}/${pdbidLower}_final.cif`, format: "mmcif" },
    ];
  }
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

  // For PDB-REDO analyses the electron density comes from an external
  // map-maker service that doesn't always have (or correctly return) a
  // map for a given entry — see core.pdb_redo_utils.get_EDM, which now
  // validates the CCP4 magic bytes server-side and 404s if it's not a
  // real map. Probe that once per scene build with a cheap HEAD request
  // so a missing/broken map just means "no density shown" instead of
  // Mol* throwing mid-load and aborting the whole structure.
  let redoMapAvailable = false;
  if (currentSource === "PDB_REDO" && currentDensityBoxes && hasAnyBox(currentDensityBoxes)) {
    redoMapAvailable = await checkPdbRedoMapAvailable(currentPdbId);
    if (!redoMapAvailable) {
      console.warn(`[VHELIBS] No usable PDB-REDO density map for ${currentPdbId} — showing structure without density.`);
    }
  }

  for (const src of sources) {
    try {
      const mvsData = buildMvsData({
        sourceUrl: src.url,
        sourceFormat: src.format,
        ligandRes: currentLigandRes,
        bsRes: currentBsRes,
        layers: layerState,
        focus,
        densityBoxes: currentDensityBoxes,
        densityAtoms: currentDensityAtoms,
        atomRadius: currentAtomRadius,
        densityLayers: densityLayerState,
        isovalue: currentIsovalue,
        focusResidue: currentFocusRes,
        mapSource: currentSource,
        redoMapAvailable,
      });
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

// HEAD-checks the PDB-REDO map endpoint without downloading the (possibly
// large) body, so we know up front whether it's safe to include a volume
// node for it in the MVS tree. Bounded with an explicit timeout: without
// one, a hung backend/proxy would leave applyMolstarScene() waiting
// indefinitely before it can even start building the scene, which reads
// to the user as a frozen viewer with no feedback.
async function checkPdbRedoMapAvailable(pdbid, timeoutMs = 8000) {
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs);
  try {
    const r = await fetch(pdbRedoMapUrl(pdbid), { method: "HEAD", signal: controller.signal });
    return r.ok;
  } catch {
    return false;
  } finally {
    clearTimeout(timer);
  }
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

function buildMvsData({
  sourceUrl, sourceFormat, ligandRes, bsRes, layers, focus,
  densityBoxes, densityAtoms, atomRadius, densityLayers, isovalue,
  focusResidue, mapSource, redoMapAvailable,
}) {
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
  // Each region (ligand / binding site / residues-to-examine) is rendered
  // from a single, already-masked density map fetched from our own
  // backend. The backend crops the source map to the region's padded box and
  // zeroes every voxel further than `atomRadius` Angstrom from every atom of
  // that region, using gemmi's native set_points_around() mask primitive
  //
  // PDB-REDO ANALYSES ARE DIFFERENT: PDB-REDO re-refines the model, so its
  // density map does not match the original entry's map served by EBI —
  // /api/density-box would silently show the wrong (original-PDB) density
  // for a PDB-REDO result. The masking endpoint takes a `source` param for
  // exactly this reason: `source=pdb_redo` masks PDB-REDO's own full map
  // (see core.pdb_redo_utils.get_EDM) instead of the standard EBI one, so
  // no separate "shared volume" plumbing is needed here any more — each
  // region just asks for its own masked map with the right source.
  const isRedo = mapSource === "PDB_REDO";

  if (densityBoxes) {
    for (const region of DENSITY_REGIONS) {
      if (!densityLayers[region.layerKey]) continue;
      if (isRedo && !redoMapAvailable) continue;

      const box = densityBoxes[region.boxKey];
      if (!box || !box.min || !box.max) continue;

      const atoms = (densityAtoms && densityAtoms[region.boxKey]) || [];
      if (!atoms.length) continue;

      const maskUrl = densityMaskUrl(currentPdbId, region.boxKey, atoms, box, atomRadius,
        isRedo ? "pdb_redo" : "pdb");
      const volume = builder.download({ url: maskUrl }).parse({ format: "map" }).volume({});

      volume.representation({ type: "isosurface", relative_isovalue: isovalue, show_wireframe: false, show_faces: true })
        .color({ color: region.color })
        .opacity({ opacity: 0.35 });
    }
  }

  return builder.getState();
}

// Builds the URL for a region's pre-masked density map. The backend crops the
// source map to `box` and zeroes every voxel further than `radius` Å
// from every atom in `atoms`, so what comes back is ready to render as a
// single isosurface.
function densityMaskUrl(pdbid, region, atoms, box, radius, source = "pdb") {
  const min = box.min.join(",");
  const max = box.max.join(",");
  const atomsParam = atoms.map(a => a.center.join(",")).join(";");
  const params = new URLSearchParams({ min, max, atoms: atomsParam, radius, source });
  return `/api/density-mask/${pdbid.toLowerCase()}/${region}?${params.toString()}`;
}

// Builds the URL for the full PDB-REDO density map (see routes.py /api/edm
// and core.pdb_redo_utils.get_EDM). Used only by checkPdbRedoMapAvailable()
// below as a cheap preflight — the actual density shown in the viewer now
// always goes through densityMaskUrl() above, for both sources.
function pdbRedoMapUrl(pdbid) {
  return `/api/edm/${pdbid.toLowerCase()}?source=pdb_redo`;
}

// ── Quality review panel (3D Viewer, right-hand column) ────
// Lets the user override the computed Good/Dubious/Bad call for each
// "component to examine" (residues_to_examine) and for the ligand/binding
// site overall. Nothing here touches the underlying analysis — it's a
// manual annotation layer that only reaches the Results tab when the user
// clicks "Confirm changes" (see reviewConfirmBtn below).

let reviewContext = null;
// Shape while active:
// {
//   pdbid: lowercase pdbid,
//   ligandIndex: index into results[i].ligands for this structure,
//   original: { residueQualities: {resKey: "Good"|"Dubious"|"Bad"|null}, ligandQuality, bsQuality },
//   working:  same shape — the in-progress edits, discarded by Reset, written back by Confirm
// }

function computeInitialReviewState(reviewData) {
  const residueQualities = {};
  currentRteRes.forEach(res => {
    residueQualities[res] = (reviewData && reviewData.residueQualities && reviewData.residueQualities[res]) || null;
  });
  return {
    residueQualities,
    ligandQuality: (reviewData && reviewData.ligandQuality) || null,
    bsQuality: (reviewData && reviewData.bsQuality) || null,
  };
}

function showReviewPlaceholder() {
  reviewContext = null;
  reviewPanelPlaceholder.classList.remove("hidden");
  reviewComponentsSection.classList.add("hidden");
  reviewOverallSection.classList.add("hidden");
  reviewActions.classList.add("hidden");
  reviewStatus.textContent = "";
}

function setupReviewContext(pdbid, reviewData) {
  if (!reviewData || Number.isNaN(reviewData.ligandIndex)) {
    showReviewPlaceholder();
    return;
  }
  const state = computeInitialReviewState(reviewData);
  reviewContext = {
    pdbid: pdbid.toLowerCase(),
    ligandIndex: reviewData.ligandIndex,
    original: JSON.parse(JSON.stringify(state)),
    working: JSON.parse(JSON.stringify(state)),
  };
  reviewPanelPlaceholder.classList.add("hidden");
  reviewComponentsSection.classList.remove("hidden");
  reviewOverallSection.classList.remove("hidden");
  reviewActions.classList.remove("hidden");
  reviewStatus.textContent = "";
  renderReviewPanel();
}

// Builds a fresh Good/Dubious/Bad button group. `compact` uses single-letter
// labels so rows stay one line wide in the narrow 260px review column.
function buildQualityBtnGroup(currentValue, onSelect, compact) {
  const group = document.createElement("div");
  group.className = "quality-btn-group" + (compact ? " compact" : "");
  ["Good", "Dubious", "Bad"].forEach(q => {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = `quality-btn quality-btn-${q.toLowerCase()}`;
    btn.textContent = compact ? q[0] : q;
    btn.title = q;
    btn.classList.toggle("active", currentValue === q);
    btn.addEventListener("click", () => onSelect(q));
    group.appendChild(btn);
  });
  return group;
}

// Updates the active/inactive styling on a *static* (already-built)
// quality-btn-group — used for the Ligand/Binding site overall selectors,
// which live permanently in index.html rather than being rebuilt each render.
function syncQualityBtnGroup(container, currentValue) {
  container.querySelectorAll(".quality-btn").forEach(btn => {
    btn.classList.toggle("active", btn.dataset.quality === currentValue);
  });
}

function renderReviewPanel() {
  if (!reviewContext) return;
  const { working } = reviewContext;

  reviewComponentsList.innerHTML = "";
  if (!currentRteRes.length) {
    reviewComponentsEmpty.classList.remove("hidden");
  } else {
    reviewComponentsEmpty.classList.add("hidden");
    currentRteRes.forEach(res => {
      const row = document.createElement("div");
      row.className = "review-component-row";

      const label = document.createElement("button");
      label.type = "button";
      label.className = "review-component-label";
      label.textContent = res;
      label.title = "Click to focus this component in the 3D view";
      label.classList.toggle("active", currentFocusRes === res);
      label.addEventListener("click", () => {
        if (!residueSelector(res) || !currentPdbId || isLoadingStructure || !viewerInstance) return;
        currentFocusRes = res;
        reviewComponentsList.querySelectorAll(".review-component-label").forEach(b => {
          b.classList.toggle("active", b.textContent === res);
        });
        applyMolstarScene({ focus: false, keepCamera: false })
          .catch(err => console.error("[VHELIBS] Failed to focus residue:", err));
      });

      // Clicking Good/Dubious/Bad only ever edits the in-progress `working`
      // state — the component stays listed here regardless of which button
      // is active, per spec ("if marked GOOD it should show as good, but
      // without disappearing from the components to examine list").
      const group = buildQualityBtnGroup(working.residueQualities[res], q => {
        working.residueQualities[res] = q;
        renderReviewPanel();
      }, true);

      row.appendChild(label);
      row.appendChild(group);
      reviewComponentsList.appendChild(row);
    });
  }

  syncQualityBtnGroup(reviewLigandQualityBtns, working.ligandQuality);
  syncQualityBtnGroup(reviewBsQualityBtns, working.bsQuality);
}

reviewLigandQualityBtns.querySelectorAll(".quality-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    if (!reviewContext) return;
    reviewContext.working.ligandQuality = btn.dataset.quality;
    renderReviewPanel();
  });
});
reviewBsQualityBtns.querySelectorAll(".quality-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    if (!reviewContext) return;
    reviewContext.working.bsQuality = btn.dataset.quality;
    renderReviewPanel();
  });
});

reviewResetBtn.addEventListener("click", () => {
  if (!reviewContext) return;
  reviewContext.working = JSON.parse(JSON.stringify(reviewContext.original));
  renderReviewPanel();
  reviewStatus.textContent = "Reverted to the last confirmed classification.";
  setTimeout(() => { if (reviewStatus.textContent.startsWith("Reverted")) reviewStatus.textContent = ""; }, 3000);
});

// Re-renders the Results tab after a Confirm, without discarding whatever
// quality filters the user had already toggled off — renderResults() on
// its own resets every filter toggle back to active, which would be a
// jarring side effect of simply reviewing one ligand.
function refreshResultsAfterEdit(results) {
  const inactiveToggles = Array.from(filterToggles).filter(btn => !toggleIsActive(btn));
  renderResults(results);
  inactiveToggles.forEach(btn => btn.classList.add("inactive"));
  applyResultsFilter();
}

reviewConfirmBtn.addEventListener("click", () => {
  if (!reviewContext) return;
  if (!lastAnalysisResults) {
    alert("There are no results loaded to update.");
    return;
  }
  const r = lastAnalysisResults.find(res => !res.error && (res.pdbid || "").toLowerCase() === reviewContext.pdbid);
  const l = r && r.ligands && r.ligands[reviewContext.ligandIndex];
  if (!l) {
    alert("Could not find this ligand entry in the current results — they may have been re-run since this view was opened.");
    return;
  }

  l.ligand_quality = reviewContext.working.ligandQuality || l.ligand_quality;
  l.binding_site_quality = reviewContext.working.bsQuality || l.binding_site_quality;
  l.residue_qualities = Object.assign({}, l.residue_qualities, reviewContext.working.residueQualities);
  l._manuallyEdited = true;

  // The just-confirmed state becomes the new baseline: Reset from here on
  // reverts to *this*, not to the original computed classification.
  reviewContext.original = JSON.parse(JSON.stringify(reviewContext.working));

  refreshResultsAfterEdit(lastAnalysisResults);
  reviewStatus.textContent = "Changes applied to Results.";
  setTimeout(() => { if (reviewStatus.textContent === "Changes applied to Results.") reviewStatus.textContent = ""; }, 3000);
});

// Wires a group of checkboxes to a state object + scene rebuild.
// Because Mol* diffs the state tree, rebuilding the scene on every
// toggle only touches the representations that actually changed.
function wireLayerCheckboxes(ids, keys, state, failureLabel) {
  ids.forEach((id, i) => {
    const key = keys[i];
    document.getElementById(id).addEventListener("change", function () {
      state[key] = this.checked;
      if (!currentPdbId || isLoadingStructure || !viewerInstance) return;
      applyMolstarScene({ focus: false, keepCamera: true })
        .catch(err => console.error(`[VHELIBS] Failed to update ${failureLabel}:`, err));
    });
  });
}

// Structure-layer checkboxes (protein cartoon / ligand / binding site).
wireLayerCheckboxes(
  ["chkProtein", "chkLigand", "chkBS"],
  ["protein", "ligand", "bs"],
  layerState,
  "layer visibility"
);

// Density-region checkboxes (ligand / binding site / residues to examine).
wireLayerCheckboxes(
  ["chkDensityLigand", "chkDensityBS", "chkDensityRTE"],
  ["ligand", "bs", "rte"],
  densityLayerState,
  "density visibility"
);

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

// Atom mask radius slider — sent to, where it's quantized to a 0.25 Å step and
// used server-side to mask the density around each atom. Larger radius =
// more of the surrounding density is kept per atom; smaller = tighter to
// the atom itself. Debounced the same way as the isovalue slider.
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
  // Defensive fallback only: the backend (core.pdb_atom.format_reskey) now
  // always inserts an explicit space between chain and residue number, so
  // this should no longer be reachable in practice. Kept in case a residue
  // key ever arrives without that separator (e.g. malformed/legacy data),
  // so a parse failure here degrades to "residue not selectable" rather
  // than silently falling back to "show every ligand in the structure" in
  // buildMvsData.
  const m = res.match(/^(\S*?)\s*(\D*)(-?\d+)(\D*)$/);
  if (m) {
    const chain  = m[2];
    const resnum = parseInt(m[3], 10);
    if (chain && !Number.isNaN(resnum)) {
      return { auth_asym_id: chain, auth_seq_id: resnum };
    }
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

// Like esc(), but also escapes single quotes.
function escAttr(s) {
  return esc(s).replace(/'/g, "&#39;");
}
