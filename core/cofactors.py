# -*- coding: utf-8 -*-
#
#   Copyright 2012-2024 Adrià Cereto Massagué
#   Migrated to web version
#
import csv

metals = {
    '3CO': 'Cobalt', 'NAW': 'SODIUM ION', 'MO3': 'Magnesium Ion, 3 Waters coordinated',
    'MO2': 'Magnesium Ion, 2 Waters coordinated', 'MO1': 'Magnesium Ion, 1 Water coordinated',
    'MO6': 'Magnesium Ion, 6 Waters coordinated', 'MO5': 'Magnesium Ion, 5 Waters coordinated',
    'MO4': 'Magnesium Ion, 4 Waters coordinated', 'PT': 'Platinum ion', 'PB': 'LEAD (II) ION',
    'ZO3': 'Zinc', 'MH2': 'Manganese', 'MH3': 'Manganese', 'ZN': 'Zinc', 'K': 'POTASSIUM ION',
    'MN6': 'Manganese', 'MG': 'MAGNESIUM ION', 'MN': 'Manganese', 'O4M': 'Manganese',
    'FE': 'Fe(2+)', 'CD3': 'CADMIUM ION', 'CD5': 'CADMIUM ION', 'AU3': 'GOLD',
    'CU': 'COPPER ION', 'CO': 'COBALT ION', 'CL': 'CHLORIDE ION', 'CA': 'Calcium',
    'CE': 'CERIUM (III) ION', 'CD': 'CADMIUM ION', 'SE': 'Selenium', 'YB': 'Yterbium ion',
    'CS': 'CESIUM ION', 'LI': 'LITHIUM ION', '3NI': 'NICKEL (III) ION', 'NH4': 'Ammonium Ion',
    'IOD': 'IODIDE ION', '6MO': 'Molybdenum', 'RU': 'Rutenium ion', 'BR': 'BROMIDE ION',
    'KO4': 'POTASSIUM ION, 4 WATERS COORDINATED', '4MO': 'Molybdenum', 'BA': 'Barium ion',
    'ZNO': 'Zinc', 'NAO': 'SODIUM ION', 'V': 'Vanadium', 'HG': 'Mercury (II) Bound to Cys 206',
    'NA2': 'SODIUM ION', 'ZN2': 'Zinc', 'NA6': 'SODIUM ION', 'NA5': 'SODIUM ION',
    'CD1': 'CADMIUM ION', 'CR': 'CHROMIUM ION', 'ZN3': 'Zinc', 'CUA': 'DINUCLEAR COPPER ION',
    '1CU': 'COPPER ION, 1 WATER COORDINATED', 'MN3': 'Manganese', 'MN5': 'Manganese',
    'MW3': 'Manganese', 'MW2': 'Manganese', 'MW1': 'Manganese',
    '2OF': 'FERROUS ION, 2 WATERS COORDINATED', 'NI': 'Nickel ion', 'NA': 'SODIUM ION',
    'FE2': 'FE (II) ION', 'AG': 'SILVER ION', 'AL': 'ALUMINUM ION', 'AU': 'GOLD',
    'W': 'Tungsten', 'EU': 'Europium ion', 'OS': 'Osmium ion', 'HO': 'Holmium ion',
    'IR': 'Iridium(4+)', 'IR3': 'IRIDIUM (III) ION', 'U1': 'Uranium', 'SB': 'ANTIMONY (III) ION',
    'PR': 'PRASEODYMIUM ION', 'GD': 'GADOLINIUM ATOM', 'SM': 'SAMARIUM (III) ION',
    'LU': 'LUTETIUM (III) ION', 'TB': 'TERBIUM(III) ION', 'RB': 'RUBIDIUM ION',
    'Y1': 'yttrium(+2) cation', 'TL': 'THALLIUM (I) ION', 'PD': 'PALLADIUM ION',
    'SO3': 'SULFITE ION', 'SR': 'Strontium ion', 'C1O': 'CU- O LINKAGE',
    'O': 'Oxygen atom', 'F': 'FLUORIDE',
}

ligand_blacklist = {
    '': '', 'EOH': 'Ethanol', 'HEM': 'Protoprphyrin IX containing FE', 'SIN': 'Succinic acid',
    'SPV': 'SULFOPYRUVATE', '543': 'CALCIUM ION, 6 WATERS PLUS ETHANOL COORDINATED',
    'SF4': 'Iron Sulfur cluster', '7QM': 'MENAQUINONE- 7 (ALTERED)', 'MQ7': 'MENAQUINONE- 7',
    'MQ9': 'MENAQUINONE 9', 'MQ8': 'MENAQUINONE 8', 'PQQ': 'PYRROLOQUINOLINE QUINONE',
    'UQ1': 'UBIQUINONE-1', 'UQ2': 'UBIQUINONE-2', 'UQ7': 'UBIQUINONE-7',
    'EGL': 'Ethylene Glycol', 'BEF': 'BERYLLIUM TRIFLUORIDE ION', 'DMS': 'Dimethyl Sulfoxide',
    'NO2': 'NITRITE ION', 'NO3': 'NITRATE ION', 'CIT': 'Citric Acid',
    'PI': 'HYDROGENPHOSPHATE ION', 'S': 'sulfur', 'MD': 'molybdopterin guanine dinucleotide',
    'LCP': 'PERCHLORATE ION', 'MOO': 'MOLYBDATE ION', 'LCO': 'CHLORATE ION',
    'COS': 'Coenzyme A persulfide', 'PEP': 'PHOSPHOENOLPYRUVATE', 'PER': 'PEROXIDE ION',
    'COA': 'COENZYME A', 'SMM': 'S-ADENOSYLMETHIONINE METHYL ESTER', 'BCT': 'Bicarbonate',
    'AUC': 'GOLD CYANIDE', 'FMN': 'FMN', 'HEA': 'HEM A', 'HEB': 'HEM B', 'HEC': 'HEM C',
    'NDP': 'NADPH', 'HEO': 'Heme o', 'CO3': 'Carbonate Ion', 'PLP': 'PYRIDOXAL- 5- PHOSPHATE',
    'ALF': 'TETRAFLUOROALUMINATE ION', 'EDO': '1,2-Ethanediol', 'FS1': 'IRON/SULFUR CLUSTER',
    'FS3': 'Fe3-S4 Cluster', 'FS4': 'Iron/Sulfur Cluster', 'DPM': 'Dipyrromethane',
    'SOH': 'HYDROGEN SULFATE', 'PPR': 'PHOSPHONOPYRUVATE', 'GAI': 'Guanidine', 'FAD': 'FAD',
    'SO4': 'Sulfate Ion', 'NH2': 'Amino group', 'ARS': 'Arsenic', 'NAD': 'NAD', 'NAH': 'NAD',
    'NAI': 'NADH', 'NAP': 'NADPH', 'GHA': 'MIO', 'THF': '5-HYDROXYMETHYLENE-6-HYDROFOLIC ACID',
    'GSH': 'Glutathione', 'SAM': 'S-ADENOSYLMETHIONINE', 'BF4': 'BERYLLIUM TETRAFLUORIDE ION',
    'TPQ': '5-(2-CARBOXY-2-AMINOETHYL)-2-HYDROXY-1,4-BENZOQUINONE', 'TPP': 'Thiamine diphosphate',
    'LPA': 'LIPOIC ACID', 'OH': 'HYDROXIDE ION', 'PYR': 'Pyruvic acid', 'BO4': 'BORATE ION',
    'FLC': 'CITRATE ANION', 'FUM': 'Fumaric acid', 'PDP': 'PYRIDOXAL- 5- DIPHOSPHATE',
    'BME': 'Betamercaptoethanol', 'WO4': 'Tungstate ion', 'YT3': 'Yttrium ion', 'DHE': 'Heme d',
    'OXY': 'Oxygen molecule', '2HP': 'DIHYDROGENPHOSPHATE ION', 'PO3': 'PHOSPHITE ION',
    'OXL': 'OXALATE ION', 'TRS': '2-AMINO-2-HYDROXYMETHYL-PROPANE-1,3-DIOL', 'BTN': 'Biotin',
    'ACE': 'Acetyl Group', 'SUL': 'Sulfate Anion', 'B12': 'Cobalamin',
    'ATP': 'ADENOSINE-5-TRIPHOSPHATE', 'ACT': 'Acetate Ion', 'U10': 'UBIQUINONE-10',
    'ACY': 'Acetic Acid', 'ASC': 'ASCORBIC ACID', 'SEO': 'Beta-Mercaptoethanol 102LA4',
    'ADP': 'ADENOSINE-5-DIPHOSPHATE', 'UNX': 'Unknown', 'CUZ': '(MU-4-SULFIDO)-TETRA-NUCLEAR COPPER ION',
    'XE': 'Xenon', 'GTP': 'Guanosine-5-triphosphate', 'GTT': 'Glutathione',
    'FES': 'Fe2/S2 (Inorganic) Cluster', 'FRU': 'Fructose', 'GMP': 'Guanosine',
    'H2S': 'HYDROSULFURIC ACID', 'AZI': 'AZIDE ION', 'BPV': 'Bromopyruvate',
    'MSE': 'Selenmethionine', 'BIO': 'BIOPTERIN', 'CYN': 'CYANIDE ION', 'GOL': 'Glycerol',
    'PO4': 'Phosphate Ion', 'MLI': 'MALONATE ION', 'MLT': 'MALATE ION', 'RBF': 'Riboflavin',
    'SRM': 'Siroheme', 'UQ8': 'UBIQUINONE-8', 'HOH': 'Water', 'NAG': 'n-acetylglucosamine',
    'EPE': 'Hepes', 'PEG': 'poly(ethylene glycol)', '0U': 'L-nucleotide', '6HT': 'L-nucleotide',
    '0C': 'L-nucleotide', '0G': 'L-nucleotide',
    'PG0': 'Imported from Twilight', 'ETX': 'Imported from Twilight',
    'P4C': 'Imported from Twilight', '1PG': 'Imported from Twilight',
    'DOD': 'Imported from Twilight', '1PE': 'Imported from Twilight',
    'BEZ': 'Imported from Twilight', '2PE': 'Imported from Twilight',
    'TME': 'Imported from Twilight', 'P6G': 'Imported from Twilight',
    '3PO': 'Imported from Twilight', '7PE': 'Imported from Twilight',
    'DTT': 'Imported from Twilight', 'NHE': 'Imported from Twilight',
    'PGE': 'Imported from Twilight', 'PGO': 'Imported from Twilight',
    '8PE': 'Imported from Twilight', '12P': 'Imported from Twilight',
    '13P': 'Imported from Twilight', 'PGR': 'Imported from Twilight',
    '9PE': 'Imported from Twilight', 'P33': 'Imported from Twilight',
    '3PG': 'Imported from Twilight', 'TLA': 'Imported from Twilight',
    'MLA': 'Imported from Twilight', 'IMD': 'Imported from Twilight',
    'FMT': 'Imported from Twilight', 'HEZ': 'Imported from Twilight',
    'ACN': 'Imported from Twilight', 'PDO': 'Imported from Twilight',
    'PG4': 'Imported from Twilight', 'PG5': 'Imported from Twilight',
    'NH3': 'Imported from Twilight', 'MRD': 'Imported from Twilight',
    '15P': 'Imported from Twilight', 'MPD': 'Imported from Twilight',
    'IPA': 'Imported from Twilight', 'PE8': 'Imported from Twilight',
    'PE9': 'Imported from Twilight', 'HTO': 'Imported from Twilight',
    'PE5': 'Imported from Twilight', 'PE6': 'Imported from Twilight',
    'PE7': 'Imported from Twilight', 'TBU': 'Imported from Twilight',
    'PE1': 'Imported from Twilight', 'PE2': 'Imported from Twilight',
    'PE3': 'Imported from Twilight', 'CE1': 'Imported from Twilight',
    'PE4': 'Imported from Twilight', 'MOH': 'Imported from Twilight',
    'PG6': 'Imported from Twilight', 'MES': 'Imported from Twilight',
    'P22': 'Imported from Twilight',
}

DEFAULT_METALS = dict(metals)
DEFAULT_LIGAND_BLACKLIST = dict(ligand_blacklist)


def get_default_entries():
    """Return the built-in blacklist/metal entries as a flat list, for the
    frontend to render as togglable checkboxes.

    Returns:
        list: A list of ``{"code": str, "name": str, "category": str}``
        dicts, one per entry in :data:`DEFAULT_LIGAND_BLACKLIST` (category
        ``"blacklist"``) and :data:`DEFAULT_METALS` (category ``"metal"``).
        The empty-string entry in ``DEFAULT_LIGAND_BLACKLIST`` (used
        internally for unnamed components) is skipped, since it has nothing
        meaningful to show or toggle in the UI.
    """
    entries = []
    for code, name in DEFAULT_LIGAND_BLACKLIST.items():
        if not code:
            continue
        entries.append({"code": code, "name": name, "category": "blacklist"})
    for code, name in DEFAULT_METALS.items():
        entries.append({"code": code, "name": name, "category": "metal"})
    entries.sort(key=lambda e: (e["category"], e["code"]))
    return entries


def build_effective_lists(disabled_codes=None, custom_entries=None, replace=None):
    """Compute the metals/ligand_blacklist dicts to use for one analysis run.

    Builds fresh dict instances rather than mutating the module-level
    ``metals``/``ligand_blacklist`` globals, so concurrent requests/jobs
    never interfere with each other's blacklist customization.

    Args:
        disabled_codes (iterable, optional): Component codes (e.g.
            ``"HOH"``, ``"MG"``) to remove from the built-in defaults before
            analysis, i.e. entries the user unchecked in the UI. Ignored
            when ``replace`` is given, since starting from an uploaded file
            already means "use exactly these entries". Case-insensitive.
        custom_entries (iterable, optional): Extra entries to add on top,
            each a dict with keys ``"code"``, ``"name"`` (optional), and
            ``"category"`` (``"blacklist"`` or ``"metal"``, defaults to
            ``"blacklist"``). These are applied last, so a custom entry can
            reintroduce a code that was disabled/omitted above.
        replace (dict, optional): If given with a non-empty ``"metals"``
            and/or ``"ligand_blacklist"`` key, those dicts fully replace the
            built-in defaults as the starting point (e.g. from an uploaded
            file), instead of starting from :data:`DEFAULT_METALS` /
            :data:`DEFAULT_LIGAND_BLACKLIST`.

    Returns:
        tuple: ``(metals_dict, ligand_blacklist_dict)`` — new dict
        instances safe for a single caller to use and mutate freely.

    Raises:
        None
    """
    disabled = {str(c).strip().upper() for c in (disabled_codes or []) if str(c).strip()}

    if replace and (replace.get("metals") or replace.get("ligand_blacklist")):
        eff_metals = {str(k).strip().upper(): v for k, v in (replace.get("metals") or {}).items()}
        eff_blacklist = {str(k).strip().upper(): v for k, v in (replace.get("ligand_blacklist") or {}).items()}
    else:
        eff_metals = dict(DEFAULT_METALS)
        eff_blacklist = dict(DEFAULT_LIGAND_BLACKLIST)
        for code in disabled:
            eff_metals.pop(code, None)
            eff_blacklist.pop(code, None)

    for entry in (custom_entries or []):
        code = str(entry.get("code", "")).strip().upper()
        if not code:
            continue
        name = str(entry.get("name", "")).strip() or code
        if entry.get("category") == "metal":
            eff_metals[code] = name
        else:
            eff_blacklist[code] = name

    return eff_metals, eff_blacklist


def parse_uploaded_list(text):
    """Parse a user-uploaded blacklist file into metals/blacklist dicts.

    Supports two formats:

    1. The native :func:`dump_lists` CSV format, with ``[Blacklist]`` and
       ``[Non-propagating]`` section headers followed by ``code,name`` rows
       (this is exactly what "Load from file…" -> a file previously saved
       with :func:`dump_lists` looks like).
    2. A simple format with no section headers: one entry per line, either
       ``CODE`` alone or ``CODE,Name`` / ``CODE;Name`` / ``CODE<tab>Name``.
       Every entry parsed this way is treated as a ``"blacklist"`` category
       entry, since that's what most users mean by "ligands to ignore";
       codes meant to be treated as metals can still be added individually
       via the category selector in the UI, or included under an explicit
       ``[Non-propagating]`` header.

    Blank lines and lines starting with ``#`` are ignored in both formats.

    Args:
        text (str): Raw file content to parse.

    Returns:
        tuple: ``(metals_dict, ligand_blacklist_dict)``. Both are empty
        dicts if ``text`` contains no parseable rows.

    Raises:
        None
    """
    new_m, new_lb = {}, {}
    has_sections = "[Blacklist]" in text or "[Non-propagating]" in text
    d = new_lb if not has_sections else None

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line == "[Blacklist]":
            d = new_lb
            continue
        if line == "[Non-propagating]":
            d = new_m
            continue
        if d is None:
            continue

        for sep in (",", ";", "\t"):
            if sep in line:
                parts = [p.strip() for p in line.split(sep, 1)]
                break
        else:
            parts = [line, ""]

        code = parts[0].strip().upper()
        if not code:
            continue
        name = parts[1].strip() if len(parts) > 1 and parts[1].strip() else code
        d[code] = name

    return new_m, new_lb


def update_lists(new_m=metals, new_lb=ligand_blacklist):
    """Replace the module-level cofactor lookup dictionaries.

    Overwrites the global ``metals`` and ``ligand_blacklist`` dictionaries
    with the provided replacements.

    Args:
        new_m (dict, optional): New mapping of metal/ion codes to their
            descriptive names. Defaults to the current ``metals`` dict.
        new_lb (dict, optional): New mapping of blacklisted ligand codes to
            their descriptive names. Defaults to the current
            ``ligand_blacklist`` dict.

    Returns:
        None: This function does not return a value; it mutates the
        module-level globals in place.

    Raises:
        None
    """
    global metals, ligand_blacklist
    metals = new_m
    ligand_blacklist = new_lb


def dump_lists(fname='notligands'):
    """Write the ``ligand_blacklist`` and ``metals`` dictionaries to a CSV file.

    The CSV file is written with two sections, each preceded by a header
    row: ``[Blacklist]`` for the ``ligand_blacklist`` entries and
    ``[Non-propagating]`` for the ``metals`` entries. Each entry is written
    as a ``key, value`` row.

    Args:
        fname (str, optional): Path or base name of the output CSV file.
            If it does not already end with ``.csv``, that extension is
            appended automatically. Defaults to ``'notligands'``.

    Returns:
        int: Always returns ``0`` on successful completion.

    Raises:
        OSError: If the file cannot be opened or written to (e.g. due to
            invalid path or insufficient permissions).
    """
    if not fname.endswith('.csv'):
        fname += '.csv'
    with open(fname, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['[Blacklist]'])
        for k, v in ligand_blacklist.items():
            w.writerow([k, v])
        w.writerow(['[Non-propagating]'])
        for k, v in metals.items():
            w.writerow([k, v])
    return 0


def load_lists(fname):
    """Load cofactor lists from a CSV file and update the global dictionaries.

    Parses a CSV file structured with ``[Blacklist]`` and
    ``[Non-propagating]`` section headers, rebuilds the corresponding
    dictionaries from the rows found under each section, and applies them
    via :func:`update_lists`.

    Args:
        fname (str): Path to the CSV file to read. The file must contain
            section headers (``[Blacklist]``, ``[Non-propagating]``)
            followed by rows of ``key, value`` pairs.

    Returns:
        int: Always returns ``0`` on successful completion.

    Raises:
        OSError: If the file cannot be opened or read (e.g. it does not
            exist or is not accessible).
        csv.Error: If the file content cannot be parsed as valid CSV.
        IndexError: If a data row is malformed and does not contain both
            a key and a value column.
    """
    new_m, new_lb = {}, {}
    d = None
    with open(fname, 'r', newline='') as f:
        for row in csv.reader(f):
            if row:
                if row[0] == '[Blacklist]':
                    d = new_lb
                elif row[0] == '[Non-propagating]':
                    d = new_m
                elif d is not None:
                    d[row[0]] = row[1]
    update_lists(new_m, new_lb)
    return 0
