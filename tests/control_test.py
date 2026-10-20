"""
Validation of Pólya Pipeline over 32 aromatic molecules → tests/validation_report.txt

Usage (from the Pólya Pipeline folder):
    python tests/control_test.py

For every molecule it checks Pólya's count against the enumeration, an independent
brute-force count of distinct compounds, hand-calculated values, the IUPAC numbering,
the GUI drawing and the .xyz/.gjf files. Files are written to a temporary folder, so
generated_files/ is not touched. Takes about a minute.
"""
import contextlib, io, sys, tempfile, time, datetime, platform
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
SCRATCH = Path(tempfile.mkdtemp(prefix="polya_validation_"))
REPORT = Path(__file__).resolve().parent / "validation_report.txt"
sys.path.insert(0, str(PROJECT))

import numpy, rdkit, PIL
import derivatives as der
from rdkit import Chem

der.OUTPUT_DIR = SCRATCH / "generated_files"     # keep the real folder untouched
PREVIEW_DIR = SCRATCH / "previews"; PREVIEW_DIR.mkdir(exist_ok=True)
PREVIEW_SIZE, PREVIEW_BG = (300, 130), (1.0, 1.0, 1.0)

# (group, name, SMILES, substituents, max_subst, expected Pólya total, expected compounds)
# Expected values are worked out by hand (see notes at the end); None = not given.
CASES = [
    ("Monocyclic, 3 C",  "cyclopropenyl cation",   "[cH+]1cc1",             "N",   None, 4,   3),
    ("Monocyclic, 5 C",  "cyclopentadienyl anion", "[cH-]1cccc1",           "N",   None, 8,   8),
    ("Monocyclic, 6 C",  "benzene",                "c1ccccc1",              "N",   None, 13,  13),
    ("Monocyclic, 6 C",  "benzene (N,O+)",         "c1ccccc1",              "N,[O+]", None, 92, 92),
    ("Monocyclic, 7 C",  "tropylium cation",       "[cH+]1cccccc1",         "N",   None, 18,  17),
    ("Monocyclic, 18 C", "[18]annulene",           "c1ccccccccccccccccc1",  "N",   3,    None, None),
    ("Heteroaromatic",   "pyridine",               "n1ccccc1",              "N",   None, 20,  12),
    ("Heteroaromatic",   "pyridine (P)",           "n1ccccc1",              "P",   None, 20,  20),
    ("Heteroaromatic",   "pyrimidine",             "n1cnccc1",              "N",   None, 12,  9),
    ("Heteroaromatic",   "pyrimidine (P)",         "n1cnccc1",              "P",   None, 12,  12),
    ("Heteroaromatic",   "thiophene",              "s1cccc1",               "N",   None, 10,  10),
    ("Heteroaromatic",   "furan",                  "o1cccc1",               "N",   None, 10,  10),
    ("Heteroaromatic",   "pyrrole",                "[nH]1cccc1",            "N",   None, 10,  10),
    ("Substituted ring", "toluene",                "Cc1ccccc1",             "N",   None, 20,  20),
    ("Substituted ring", "toluene (CH3 numbered)", "[CH3:6]c1[cH:1][cH:2][cH:3][cH:4][cH:5]1", "N", None, 40, 40),
    ("Substituted ring", "phenol",                 "Oc1ccccc1",             "N",   None, 20,  20),
    ("Substituted ring", "chlorobenzene",          "Clc1ccccc1",            "N",   None, 20,  20),
    ("Substituted ring", "nitrobenzene",           "O=[N+]([O-])c1ccccc1",  "N",   None, 20,  20),
    ("Substituted ring", "benzonitrile",           "N#Cc1ccccc1",           "N",   None, 20,  20),
    ("Substituted ring", "p-xylene",               "Cc1ccc(C)cc1",          "N",   None, 7,   7),
    ("Substituted ring", "mesitylene",             "Cc1cc(C)cc(C)c1",       "N",   None, 4,   4),
    ("Substituted ring", "o-dichlorobenzene",      "Clc1ccccc1Cl",          "N",   None, 10,  10),
    ("Fused / polycyclic", "naphthalene",          "c1ccc2ccccc2c1",        "N",   None, 76,  76),
    ("Fused / polycyclic", "naphthalene (N,O+)",   "c1ccc2ccccc2c1",        "N,[O+]", None, 1701, 1701),
    ("Fused / polycyclic", "azulene",              "c1ccc2cccc2cc1",        "N",   None, None, None),
    ("Fused / polycyclic", "indole",               "c1ccc2[nH]ccc2c1",      "N",   None, 64,  64),
    ("Fused / polycyclic", "quinoline",            "c1ccc2ncccc2c1",        "N",   None, 128, None),
    ("Fused / polycyclic", "quinoline (P)",        "c1ccc2ncccc2c1",        "P",   None, 128, 128),
    ("Fused / polycyclic", "anthracene",           "c1ccc2cc3ccccc3cc2c1",  "N",   None, None, None),
    ("Fused / polycyclic", "phenanthrene",         "c1ccc2c(c1)ccc1ccccc12", "N",  None, None, None),
    ("Fused / polycyclic", "biphenyl",             "c1ccc(cc1)-c1ccccc1",   "N",   None, None, None),
    ("Fused / polycyclic", "pyrene",               "c1cc2ccc3cccc4ccc(c1)c2c34", "N", None, None, None),
]
# IUPAC numbering expected for the parent: (position labels, locants of the 1-substituted
# derivatives). Locants of the monosubstituted rows are the textbook ones (e.g. naphthalene
# 1 = α, 2 = β; anthracene 1, 2, 9).
NUMBERING = {
    "benzene":           ("1 2 3 4 5 6", {"1"}),
    "pyridine (P)":      ("2 3 4 5 6", {"2", "3", "4"}),
    "pyrimidine (P)":    ("2 4 5 6", {"2", "4", "5"}),
    "thiophene":         ("2 3 4 5", {"2", "3"}),
    "furan":             ("2 3 4 5", {"2", "3"}),
    "pyrrole":           ("2 3 4 5", {"2", "3"}),
    "toluene":           ("2 3 4 5 6", {"2", "3", "4"}),
    "phenol":            ("2 3 4 5 6", {"2", "3", "4"}),
    "o-dichlorobenzene": ("3 4 5 6", {"3", "4"}),
    "p-xylene":          ("2 3 5 6", {"2"}),
    "mesitylene":        ("2 4 6", {"2"}),
    "naphthalene":       ("1 2 3 4 5 6 7 8", {"1", "2"}),
    "azulene":           ("1 2 3 4 5 6 7 8", {"1", "2", "4", "5", "6"}),
    "indole":            ("2 3 4 5 6 7", {"2", "3", "4", "5", "6", "7"}),
    "quinoline (P)":     ("2 3 4 5 6 7 8", {"2", "3", "4", "5", "6", "7", "8"}),
    "anthracene":        ("1 2 3 4 5 6 7 8 9 10", {"1", "2", "9"}),
    "phenanthrene":      ("1 2 3 4 5 6 7 8 9 10", {"1", "2", "3", "4", "9"}),
    "pyrene":            ("1 2 3 4 5 6 7 8 9 10", {"1", "2", "4"}),
    "biphenyl":          ("2 2′ 3 3′ 4 4′ 5 5′ 6 6′", {"2", "3", "4"}),
    "[18]annulene":      (" ".join(str(i) for i in range(1, 19)), {"1"}),
}

XYZ_LIMIT = 100       # write .xyz/.gjf files only when the run has at most this many derivatives
CUSTOM_ROUTE = {"benzene (N,O+)": "#m062x/def2tzvp opt freq"}   # other runs: default route


def check_gjf(name, subst, rows, route):
    """Checks every .gjf of a run: layout, charge/mult, atoms and same geometry as the .xyz."""
    folder = der.safe_filename(name)
    tag = der.safe_filename("-".join(s.strip() for s in subst.split(",")))
    problems, n = [], 0
    for row in rows:
        smi = row[-1]
        if "invalid" in smi:
            continue
        stem = f"{folder}_{tag}_{row[0]:03d}"
        gjf = der.OUTPUT_DIR / "gaussian" / folder / f"{stem}.gjf"
        xyz = der.OUTPUT_DIR / "xyz" / folder / f"{stem}.xyz"
        if not gjf.exists():
            problems.append(f"{stem}: missing"); continue
        text = gjf.read_text(encoding="utf-8")
        L = text.split("\n")
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        q = Chem.GetFormalCharge(mol)
        m = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) + 1
        atoms = L[5:5 + mol.GetNumAtoms()]
        ok = (L[0] == route and L[1] == "" and L[2].startswith(stem) and smi in L[2]
              and L[3] == "" and L[4] == f"{q} {m}"
              and all(len(a.split()) == 4 for a in atoms)
              and sorted(a.split()[0] for a in atoms) == sorted(at.GetSymbol() for at in mol.GetAtoms())
              and L[5 + mol.GetNumAtoms()] == "" and text.endswith("\n\n"))
        if ok and xyz.exists():
            xl = xyz.read_text(encoding="utf-8").splitlines()[2:]
            ok = all(abs(float(a) - float(b)) < 1e-4
                     for g, x in zip(atoms, xl) for a, b in zip(g.split()[1:], x.split()[1:]))
        if not ok:
            problems.append(f"{stem}: wrong content")
        n += 1
    return n, problems


def constitution_key(smi):
    """Resonance-independent identity: elements + connectivity + H count per atom.

    Two SMILES of the same ion that only differ in where the charge is drawn give the
    same key; different constitutional isomers or tautomers never do.
    """
    m = Chem.MolFromSmiles(smi)
    rw = Chem.RWMol(m)
    total_charge = Chem.GetFormalCharge(m)
    for a in rw.GetAtoms():
        h = a.GetTotalNumHs()
        a.SetIsotope(h + 1)                 # H count as a label, always written
        a.SetAtomMapNum(0)
        a.SetFormalCharge(0); a.SetIsAromatic(False)
        a.SetNoImplicit(True); a.SetNumExplicitHs(0)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    return Chem.MolToSmiles(rw, canonical=True) + f" q={total_charge}"


def brute_force(smiles, subst_text, max_subst):
    """Independent count: distinct compounds over EVERY substitution pattern (no symmetry)."""
    mol = Chem.MolFromSmiles(smiles)
    idx_pos = [i for _, i in der.substitution_positions(mol)[0]]   # same positions as the program
    subst = [s.strip() for s in subst_text.split(",")]
    het = [der.read_heteroatom(s) for s in subst]
    unique = set()
    for p in der.candidate_patterns(len(idx_pos), len(subst), max_subst):
        smi = der.build_smiles(mol, idx_pos, p, het)
        if smi is not None:
            unique.add(constitution_key(smi))
    return len(unique)


def preview_ok(smiles, name):
    img = der.smiles_image(smiles, size=PREVIEW_SIZE, background=PREVIEW_BG).convert("RGB")
    img.save(PREVIEW_DIR / f"{der.safe_filename(name)}.png")
    bg = tuple(round(c * 255) for c in PREVIEW_BG)
    drawn = sum(1 for px in (img.get_flattened_data() if hasattr(img, "get_flattened_data") else img.getdata()) if max(abs(a - b) for a, b in zip(px, bg)) > 30)
    return img.size == PREVIEW_SIZE and drawn > 200


summary, logs = [], []
for group, name, smi, subst, max_subst, exp_polya, exp_comp in CASES:
    t0 = time.perf_counter()
    rec = dict(group=group, name=name, smi=smi, subst=subst, max=max_subst,
               exp_polya=exp_polya, exp_comp=exp_comp)
    buf = io.StringIO()
    try:
        rec["preview"] = preview_ok(smi, name)
        with contextlib.redirect_stdout(io.StringIO()):
            rows = der.run_pipeline(smi, subst, name=der.safe_filename(name),
                                    max_subst=max_subst, generate_xyz_files=False)
        write_xyz = len(rows) <= XYZ_LIMIT
        route = CUSTOM_ROUTE.get(name, der.DEFAULT_ROUTE)
        with contextlib.redirect_stdout(buf):
            rows = der.run_pipeline(smi, subst, name=der.safe_filename(name),
                                    max_subst=max_subst, generate_xyz_files=write_xyz,
                                    generate_gaussian_files=write_xyz, gaussian_route=route)
        out = buf.getvalue()
        lines = out.splitlines()
        rec["rows"] = len(rows)
        rec["invalid"] = sum("invalid" in r[-1] for r in rows)
        valid = rec["rows"] - rec["invalid"]
        rec["distinct"] = len({constitution_key(r[-1]) for r in rows if "invalid" not in r[-1]})
        totals = [l for l in lines if l.strip().startswith("Total")]
        rec["polya"] = int(totals[0].split("|")[1])
        # patterns enumerated = rows before removing duplicates
        rec["patterns"] = int(totals[1].split("|")[1]) if len(totals) > 1 else rec["rows"]
        rec["dups"] = rec["patterns"] - rec["rows"]
        rec["auto_num"] = any("numbered automatically" in l for l in lines)
        rec["ci"] = ("✓" if "✓ Cycle indices match" in out else
                     "ℹ" if "less symmetric than the graph" in out else
                     "skip" if "[3D]     skipped" in out else "⚠")
        g3 = next(l for l in lines if "[3D]" in l)
        gg = next(l for l in lines if "[Graph]" in l)
        rec["G3d"] = g3.split("|G| =")[1].strip() if "|G|" in g3 else "—"
        rec["ngr"], rec["Ggr"] = gg.split("n =")[1].split()[0], gg.split("|G| =")[1].strip()
        rec["warnings"] = sum("⚠" in l for l in lines)
        rec["xyz"] = ("—" if not write_xyz else
                      "FAIL" if "FAILED" in out else
                      next(l for l in lines if "geometries written" in l).split("✓")[1].split()[0])
        if write_xyz:
            n_gjf, gjf_problems = check_gjf(der.safe_filename(name), subst, rows, route)
            rec["gjf"] = str(n_gjf) if not gjf_problems else "FAIL"
            rec["gjf_problems"] = gjf_problems
        else:
            rec["gjf"] = "—"
        rec["brute"] = brute_force(smi, subst, max_subst)
        pos, rec["scheme"] = der.substitution_positions(Chem.MolFromSmiles(smi))
        rec["labels"] = " ".join(str(l) for l, _ in pos)
        rec["mono"] = sorted({r[-3].split(": ")[1] for r in rows if r[1] == 1},
                             key=der.locant_key)
        if name in NUMBERING:
            exp_labels, exp_mono = NUMBERING[name]
            rec["num_ok"] = (rec["scheme"] == "iupac" and rec["labels"] == exp_labels
                             and set(rec["mono"]) == exp_mono)
        else:
            rec["num_ok"] = None
        checks = {
            "preview":   rec["preview"],
            "polya":     max_subst is not None or rec["polya"] == rec["patterns"],
            "unique":    rec["distinct"] == valid,
            "brute":     rec["brute"] == valid,
            "exp_polya": exp_polya is None or exp_polya == rec["polya"],
            "exp_comp":  exp_comp is None or exp_comp == valid,
            "xyz":       rec["xyz"] != "FAIL",
            "gjf":       rec["gjf"] != "FAIL" and rec["gjf"] in ("—", rec["xyz"]),
            # a warning is only expected for impossible-valence derivatives
            "warnings":  rec["warnings"] == (1 if rec["invalid"] else 0),
            "numbering": rec["num_ok"] is not False,
        }
        failed = [k for k, ok in checks.items() if not ok]
        rec["status"] = "PASS" if not failed else "FAIL (" + ", ".join(failed) + ")"
    except Exception as e:
        out = buf.getvalue() + f"\n  EXCEPTION: {type(e).__name__}: {e}\n"
        rec["status"] = f"ERROR ({type(e).__name__}: {e})"
    rec["time"] = time.perf_counter() - t0
    summary.append(rec)
    logs.append((name, smi, subst, max_subst, out))
    print(f"{rec['status']:5}  {name}  ({rec['time']:.1f} s)", flush=True)


# Invalid route must be rejected before anything is computed
try:
    with contextlib.redirect_stdout(io.StringIO()):
        der.run_pipeline("c1ccccc1", "N", name="route_test", generate_xyz_files=False,
                         generate_gaussian_files=True, gaussian_route="b3lyp/6-31G opt")
    ROUTE_TEST = "NOT rejected  ✗"
except ValueError as e:
    ROUTE_TEST = f"rejected ✓ ({e})"

def gjf_rel(name, subst, num):
    folder = der.safe_filename(der.safe_filename(name))
    tag = der.safe_filename("-".join(s.strip() for s in subst.split(",")))
    return f"{folder}/{folder}_{tag}_{num:03d}.gjf"

GJF_EXAMPLES = [   # derivative 2 = one N, 3 = one [O+] (rows sorted by substituents)
    ("benzene + 1 N (default route)", gjf_rel("benzene", "N", 2)),
    ("benzene + 1 [O+], charge +1 (custom route)", gjf_rel("benzene (N,O+)", "N,[O+]", 3)),
    ("tropylium cation + 1 N, charge +1", gjf_rel("tropylium cation", "N", 2)),
]

NOTES = """\
  1. Duplicate compounds (pyridine, pyrimidine and quinoline with N).
     Pólya's theorem counts substitution PATTERNS under the symmetry of the parent molecule. When the
     substituent is the same element as a heteroatom of the ring, different patterns can give the
     same compound (e.g. three non-equivalent ways of putting two N on pyridine all give
     1,2,4-triazine, c1cnncn1), because the "old" and the "new" N cannot be told apart.
     The program now removes these duplicates and prints both numbers (Pólya patterns and distinct
     compounds, per number of substituents). Column Dup. = patterns removed. The brute-force count
     confirms that what is left are exactly the different compounds. With P as substituent the same
     rings have no duplicates.

  2. Automatic numbering of the ring positions (toluene, p-xylene, mesitylene).
     CH groups outside the aromatic rings (CH3) are no longer substitution positions: the program
     numbers the aromatic CH itself (column Auto = ✓), prints the numbered SMILES and, in the GUI,
     writes it into the SMILES field. Toluene therefore gives 20 derivatives (5 ring positions).
     To substitute the CH3 as well, number it in the SMILES: run "toluene (CH3 numbered)" gives 40
     and its 3D check is skipped, because the CH3 position carries 3 H.

  3. 3D check (column CI). It now compares only the H atoms of the substitution positions, so the
     H of CH3, OH and NH no longer cause a mismatch (phenol, pyrrole, toluene, indole: ✓).
       ✓     3D and graph cycle indices are the same.
       ℹ     the optimised 3D geometry is less symmetric than the graph (twisted biphenyl,
             [18]annulene, azulene, charge localised by the force field in cyclopropenyl⁺ and
             cyclopentadienyl⁻). Only a note: the counting uses the graph symmetry.
       skip  a position carries more than one H.
       ⚠     the 3D symmetry does not fit the graph at all — did not happen in any run.
     Column ⚠ = warning lines printed by the program. The only warning left is for derivatives
     with impossible valence (all-N cyclopropenyl and tropylium), which must be reported.

  4. Charged rings (cyclopropenyl⁺, cyclopentadienyl⁻, tropylium⁺). The SMILES puts the charge
     on one atom, so the same ion can be written in several ways (resonance forms). With "Respect
     formal charges" off (default) all ring CH are equivalent, which is the chemically correct
     choice, and the counts match the hand calculation. The brute-force check compares compounds
     by constitution (elements, connectivity, H count), so resonance forms count once.

  5. Hand-calculated values (Pólya, 2 colours unless noted):
       benzene D6:            (2^6 + 2·2 + 2·2^2 + 2^3 + 3·2^4 + 3·2^3)/12 = 13
       benzene N,[O+] (3 col.): (3^6 + 2·3 + 2·3^2 + 3^3 + 3·3^4 + 3·3^3)/12 = 92
       naphthalene D2h:       (2^8 + 3·2^4)/4 = 76          N,[O+]: (3^8 + 3·3^4)/4 = 1701
       cyclopropenyl D3:      (2^3 + 2·2 + 3·2^2)/6 = 4     (3 compounds: all-N is impossible)
       cyclopentadienyl D5:   (2^5 + 4·2 + 5·2^3)/10 = 8
       tropylium D7:          (2^7 + 6·2 + 7·2^4)/14 = 18   (17 compounds: all-N is impossible)
       pyridine C2v (5 pos.): (2^5 + 2^3)/2 = 20     pyrimidine (4 pos.): (2^4 + 2^3)/2 = 12
         with N they are the diazines/triazines/…: 12 and 9 different compounds
       thiophene / furan / pyrrole / o-dichlorobenzene (4 pos.): (2^4 + 2^2)/2 = 10
       toluene, phenol, chlorobenzene, nitrobenzene, benzonitrile (5 pos.): (2^5 + 2^3)/2 = 20
       toluene with the CH3 numbered (6 pos.): (2^6 + 2^4)/2 = 40
       p-xylene (4 ring pos., D2h): (2^4 + 3·2^2)/4 = 7
       mesitylene (3 ring pos., D3h): (2^3 + 2·2 + 3·2^2)/6 = 4
       indole (6 pos., no symmetry): 2^6 = 64     quinoline (7 pos.): 2^7 = 128

  6. The GUI preview images (check "Prev.") were rendered at 300x130 px, as in the GUI.
     The .xyz, .gjf and CSV files of this test were written to a temporary folder, so the
     generated_files folder of the project was not modified.""".splitlines()


# ── Write the report ───────────────────────────────────────────────────────────
L = []
w = L.append
w("=" * 120)
w("  PÓLYA PIPELINE — VALIDATION REPORT")
w("=" * 120)
w(f"  Date:     {datetime.datetime.now():%Y-%m-%d %H:%M}")
w(f"  Python {platform.python_version()} · numpy {numpy.__version__} · rdkit {rdkit.__version__} · pillow {PIL.__version__}")
w(f"  Program:  Pólya Pipeline {der.__version__} — derivatives.py (run_pipeline) + GUI drawing (smiles_image)")
w("")
w("  What is checked for each molecule (PASS = all of them hold):")
w("    preview    the 2D drawing of the GUI renders (correct size, structure actually drawn).")
w("    polya      the total predicted by Pólya's theorem = number of patterns enumerated")
w("               (skipped when Max. subst. is used, since Pólya then counts all degrees).")
w("    unique     every derivative listed is a different compound (duplicates were removed).")
w("    brute      independent count: every possible substitution pattern is built and the different")
w("               compounds are counted, without using any symmetry. Must equal the derivatives listed.")
w("    exp_polya  Pólya total = value worked out by hand (where given).")
w("    exp_comp   derivatives listed = number of different compounds worked out by hand (where given).")
w("    xyz        3D geometries generated for every derivative (only for runs with ≤ 100 derivatives).")
w("    warnings   a ⚠ warning is printed only when it is needed (impossible valence).")
w("    numbering  the positions get the IUPAC locants of the parent, and every monosubstituted")
w("               derivative is listed with its textbook locant (where an expectation is given).")
w("")
w("-" * 120)
w("  SUMMARY")
w("-" * 120)
hdr = (f"  {'Molecule':<24}{'Subst.':<8}{'Max':>4}{'Auto':>5}{'n':>4}{'|G| gr':>8}{'|G| 3D':>8}{'CI':>6}"
       f"{'Pólya':>7}{'Dup.':>6}{'Found':>7}{'Brute':>7}{'Exp.':>7}{'Inval.':>7}{'⚠':>3}  {'XYZ':>5}{'GJF':>5}{'Prev.':>6}  {'Result'}")
w(hdr); w("  " + "─" * (len(hdr) - 2))
current = None
for r in summary:
    if r["group"] != current:
        current = r["group"]; w(f"  [{current}]")
    if r["status"].startswith("ERROR"):
        w(f"  {r['name']:<24}{r['subst']:<8}  {r['status']}"); continue
    w(f"  {r['name']:<24}{r['subst']:<8}{(r['max'] or '—')!s:>4}{'✓' if r['auto_num'] else '':>5}{r['ngr']:>4}"
      f"{r['Ggr']:>8}{r['G3d']:>8}{r['ci']:>6}{r['polya']:>7}{r['dups']:>6}{r['rows']:>7}{r['brute']:>7}"
      f"{(r['exp_comp'] or '—')!s:>7}{r['invalid']:>7}{r['warnings']:>3}  {r['xyz']:>5}{r['gjf']:>5}{'✓' if r['preview'] else '✗':>6}  {r['status']}")
n_pass = sum(r["status"] == "PASS" for r in summary)
w("")
w(f"  {n_pass} / {len(summary)} runs PASS    ·    total time {sum(r['time'] for r in summary):.1f} s")
w("")
w("  Columns: Auto = ring positions numbered automatically; n = substitution positions; |G| = order of")
w("  the symmetry group (graph / 3D geometry); CI = 3D check (see note 3); Pólya = total predicted by")
w("  Pólya's theorem; Dup. = patterns removed because they repeat a compound; Found = derivatives")
w("  listed; Brute = different compounds by brute force; Exp. = different compounds by hand;")
w("  Inval. = derivatives with impossible valence; ⚠ = warning lines; XYZ = .xyz files written;")
w("  GJF = Gaussian 16 input files written and checked (see the section below).")
w("  With Max = 3, Pólya still shows the total for all numbers of substituents.")
w("")
w("-" * 120)
w("  GAUSSIAN 16 INPUT FILES (.gjf)")
w("-" * 120)
w("  Generated together with the .xyz files (runs with ≤ 100 derivatives). Every file was checked for:")
w("    · line 1 = route, line 2 blank, line 3 = title (file name + SMILES), line 4 blank;")
w("    · line 5 = 'charge multiplicity', equal to the formal charge and spin of the derivative's SMILES;")
w("    · one line per atom (symbol x y z), the same atoms as the molecule with its H atoms;")
w("    · coordinates identical to the .xyz file of the same derivative (same geometry);")
w("    · a blank line after the coordinates and at the end of the file (required by Gaussian).")
n_gjf = sum(int(r["gjf"]) for r in summary if r.get("gjf", "—").isdigit())
bad = [p for r in summary for p in r.get("gjf_problems", [])]
w(f"  Result: {n_gjf} .gjf files checked, {len(bad)} with problems." + ("" if not bad else "  " + "; ".join(bad[:10])))
w(f"  Default route: {der.DEFAULT_ROUTE}   ·   custom route tested in 'benzene (N,O+)': "
  f"{CUSTOM_ROUTE['benzene (N,O+)']}")
w(f"  Route without '#' (e.g. 'b3lyp/6-31G opt'): {ROUTE_TEST}")
for label, rel in GJF_EXAMPLES:
    p = der.OUTPUT_DIR / "gaussian" / rel
    w("")
    w(f"  Example — {label}:  generated_files/gaussian/{rel}")
    w("  " + "┄" * 70)
    for line in p.read_text(encoding="utf-8").split("\n"):
        w("  │ " + line)
    w("  " + "┄" * 70)
w("")
w("-" * 120)
w("  POSITION NUMBERING (IUPAC)")
w("-" * 120)
w("  Positions are numbered with the IUPAC locants of the parent ring system: monocycles by the")
w("  lowest-locant rules (senior heteroatom = 1, then heteroatoms, then existing substituents),")
w("  fused parents from a table of standard numberings, separate ring systems with primes (2′…).")
w("  A ring system outside the table is reported as 'program' numbering instead of guessing.")
w("  Each derivative is listed with the pattern that gives its substituents the lowest locants.")
w("")
w(f"  {'Molecule':<24}{'Numbering':<10}{'Positions':<44}{'1-substituted at':<22}Check")
w("  " + "─" * 108)
for r in summary:
    if "labels" not in r:
        continue
    mark = "—" if r["num_ok"] is None else ("✓" if r["num_ok"] else "✗")
    w(f"  {r['name']:<24}{r['scheme']:<10}{r['labels'][:43]:<44}{', '.join(r['mono'])[:21]:<22}{mark}")
w("")
w("  Extra check: 1,5-naphthyridine (c1cnc2cccnc2c1) is not in the table → "
  f"'{der.substitution_positions(Chem.MolFromSmiles('c1cnc2cccnc2c1'))[1]}' numbering, as intended.")
w("")
w("-" * 120)
w("  NOTES")
w("-" * 120)
for line in NOTES:
    w(line)
w("")
w("=" * 120)
w("  FULL PROGRAM OUTPUT FOR EACH RUN")
w("=" * 120)
for name, smi, subst, max_subst, out in logs:
    w("")
    w(f"### {name}   SMILES = {smi}   substituents = {subst}   max. subst. = {max_subst or 'all'}")
    w(out.rstrip())
REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")
print(f"\n{sum(r['status'] == 'PASS' for r in summary)} / {len(summary)} runs PASS  ·  report: {REPORT}")
