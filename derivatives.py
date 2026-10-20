#!/usr/bin/env python3
"""
derivatives.py — Pólya Pipeline engine
======================================
Counts and generates the derivatives of an aromatic molecule obtained by replacing
ring CH groups with other atoms, using Pólya's enumeration theorem:
  SMILES → 3D geometry → symmetry → cycle index → counting → derivative SMILES → XYZ / Gaussian

The only input is the SMILES of the parent molecule. The 3D geometry is generated
internally with ETKDG (RDKit) and is used only for symmetry detection.

Results are printed to the screen as a table and saved as CSV in
generated_files/tables/. The .xyz files go to generated_files/xyz/<name>/ and the
Gaussian 16 input files (.gjf, optional) to generated_files/gaussian/<name>/.
All the folders are inside the Pólya Pipeline directory.

Usage:
    python derivatives.py --smiles "c1ccccc1" --subst N
    python derivatives.py --smiles "c1ccc2ccccc2c1" --subst N,P
    python derivatives.py --smiles "c1ccccc1" --subst N --max-subst 2 --no-xyz
    python derivatives.py --smiles "c1ccccc1" --name benzene --subst N
    python derivatives.py --smiles "c1ccccc1" --subst N --gaussian
    python derivatives.py --draw "n1ccccc1"

Requires: numpy, rdkit, pillow

Copyright (c) 2026 Jorge Tena Martínez. Released under the MIT License (see LICENSE).
"""

__version__ = "1.0.0"

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb, prod
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# All output goes to the Pólya Pipeline directory, regardless of the working directory
OUTPUT_DIR = Path(__file__).resolve().parent / "generated_files"

# Above this many candidate patterns the enumeration becomes slow
LARGE_ENUMERATION = 2_000_000


def safe_filename(text):
    """Replaces characters that are not allowed in file/folder names."""
    return re.sub(r'[^\w.+-]', '_', text).strip('._') or "molecule"


# ══════════════════════════════════════════════════════════════════════════════
# 1. GEOMETRY AND SYMMETRY
# ══════════════════════════════════════════════════════════════════════════════

def mol_to_atoms(mol):
    """Converts an RDKit mol with a conformer into a list [(symbol, np.array)]."""
    conf = mol.GetConformer()
    atoms = []
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        atoms.append((a.GetSymbol(), np.array([p.x, p.y, p.z])))
    return atoms


def molecular_plane_normal(positions):
    c = positions.mean(axis=0)
    _, _, Vt = np.linalg.svd(positions - c)
    n = Vt[-1]
    return n / np.linalg.norm(n)


def plane_basis(normal):
    n = normal / np.linalg.norm(normal)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(n, ref)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    u = np.cross(n, ref); u /= np.linalg.norm(u)
    v = np.cross(n, u);   v /= np.linalg.norm(v)
    return u, v


def rotation_matrix(axis, angle):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def reflection_matrix(normal):
    n = normal / np.linalg.norm(normal)
    return np.eye(3) - 2 * np.outer(n, n)


def apply_op(positions, matrix, center):
    return (matrix @ (positions - center).T).T + center


def find_permutation(original, transformed, tol):
    n = len(original)
    perm = {}
    for i, tp in enumerate(transformed):
        dists = np.linalg.norm(original - tp, axis=1)
        j = int(np.argmin(dists))
        if dists[j] < tol:
            perm[i] = j
        else:
            return None
    if len(set(perm.values())) != n:
        return None
    return [perm[i] for i in range(n)]


def atoms_preserved(orig_pos, orig_sym, trans_pos, tol):
    for tp, ts in zip(trans_pos, orig_sym):
        dists = np.linalg.norm(orig_pos - tp, axis=1)
        j = int(np.argmin(dists))
        if dists[j] > tol or orig_sym[j] != ts:
            return False
    return True


def perm_preserves_distances(h_pos, perm, dtol=0.035):
    """Rejects spurious permutations by checking the H-H distance matrix."""
    n = len(h_pos)
    for i in range(n):
        for j in range(i + 1, n):
            d_orig = np.linalg.norm(h_pos[i] - h_pos[j])
            d_perm = np.linalg.norm(h_pos[perm[i]] - h_pos[perm[j]])
            if abs(d_orig - d_perm) > dtol:
                return False
    return True


def detect_symmetry_operations(atoms, h_pos, tol=0.15, use_reflections=True):
    """Returns the list of point-group permutations acting on the H positions."""
    n_h = len(h_pos)
    if n_h == 0:
        return []
    all_sym = [sym for sym, _ in atoms]
    all_pos = np.array([pos for _, pos in atoms])
    heavy_pos = np.array([pos for sym, pos in atoms if sym != 'H'])
    center = heavy_pos.mean(axis=0)
    principal = molecular_plane_normal(all_pos)
    operations = []

    for k in range(1, n_h + 1):
        angle_unit = 2 * np.pi / k
        for m in range(k):
            angle = angle_unit * m
            R = rotation_matrix(principal, angle)
            th = apply_op(h_pos, R, center)
            perm = find_permutation(h_pos, th, tol)
            if perm is None:
                continue
            ta = apply_op(all_pos, R, center)
            if atoms_preserved(all_pos, all_sym, ta, tol):
                if perm not in operations and perm_preserves_distances(h_pos, perm):
                    operations.append(perm)

    if not use_reflections:
        return operations

    u, v = plane_basis(principal)
    h_centered = h_pos - center
    h_angles = [np.arctan2(np.dot(p, v), np.dot(p, u)) for p in h_centered]
    trial_angles = list(h_angles)
    sorted_a = sorted(h_angles)
    for i in range(len(sorted_a)):
        mid = (sorted_a[i] + sorted_a[(i + 1) % len(sorted_a)]) / 2
        trial_angles.append(mid)
        trial_angles.append(mid + np.pi)

    for theta in trial_angles:
        d = np.cos(theta) * u + np.sin(theta) * v
        plane_n = np.cross(principal, d)
        if np.linalg.norm(plane_n) < 1e-10:
            continue
        plane_n /= np.linalg.norm(plane_n)
        Ref = reflection_matrix(plane_n)
        th = apply_op(h_pos, Ref, center)
        perm = find_permutation(h_pos, th, tol)
        if perm is None:
            continue
        ta = apply_op(all_pos, Ref, center)
        if atoms_preserved(all_pos, all_sym, ta, tol):
            if perm not in operations and perm_preserves_distances(h_pos, perm):
                operations.append(perm)

    return operations


# ══════════════════════════════════════════════════════════════════════════════
# 2. POLYNOMIALS AND PÓLYA'S THEOREM
# ══════════════════════════════════════════════════════════════════════════════

VARIABLES = "xyzwuvtr"


def poly_mul(p, q):
    r = defaultdict(int)
    for ea, ca in p.items():
        for eb, cb in q.items():
            r[tuple(a + b for a, b in zip(ea, eb))] += ca * cb
    return dict(r)


def poly_pow(p, k):
    n = len(next(iter(p)))
    r = {(0,) * n: 1}
    base = p
    while k:
        if k & 1:
            r = poly_mul(r, base)
        base = poly_mul(base, base)
        k >>= 1
    return r


def poly_add(p, q, factor=1):
    r = defaultdict(int, p)
    for e, c in q.items():
        r[e] += factor * c
    return {e: c for e, c in r.items() if c != 0}


def s_j_poly(j, nvars):
    """s_j = 1 + x^j + y^j + ..."""
    p = {(0,) * nvars: 1}
    for v in range(nvars):
        e = [0] * nvars
        e[v] = j
        p[tuple(e)] = 1
    return p


def monomial_str(e):
    parts = []
    for v, k in zip(VARIABLES, e):
        if k == 1:
            parts.append(v)
        elif k > 1:
            parts.append(f"{v}^{k}")
    return "".join(parts)


def poly_str(p):
    terms = []
    for e in sorted(p, key=lambda e: (sum(e), [-k for k in e])):
        c, m = p[e], monomial_str(e)
        if not m:
            terms.append(str(c))
        elif c == 1:
            terms.append(m)
        else:
            terms.append(f"{c}{m}")
    return " + ".join(terms)


def cycle_structure_xyz(perm):
    """Cycle structure of a permutation (list) → dict {length: count}."""
    n = len(perm)
    visited = [False] * n
    cs = defaultdict(int)
    for start in range(n):
        if not visited[start]:
            length = 0
            j = start
            while not visited[j]:
                visited[j] = True
                j = perm[j]
                length += 1
            cs[length] += 1
    return dict(cs)


def cycle_index_from_ops(operations):
    """Computes the cycle index Z(G) from the list of permutations."""
    ci = defaultdict(int)
    for perm in operations:
        cs = cycle_structure_xyz(perm)
        key = tuple(sorted(cs.items()))
        ci[key] += 1
    return dict(ci)


def term_str(coef, cycles):
    s = " ".join(f"s{j}^{k}" if k > 1 else f"s{j}" for j, k in sorted(cycles.items()))
    return f"{coef} {s}" if coef != 1 else s


def ci_to_terms(ci_dict):
    """Converts the cycle index dict into a list [(coef, cycles_dict)]."""
    terms = []
    for cs_key, coef in ci_dict.items():
        terms.append((coef, dict(cs_key)))
    return terms


def polya_eval(terms, nvars):
    """Expands the Pólya polynomial. Returns (order, n_pos, polynomial)."""
    order = sum(c for c, _ in terms)
    positions = {sum(j * k for j, k in cycles.items()) for _, cycles in terms}
    if len(positions) != 1:
        raise ValueError("The terms do not have the same number of positions.")
    n = positions.pop()
    numerator = {}
    for coef, cycles in terms:
        p = {(0,) * nvars: 1}
        for j, k in cycles.items():
            p = poly_mul(p, poly_pow(s_j_poly(j, nvars), k))
        numerator = poly_add(numerator, p, coef)
    if any(c % order for c in numerator.values()):
        raise ValueError(f"The coefficients are not divisible by |G| = {order}.")
    return order, n, {e: c // order for e, c in numerator.items()}


# ══════════════════════════════════════════════════════════════════════════════
# 3. ENUMERATION OF DERIVATIVES AND SMILES
# ══════════════════════════════════════════════════════════════════════════════

def substitutable_positions(mol):
    """[(label, atom_index)] of the substitutable CH atoms (or those mapped with :n)."""
    mapped = [(a.GetAtomMapNum(), a.GetIdx()) for a in mol.GetAtoms() if a.GetAtomMapNum()]
    if mapped:
        return sorted(mapped)
    ch = [a.GetIdx() for a in mol.GetAtoms()
          if a.GetAtomicNum() == 6 and a.GetTotalNumHs() > 0]
    return [(k + 1, i) for k, i in enumerate(ch)]


def _ring_only_positions(mol):
    """Aromatic CH atoms when the molecule also has other CH (CH3, …), else None."""
    ch = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() > 0]
    ring = [a.GetIdx() for a in ch if a.GetIsAromatic()]
    if not ring or len(ring) == len(ch):
        return None
    return ring


# ── IUPAC locants ─────────────────────────────────────────────────────────────
# Heteroatom seniority for ring numbering (O > S > Se > Te > N > P > As > … > B)
HETERO_PRIORITY = {8: 0, 16: 1, 34: 2, 52: 3, 7: 4, 15: 5, 33: 6, 51: 7, 83: 8,
                   14: 9, 32: 10, 50: 11, 82: 12, 5: 13, 80: 14}

# Fused parents with their IUPAC numbering (atom map = locant; fusion atoms unmapped)
FUSED_TEMPLATES = [
    ("naphthalene",    "[c:1]1[c:2][c:3][c:4]c2[c:5][c:6][c:7][c:8]c12"),
    ("quinoline",      "[n:1]1[c:2][c:3][c:4]c2[c:5][c:6][c:7][c:8]c12"),
    ("isoquinoline",   "[c:1]1[n:2][c:3][c:4]c2[c:5][c:6][c:7][c:8]c12"),
    ("azulene",        "[c:1]1[c:2][c:3]c2[c:4][c:5][c:6][c:7][c:8]c12"),
    ("indole",         "[n:1]1[c:2][c:3]c2[c:4][c:5][c:6][c:7]c12"),
    ("benzofuran",     "[o:1]1[c:2][c:3]c2[c:4][c:5][c:6][c:7]c12"),
    ("benzothiophene", "[s:1]1[c:2][c:3]c2[c:4][c:5][c:6][c:7]c12"),
    ("anthracene",     "[c:1]1[c:2][c:3][c:4]c2[c:10]c3[c:5][c:6][c:7][c:8]c3[c:9]c12"),
    ("phenanthrene",   "[c:1]1[c:2][c:3][c:4]c2c1[c:10][c:9]c1[c:8][c:7][c:6][c:5]c21"),
    ("pyrene",         "[c:1]1[c:2][c:3]c2[c:4][c:5]c3[c:6][c:7][c:8]c4[c:9][c:10]c1c2c34"),
]
_TEMPLATE_MOLS = None


def locant_key(label):
    """Sort key for locants: 2 < 2′ < 3 (number, primes); program numbers sort as numbers."""
    m = re.match(r"(\d+)(′*)", str(label))
    return (int(m.group(1)), len(m.group(2))) if m else (10 ** 6, 0)


def _ring_systems(mol):
    """Atom sets of the ring systems (rings that share a bond are one system)."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    systems = []
    for r in rings:
        merged = [s for s in systems if len(s["atoms"] & r) >= 2]
        new = {"atoms": set(r), "rings": [r]}
        for s in merged:
            new["atoms"] |= s["atoms"]
            new["rings"] += s["rings"]
            systems.remove(s)
        systems.append(new)
    return systems


def _substituted(mol, atoms):
    """Ring atoms of a system that carry a substituent (any heavy atom outside the system)."""
    return {i for i in atoms
            if any(n.GetAtomicNum() > 1 and n.GetIdx() not in atoms
                   for n in mol.GetAtomWithIdx(i).GetNeighbors())}


def _number_monocycle(mol, ring_order):
    """IUPAC numbering of a single ring: senior heteroatom = 1, then lowest locants to
    the heteroatoms, then to the substituents. Returns {atom index: locant}."""
    n = len(ring_order)
    subst = _substituted(mol, set(ring_order))
    prio = {i: HETERO_PRIORITY.get(mol.GetAtomWithIdx(i).GetAtomicNum())
            for i in ring_order if mol.GetAtomWithIdx(i).GetAtomicNum() != 6}
    best = None
    for start in range(n):
        for step in (1, -1):
            order = [ring_order[(start + step * k) % n] for k in range(n)]
            het = [(prio[i] if prio[i] is not None else 99, k + 1)
                   for k, i in enumerate(order) if i in prio]
            if het and (prio[order[0]] if order[0] in prio else 100) != min(p for p, _ in het):
                continue
            key = (sorted(l for _, l in het), [l for _, l in sorted(het)],
                   sorted(k + 1 for k, i in enumerate(order) if i in subst))
            if best is None or key < best[0]:
                best = (key, order)
    return {i: k + 1 for k, i in enumerate(best[1])}


def _number_fused(mol, atoms):
    """IUPAC numbering of a fused system from FUSED_TEMPLATES, or None if unknown."""
    global _TEMPLATE_MOLS
    if _TEMPLATE_MOLS is None:
        _TEMPLATE_MOLS = [Chem.MolFromSmarts(s) for _, s in FUSED_TEMPLATES]
    subst = _substituted(mol, atoms)
    for tmpl in _TEMPLATE_MOLS:
        if tmpl.GetNumAtoms() != len(atoms):
            continue
        best = None
        for match in mol.GetSubstructMatches(tmpl, uniquify=False, useChirality=False):
            if set(match) != atoms:
                continue
            loc = {match[t.GetIdx()]: t.GetAtomMapNum() for t in tmpl.GetAtoms()
                   if t.GetAtomMapNum()}
            key = sorted(loc[i] for i in subst if i in loc)
            if best is None or key < best[0]:
                best = (key, loc)
        if best:
            return best[1]
    return None


def iupac_locants(mol, atoms):
    """IUPAC locants ('2', '1′'…) for the given ring atoms, or None when any of them
    is outside a ring system the program can number (then program numbering is used)."""
    systems = [s for s in _ring_systems(mol) if s["atoms"] & set(atoms)]
    if not systems or any(not any(i in s["atoms"] for s in systems) for i in atoms):
        return None
    systems.sort(key=lambda s: min(s["atoms"]))
    labels = {}
    for primes, s in enumerate(systems):
        if len(s["rings"]) == 1:
            ring = next(r for r in mol.GetRingInfo().AtomRings() if set(r) == s["atoms"])
            loc = _number_monocycle(mol, list(ring))
        else:
            loc = _number_fused(mol, s["atoms"])
            if loc is None:
                return None
        for i, l in loc.items():
            labels[i] = f"{l}{'′' * primes}"
    if any(i not in labels for i in atoms):
        return None
    return {i: labels[i] for i in atoms}


NUMBERING_NAMES = {"iupac": "IUPAC numbering", "user": "numbered in the SMILES",
                   "program": "program numbering, not IUPAC"}


def substitution_positions(mol):
    """Substitution positions and how they are numbered.

    Returns ([(label, atom index)] sorted by label, scheme) with scheme:
      'user'    → atom maps written in the SMILES;
      'iupac'   → IUPAC locants of the parent ring system(s);
      'program' → 1, 2, 3… (ring system the program cannot number).
    Only aromatic CH are positions when the molecule also has other CH (CH3, …).
    """
    if any(a.GetAtomMapNum() for a in mol.GetAtoms()):
        return substitutable_positions(mol), "user"
    atoms = _ring_only_positions(mol) or [i for _, i in substitutable_positions(mol)]
    loc = iupac_locants(mol, atoms) if atoms else None
    if loc:
        return sorted(((loc[i], i) for i in atoms), key=lambda p: locant_key(p[0])), "iupac"
    return [(k + 1, i) for k, i in enumerate(atoms)], "program"


def number_ring_positions(mol):
    """Numbers the aromatic CH atoms with atom maps when there are other CH (CH3, …).

    The map numbers are the IUPAC locants when they are plain numbers, otherwise
    1, 2, 3… Returns (numbered mol, numbered SMILES), or None when nothing has to
    change: the SMILES already has atom maps, every CH is aromatic, or none is.
    """
    if any(a.GetAtomMapNum() for a in mol.GetAtoms()) or not _ring_only_positions(mol):
        return None
    positions, scheme = substitution_positions(mol)
    plain = all(str(l).isdigit() for l, _ in positions)
    numbered = Chem.Mol(mol)
    for k, (label, i) in enumerate(positions, 1):
        numbered.GetAtomWithIdx(i).SetAtomMapNum(int(label) if plain else k)
    return numbered, Chem.MolToSmiles(numbered, canonical=False)


def lowest_locant_pattern(pattern, group, labels):
    """The pattern of the same orbit whose substituents get the lowest locants (IUPAC)."""
    n = len(pattern)
    keys = [locant_key(l) for l in labels]

    def key(b):
        return (sorted(keys[i] for i, v in enumerate(b) if v),
                [sorted(keys[i] for i, v in enumerate(b) if v == s)
                 for s in range(1, max(b, default=0) + 1)])
    return min({tuple(pattern[g[i]] for i in range(n)) for g in group}, key=key)


def permutation_group_rdkit(mol, idx_pos, respect_charges=False):
    """Permutation group of the positions from graph automorphisms (RDKit)."""
    query = Chem.Mol(mol)
    for a in query.GetAtoms():
        a.SetAtomMapNum(0)
        if not respect_charges:
            a.SetFormalCharge(0)
    hs = [a.GetTotalNumHs() for a in mol.GetAtoms()]
    where = {at: k for k, at in enumerate(idx_pos)}
    group = set()
    for m in query.GetSubstructMatches(query, uniquify=False, maxMatches=1_000_000):
        if any(hs[i] != hs[m[i]] for i in range(len(m))):
            continue
        if all(m[at] in where for at in idx_pos):
            group.add(tuple(where[m[at]] for at in idx_pos))
    return sorted(group)


def cycle_type_rdkit(g):
    seen, ctype = set(), Counter()
    for i in range(len(g)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            j = g[j]
            length += 1
        ctype[length] += 1
    return tuple(sorted(ctype.items()))


def cycle_index_rdkit(group):
    return Counter(cycle_type_rdkit(g) for g in group)


def index_str_rdkit(index):
    terms = sorted(index.items(), key=lambda t: (
        -t[0][0][1] if t[0][0][0] == 1 else 0, t[0]))
    return " + ".join(term_str(c, dict(t)) for t, c in terms)


def count_patterns(n, nsubst, max_subst=None):
    """Number of candidate patterns that enumerate_derivatives will check."""
    if max_subst is None or max_subst >= n:
        return (nsubst + 1) ** n
    return sum(comb(n, k) * nsubst ** k for k in range(max_subst + 1))


def candidate_patterns(n, nsubst, max_subst=None):
    """Yields every pattern with at most max_subst substituted positions."""
    if max_subst is None or max_subst >= n:
        yield from product(range(nsubst + 1), repeat=n)
        return
    for k in range(max_subst + 1):
        for idx in combinations(range(n), k):
            for vals in product(range(1, nsubst + 1), repeat=k):
                a = [0] * n
                for i, v in zip(idx, vals):
                    a[i] = v
                yield tuple(a)


def enumerate_derivatives(group, n, nsubst, max_subst=None):
    """One representative (lexicographic minimum) per orbit + the orbit size."""
    derivatives = []
    for a in candidate_patterns(n, nsubst, max_subst):
        images = set()
        is_minimum = True
        for g in group:
            b = tuple(a[g[i]] for i in range(n))
            if b < a:
                is_minimum = False
                break
            images.add(b)
        if is_minimum:
            derivatives.append((a, len(images)))
    derivatives.sort(key=lambda d: (
        sum(1 for v in d[0] if v),
        [-Counter(d[0])[s] for s in range(1, nsubst + 1)],
        d[0]))
    return derivatives


def read_heteroatom(text):
    m = Chem.MolFromSmiles(text)
    if m is None or m.GetNumAtoms() != 1:
        return None
    a = m.GetAtomWithIdx(0)
    return a.GetAtomicNum(), a.GetFormalCharge()


def derivative_mol(mol, idx_pos, pattern, heteroatoms):
    """The derivative as an RDKit mol that keeps the parent's atom indices (None if impossible).

    Replacing a CH must not add H atoms. When a substituted atom carries the ring's
    charge in the SMILES (e.g. the C+ of tropylium) and the direct replacement gains an H
    or is impossible, the charge is moved to a free ring carbon first: for an ion whose
    charge is delocalised over the ring this is the same species.
    """
    direct, h_ok = _substitute(mol, idx_pos, pattern, heteroatoms, move_charges=False)
    if direct is not None and h_ok:
        return direct
    moved, moved_ok = _substitute(mol, idx_pos, pattern, heteroatoms, move_charges=True)
    if moved is not None and (moved_ok or direct is None):
        return moved
    return direct


def build_smiles(mol, idx_pos, pattern, heteroatoms):
    """Builds the SMILES of the derivative by changing the atoms at the positions."""
    m = derivative_mol(mol, idx_pos, pattern, heteroatoms)
    return None if m is None else Chem.MolToSmiles(m)


def _substitute(mol, idx_pos, pattern, heteroatoms, move_charges):
    """(derivative mol or None, True if no substituted atom gained H atoms)."""
    rw = Chem.RWMol(mol)
    substituted = [idx_pos[k] for k, s in enumerate(pattern) if s]
    if move_charges:
        free = [a.GetIdx() for a in rw.GetAtoms()
                if a.GetAtomicNum() == 6 and a.GetIsAromatic() and a.GetFormalCharge() == 0
                and a.GetTotalNumHs() > 0 and a.GetIdx() not in substituted]
        for i in substituted:
            a = rw.GetAtomWithIdx(i)
            if a.GetFormalCharge() and free:
                b = rw.GetAtomWithIdx(free.pop(0))
                b.SetNumExplicitHs(b.GetTotalNumHs())
                b.SetNoImplicit(True)
                b.SetFormalCharge(a.GetFormalCharge())
                a.SetFormalCharge(0)
    for k, s in enumerate(pattern):
        if not s:
            continue
        z, charge = heteroatoms[s - 1]
        at = rw.GetAtomWithIdx(idx_pos[k])
        at.SetAtomicNum(z)
        at.SetFormalCharge(at.GetFormalCharge() + charge)
        at.SetNumExplicitHs(0)
        at.SetNoImplicit(False)
    for a in rw.GetAtoms():
        a.SetAtomMapNum(0)
    try:
        Chem.SanitizeMol(rw)
    except Exception:
        return None, False
    h_ok = all(rw.GetAtomWithIdx(i).GetTotalNumHs()
               <= mol.GetAtomWithIdx(i).GetTotalNumHs() - 1 for i in substituted)
    return rw.GetMol(), h_ok


# ══════════════════════════════════════════════════════════════════════════════
# 4. XYZ GEOMETRY GENERATION
# ══════════════════════════════════════════════════════════════════════════════

def geometry_3d(smiles, seed=0xF00D):
    """Generates the 3D geometry of a SMILES with ETKDG + MMFF/UFF optimization."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError("could not generate the 3D geometry")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        force_field = "MMFF"
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        force_field = "UFF"
    return mol, force_field


def write_xyz(path, mol):
    """Writes an .xyz file with the header 'charge = q; mult = m'."""
    charge = Chem.GetFormalCharge(mol)
    mult = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) + 1
    conf = mol.GetConformer()
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"{mol.GetNumAtoms()}\n")
        f.write(f"charge = {charge}; mult = {mult}\n")
        for a in mol.GetAtoms():
            p = conf.GetAtomPosition(a.GetIdx())
            f.write(f"{a.GetSymbol():<2} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}\n")


DEFAULT_ROUTE = "#b3lyp/6-311++G(d,p) opt freq"


def write_gjf(path, mol, title, route=DEFAULT_ROUTE):
    """Writes a Gaussian 16 input file: route, title, 'charge mult' and coordinates."""
    charge = Chem.GetFormalCharge(mol)
    mult = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) + 1
    conf = mol.GetConformer()
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"{route}\n\n{title}\n\n{charge} {mult}\n")
        for a in mol.GetAtoms():
            p = conf.GetAtomPosition(a.GetIdx())
            f.write(f"{a.GetSymbol():<2} {p.x:17.5f}{p.y:10.5f}{p.z:10.5f}\n")
        f.write("\n\n")          # Gaussian needs a blank line after the molecule


# ══════════════════════════════════════════════════════════════════════════════
# 5. SMILES DRAWING
# ══════════════════════════════════════════════════════════════════════════════

def smiles_problem(smiles):
    """Plain-language reason why a SMILES cannot be used, or None when it is fine."""
    if not smiles.strip():
        return "Enter the SMILES of the parent molecule."
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return "Cannot read this SMILES: check brackets, ring-closure digits and atom symbols."
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        msg = str(e)
        m = re.search(r"atom # ?(\d+) (\w+)", msg)
        if "valence" in msg.lower() and m:
            return (f"Impossible valence: atom {int(m.group(1)) + 1} ({m.group(2)}) "
                    "has too many bonds.")
        if "kekulize" in msg.lower():
            return ("This aromatic ring cannot be drawn with alternating bonds: check the "
                    "lowercase atoms (pyrrole N is [nH], cyclopentadienyl is [cH-]).")
        return f"RDKit cannot use this SMILES: {msg.strip()}"
    if not substitution_positions(mol)[0]:
        return "This molecule has no CH positions to substitute."
    return None


def smiles_image(smiles, size=(500, 500), background=None, show_positions=True,
                 label_colour=(0.11, 0.31, 0.57)):
    """2D structure of a SMILES as a PIL image (rendered in memory with Cairo).

    show_positions: labels the substitution positions with the same numbers the
    results table uses (ring numbering applied automatically, as in run_pipeline).
    background: RGB tuple (0-1) for the background (default white).
    Raises ValueError for an invalid SMILES.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    labels = {i: l for l, i in substitution_positions(mol)[0]} if show_positions else {}
    return mol_image(mol, size, labels, background=background, label_colour=label_colour)


def mol_image(mol, size, labels=None, highlight=(), background=None,
              label_colour=(0.11, 0.31, 0.57), highlight_colour=(0.80, 0.87, 0.96)):
    """Draws a mol (PIL image). labels: {atom index: position label}; highlight: atoms
    drawn on a coloured disc. The 2D layout is the deterministic one of _layout_2d."""
    import io
    from PIL import Image
    from rdkit.Chem.Draw import rdMolDraw2D
    mol = Chem.Mol(mol)
    if not mol.GetNumConformers():
        _layout_2d(mol)
    for i, label in (labels or {}).items():
        mol.GetAtomWithIdx(i).SetProp("atomNote", str(label))
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)            # drawn as a note, not as "C:1"
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = d.drawOptions()
    if background is not None:
        opts.setBackgroundColour((*background, 1.0))
    opts.annotationFontScale = 0.75
    opts.setAtomNoteColour((*label_colour, 1.0))
    opts.padding = 0.08
    hl = list(highlight)
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hl,
                                       highlightAtomColors={i: highlight_colour for i in hl},
                                       highlightBonds=[])
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))


def _layout_2d(mol):
    from rdkit.Chem import rdDepictor
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)


def derivative_image(parent_smiles, positions, substituents, size, **kw):
    """Drawing of one derivative in the parent's layout, its positions numbered as the
    parent's and the replaced atoms highlighted.

    positions: the "Positions" text of the results ("N: 2,4; P: 3" or "—").
    Returns (image, SMILES). Raises ValueError when the derivative has an impossible valence.
    """
    parent = Chem.MolFromSmiles(parent_smiles)
    _layout_2d(parent)
    pos = substitution_positions(parent)[0]
    by_label = {str(l): k for k, (l, _) in enumerate(pos)}
    idx_pos = [i for _, i in pos]
    pattern = [0] * len(pos)
    for part in positions.split(";"):
        if ":" not in part:
            continue
        atom, labels = part.split(":", 1)
        s = substituents.index(atom.strip()) + 1
        for l in labels.split(","):
            pattern[by_label[l.strip()]] = s
    het = [read_heteroatom(s) for s in substituents]
    m = derivative_mol(parent, idx_pos, pattern, het)
    if m is None:
        raise ValueError("This derivative would have an impossible valence.")
    m.RemoveAllConformers()
    m.AddConformer(parent.GetConformer(), assignId=True)
    image = mol_image(m, size, {i: l for l, i in pos},
                      highlight=[idx_pos[k] for k, v in enumerate(pattern) if v], **kw)
    return image, Chem.MolToSmiles(m)


def draw_smiles(smiles):
    """Shows the 2D structure of a SMILES in a window."""
    try:
        smiles_image(smiles).show()
    except ValueError as e:
        print(f"  {e}")


# ══════════════════════════════════════════════════════════════════════════════
# 6. MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def run_pipeline(parent_smiles, subst_text, name=None, max_subst=None,
                 tol=0.15, generate_xyz_files=True,
                 respect_charges=False, generate_gaussian_files=False,
                 gaussian_route=DEFAULT_ROUTE, info=None):
    """
    info: optional dict that is filled with the run's summary (positions, symmetry
    check, Pólya total, distinct compounds, files written and their folders).

    Runs the complete pipeline from a SMILES:
      1. Picks the substitution positions (only the aromatic CH when there are other CH).
      2. Generates the 3D geometry (ETKDG + MMFF/UFF) → detects symmetry → geometric Z(G).
      3. Computes Z(G) from graph automorphisms (RDKit).
      4. Cross-checks both indices.
      5. Counts the derivatives with Pólya's theorem.
      6. Enumerates one pattern per orbit, builds its SMILES and removes duplicate compounds.
      7. Generates the .xyz and/or Gaussian 16 (.gjf) files of the derivatives.
      8. Prints the results table to the screen.
      9. Saves the results table as CSV.
    """
    # ── 0. Parse the SMILES ───────────────────────────────────────────────────
    mol = Chem.MolFromSmiles(parent_smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {parent_smiles}")
    if generate_gaussian_files and not gaussian_route.strip().startswith("#"):
        raise ValueError(f"The Gaussian route must start with '#': {gaussian_route}")

    if name is None:
        from rdkit.Chem import rdMolDescriptors
        name = rdMolDescriptors.CalcMolFormula(mol)
    info = {} if info is None else info
    info.update(name=name, numbered_smiles=None)

    print(f"\n{'═'*65}")
    print(f"  Molecule: {name}   SMILES: {parent_smiles}")
    print(f"{'═'*65}")

    # ── 1. Substitution positions ─────────────────────────────────────────────
    # CH outside the aromatic rings (CH3, …) are left out: the ring CH are numbered
    numbered = number_ring_positions(mol)
    if numbered:
        numbered_smiles = numbered[1]
        info["numbered_smiles"] = numbered_smiles
        print("\n  ℹ Only the aromatic CH are substituted (CH outside the rings, e.g. CH3,")
        print("    are left as they are). The ring positions were numbered automatically:")
        print(f"      {numbered_smiles}")
        print("    To substitute other CH too, number them yourself in the SMILES.")

    pos, numbering = substitution_positions(mol)
    idx_pos = [i for _, i in pos]
    labels = [lab for lab, _ in pos]
    n_rdkit = len(pos)
    if n_rdkit == 0:
        raise ValueError("There are no substitutable CH positions in this molecule.")
    info["numbering"] = numbering
    print(f"\n  Positions ({NUMBERING_NAMES[numbering]}): " + ", ".join(str(l) for l in labels))

    # ── 2. Symmetry from the 3D geometry (generated with ETKDG) ──────────────
    print("\n  Generating 3D geometry for symmetry detection…", end=" ", flush=True)
    parent_mol3d, parent_ff = geometry_3d(parent_smiles)
    print(f"OK ({parent_ff})")

    # Only the H atoms of the positions are compared (in position order), so H atoms
    # of CH3, OH, NH… do not count; every heavy atom must still be preserved
    atoms = mol_to_atoms(parent_mol3d)
    pos_h = [[n.GetIdx() for n in parent_mol3d.GetAtomWithIdx(i).GetNeighbors()
              if n.GetAtomicNum() == 1] for i in idx_pos]
    check_3d = all(len(h) == 1 for h in pos_h)
    if check_3d:
        heavy = [a for a in atoms if a[0] != 'H']
        h_pos = np.array([atoms[h[0]][1] for h in pos_h])
        ops = detect_symmetry_operations(heavy + [atoms[h[0]] for h in pos_h],
                                         h_pos, tol=tol)
        if not ops:
            raise ValueError("No symmetry operations were detected. "
                             "Try --tol 0.20.")
        ci_xyz = cycle_index_from_ops(ops)
        terms_xyz = ci_to_terms(ci_xyz)
        order_xyz = len(ops)
        zi_str = " + ".join(term_str(c, s) for c, s in
                            sorted(terms_xyz, key=lambda x: list(x[1].items())))

        print(f"\n  [3D]     n = {len(h_pos)} H positions  |  |G| = {order_xyz}")
        print(f"           Z(G) = ({zi_str}) / {order_xyz}")
    else:
        print("\n  [3D]     skipped (some positions carry more than one H)")

    # ── 3. Symmetry from the graph (RDKit) ────────────────────────────────────
    group = permutation_group_rdkit(mol, idx_pos, respect_charges)
    index_rdk = cycle_index_rdkit(group)
    order_rdk = len(group)

    print(f"\n  [Graph]  n = {n_rdkit} positions  |  |G| = {order_rdk}")
    print(f"           Z(G) = ({index_str_rdkit(index_rdk)}) / {order_rdk}")

    # ── 4. Cross-validation ───────────────────────────────────────────────────
    # The graph group is always used. A less symmetric 3D geometry (twisted rings,
    # force-field distortion) does not change the results, so it is only a note
    info.update(positions=n_rdkit, group_order=order_rdk, check_3d="skipped")
    if check_3d:
        norm_xyz = {t: Fraction(c, order_xyz) for t, c in ci_xyz.items()}
        norm_rdk = {tuple(sorted(t)): Fraction(c, order_rdk)
                    for t, c in index_rdk.items()}
        info["group_order_3d"] = order_xyz
        if norm_xyz == norm_rdk:
            info["check_3d"] = "match"
            print("\n  ✓ Cycle indices match (3D geometry = graph)")
        elif order_xyz < order_rdk and order_rdk % order_xyz == 0:
            info["check_3d"] = "less"
            print(f"\n  ℹ The 3D geometry is less symmetric than the graph (|G| = {order_xyz}"
                  f" vs {order_rdk}):")
            print("    twisted/flexible parts or a force-field distortion. The graph symmetry")
            print("    is used, so the results are not affected.")
        else:
            info["check_3d"] = "mismatch"
            print("\n  ⚠ The indices do NOT match — the graph one is used for enumeration.")
            print("    The 3D symmetry is not consistent with the graph: try another --tol.")

    # Use the graph group for enumeration (it is the canonical one for SMILES)
    terms_eval = [(c, dict(t)) for t, c in index_rdk.items()]

    # ── 5. Pólya counting ─────────────────────────────────────────────────────
    subst_list = [s.strip() for s in subst_text.split(",") if s.strip()]
    nvars = len(subst_list)
    heteroatoms = [read_heteroatom(s) for s in subst_list]
    if not all(heteroatoms):
        raise ValueError(f"Unrecognized substituent in: {subst_text}")

    order_p, n_p, P = polya_eval(terms_eval, nvars)
    by_degree = defaultdict(int)
    for e, c in P.items():
        by_degree[sum(e)] += c

    print(f"\n{'─'*65}")
    print(f"  Substituents: {', '.join(subst_list)}")
    print(f"  P({', '.join(VARIABLES[:nvars])}) = {poly_str(P)}")
    print(f"\n  {'# subst.':>10} | {'Pólya count':>20}")
    print(f"  {'─'*10}-+-{'─'*20}")
    for g in sorted(by_degree):
        print(f"  {g:>10} | {by_degree[g]:>20}")
    print(f"  {'─'*10}-+-{'─'*20}")
    total = sum(P.values())
    info.update(substituents=subst_list, polya_total=total, max_subst=max_subst)
    print(f"  {'Total':>10} | {total:>20}")
    print(f"{'─'*65}")

    # ── 6. Enumeration and SMILES construction ────────────────────────────────
    n_candidates = count_patterns(n_rdkit, nvars, max_subst)
    if n_candidates > LARGE_ENUMERATION:
        print(f"\n  ⚠ {n_candidates:,} patterns to check — this may take a long time."
              "\n    Use --max-subst to limit the number of substituents.")
    print(f"\n  Enumerating derivatives…", end=" ", flush=True)
    derivatives = enumerate_derivatives(group, n_rdkit, nvars, max_subst)
    # Each derivative is shown with the pattern that gives its substituents the lowest locants
    derivatives = [(lowest_locant_pattern(p, group, labels), orbit) for p, orbit in derivatives]
    derivatives.sort(key=lambda d: (
        sum(1 for v in d[0] if v),
        [-Counter(d[0])[s] for s in range(1, nvars + 1)],
        sorted(locant_key(labels[i]) for i, v in enumerate(d[0]) if v)))
    print(f"{len(derivatives)} found.")

    # A substituent equal to a ring heteroatom (N on pyridine…) makes different
    # patterns give the same compound: only the first row is kept and the
    # patterns of the removed ones are added to its "Patterns" column
    rows = []
    by_smiles = {}
    invalid = 0
    for pattern, orbit in derivatives:
        smi = build_smiles(mol, idx_pos, pattern, heteroatoms)
        if smi is None:
            invalid += 1
            smi = "invalid structure (valence)"
        elif smi in by_smiles:
            by_smiles[smi][-2] += orbit
            continue
        counts = [sum(1 for v in pattern if v == s) for s in range(1, nvars + 1)]
        positions_str = "; ".join(
            f"{subst_list[s - 1]}: " + ",".join(
                str(labels[i]) for i, v in enumerate(pattern) if v == s)
            for s in range(1, nvars + 1) if s in pattern) or "—"
        row = [len(rows) + 1, sum(counts), *counts, positions_str, orbit, smi]
        rows.append(row)
        if "invalid" not in smi:
            by_smiles[smi] = row

    n_dup = len(derivatives) - len(rows)
    info.update(patterns=len(derivatives), derivatives=len(rows),
                distinct=len(rows) - invalid, invalid=invalid, duplicates=n_dup,
                files_written=0, file_failures=0, xyz_dir=None, gjf_dir=None)
    if invalid:
        print(f"  ⚠ {invalid} structures with incompatible valence.")
    if n_dup:
        print(f"  ℹ {n_dup} patterns give a compound that was already listed (a substituent")
        print("    is the same element as a heteroatom of the ring). Duplicates removed:")
        patterns_deg = Counter(sum(1 for v in p if v) for p, _ in derivatives)
        compounds_deg = Counter(row[1] for row in rows)
        print(f"\n  {'# subst.':>10} | {'Pólya patterns':>14} | {'Distinct compounds':>18}")
        print(f"  {'─'*10}-+-{'─'*14}-+-{'─'*18}")
        for g in sorted(patterns_deg):
            print(f"  {g:>10} | {patterns_deg[g]:>14} | {compounds_deg[g]:>18}")
        print(f"  {'─'*10}-+-{'─'*14}-+-{'─'*18}")
        print(f"  {'Total':>10} | {len(derivatives):>14} | {len(rows):>18}")
    else:
        print("  ✓ All derivatives are different compounds.")

    # ── 7. XYZ / Gaussian file generation ─────────────────────────────────────
    # The substituent tag keeps runs with different substituents from overwriting each other.
    # Each geometry is optimised once and written to every format requested.
    folder = safe_filename(name)
    tag = safe_filename("-".join(subst_list))
    xyz_dir = OUTPUT_DIR / "xyz" / folder
    gjf_dir = OUTPUT_DIR / "gaussian" / folder
    if generate_xyz_files or generate_gaussian_files:
        print("\n  Generating 3D geometries…")
        for d, ext, wanted in ((xyz_dir, ".xyz", generate_xyz_files),
                               (gjf_dir, ".gjf", generate_gaussian_files)):
            if wanted:
                d.mkdir(parents=True, exist_ok=True)
                print(f"    {ext} → {d}")
        if generate_gaussian_files:
            print(f"    Gaussian route: {gaussian_route.strip()}")
        done, uff_count, failures = 0, 0, []
        for row in rows:
            num = row[0]
            smi = row[-1]
            if "invalid" in smi:
                continue
            stem = f"{folder}_{tag}_{num:03d}"
            try:
                mol3d, force_field = geometry_3d(smi)
                if generate_xyz_files:
                    write_xyz(str(xyz_dir / f"{stem}.xyz"), mol3d)
                if generate_gaussian_files:
                    write_gjf(str(gjf_dir / f"{stem}.gjf"), mol3d,
                              f"{stem}  {smi}", gaussian_route.strip())
                done += 1
                uff_count += force_field == "UFF"
            except Exception as e:
                failures.append(f"{num}: {e}")
        kinds = " + ".join(k for k, wanted in ((".xyz", generate_xyz_files),
                                               (".gjf", generate_gaussian_files)) if wanted)
        info.update(files_written=done, file_failures=len(failures),
                    xyz_dir=str(xyz_dir) if generate_xyz_files else None,
                    gjf_dir=str(gjf_dir) if generate_gaussian_files else None)
        msg = f"  ✓ {done} geometries written ({kinds})"
        if uff_count:
            msg += f"  ({uff_count} optimized with UFF)"
        print(msg)
        for f in failures:
            print(f"    ⚠ FAILED {f}")

    # ── 8. Results table ──────────────────────────────────────────────────────
    # Dynamic column widths
    w_smi = max((len(f[-1]) for f in rows), default=6)
    w_pos = max((len(f[-3]) for f in rows), default=10)
    w_smi = min(w_smi, 80)
    w_pos = min(w_pos, 40)

    # Header
    subst_cols = "".join(f"  {'# '+s:>6}" for s in subst_list)
    print(f"\n  {'#':>4}  {'Subs':>4}{subst_cols}  {'Patterns':>8}  "
          f"{'Positions':<{w_pos}}  {'SMILES'}")
    print(f"  {'─'*4}  {'─'*4}{('  ' + '─'*6) * nvars}  {'─'*8}  "
          f"{'─'*w_pos}  {'─'*w_smi}")

    for row in rows:
        num, n_tot, *counts, pos_str, orbit, smi = row
        counts_str = "".join(f"  {c:>6}" for c in counts)
        smi_trunc = smi[:w_smi] if len(smi) <= w_smi else smi[:w_smi - 1] + "…"
        pos_trunc = pos_str[:w_pos] if len(pos_str) <= w_pos else pos_str[:w_pos - 1] + "…"
        print(f"  {num:>4}  {n_tot:>4}{counts_str}  {orbit:>8}  "
              f"{pos_trunc:<{w_pos}}  {smi_trunc}")

    print(f"  {'─'*4}  {'─'*4}{('  ' + '─'*6) * nvars}  {'─'*8}  "
          f"{'─'*w_pos}  {'─'*w_smi}")

    # ── 9. Save the table as CSV (opens directly in Excel) ────────────────────
    csv_dir = OUTPUT_DIR / "tables"
    csv_dir.mkdir(parents=True, exist_ok=True)
    suffix = f"_max{max_subst}" if max_subst is not None else ""
    csv_path = csv_dir / f"{folder}_{tag}{suffix}.csv"
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["#", "Subs", *[f"# {s}" for s in subst_list],
                    "Positions", "Patterns", "SMILES"])
        w.writerows(rows)
    info["csv_path"] = str(csv_path)
    print(f"\n  ✓ Table saved to {csv_path}")

    print(f"{'═'*65}\n")
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# 7. COMMAND-LINE INTERFACE
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description=(
            f"Pólya Pipeline {__version__}: SMILES → symmetry → cycle index "
            "→ derivatives (SMILES + XYZ + Gaussian)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python derivatives.py --smiles 'c1ccccc1' --subst N\n"
            "  python derivatives.py --smiles 'c1ccc2ccccc2c1' --subst N,P\n"
            "  python derivatives.py --smiles 'c1ccccc1' --subst N --max-subst 2 --no-xyz\n"
            "  python derivatives.py --smiles 'c1ccccc1' --name benzene --subst N\n"
            "  python derivatives.py --smiles 'c1ccccc1' --subst N --gaussian\n"
            "  python derivatives.py --smiles 'c1ccccc1' --subst N --gaussian "
            "--route '#m062x/def2tzvp opt freq'\n"
            "  python derivatives.py --draw 'n1ccccc1'\n"
        ),
    )
    ap.add_argument("--smiles", "-s", default=None,
                    help="SMILES of the parent molecule (required)")
    ap.add_argument("--name", "-n", default=None,
                    help="Name of the molecule (default: molecular formula)")
    ap.add_argument("--subst", "-u", default="N",
                    help="Comma-separated substituents [N]")
    ap.add_argument("--max-subst", type=int, default=None,
                    help="Maximum number of substituents per derivative")
    ap.add_argument("--tol", "-t", type=float, default=0.15, metavar="Å",
                    help="Tolerance for symmetry detection [0.15]")
    ap.add_argument("--no-xyz", action="store_true",
                    help="Do not generate .xyz files for the derivatives")
    ap.add_argument("--gaussian", action="store_true",
                    help="Also generate Gaussian 16 input files (.gjf)")
    ap.add_argument("--route", default=DEFAULT_ROUTE,
                    help=f"Route line of the Gaussian files [{DEFAULT_ROUTE}]")
    ap.add_argument("--respect-charges", action="store_true",
                    help="Do not ignore formal charges in the automorphisms")
    ap.add_argument("--draw", metavar="SMILES",
                    help="Draw this SMILES and exit")
    ap.add_argument("--version", action="version",
                    version=f"Pólya Pipeline {__version__}")
    args = ap.parse_args()

    if args.draw:
        draw_smiles(args.draw)
        return

    smiles = args.smiles
    if smiles is None:
        try:
            smiles = input("SMILES of the parent molecule: ").strip()
        except EOFError:
            ap.error("Provide the SMILES with --smiles.")

    try:
        run_pipeline(
            parent_smiles=smiles,
            subst_text=args.subst,
            name=args.name,
            max_subst=args.max_subst,
            tol=args.tol,
            generate_xyz_files=not args.no_xyz,
            respect_charges=args.respect_charges,
            generate_gaussian_files=args.gaussian,
            gaussian_route=args.route,
        )
    except (ValueError, OSError) as e:
        print(f"\n  Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
