#!/usr/bin/env python3
"""
derivatives_gui.py — Pólya Pipeline window
==========================================
Graphical interface for the Pólya Pipeline engine (derivatives.py).

Usage:
    python derivatives_gui.py          # opens the window
    (or double-click PolyaPipeline.bat on Windows)

Keys: Enter runs, Esc cancels a run (or goes back to the parent molecule),
F1 opens the user guide, Ctrl+C copies the SMILES of the selected derivatives.

Requires: tkinter (included with Python), numpy, rdkit, pillow
Must be in the same folder as derivatives.py

Copyright (c) 2026 Jorge Tena Martínez. Released under the MIT License (see LICENSE).
"""

import ctypes
import multiprocessing as mp
import os
import queue
import subprocess
import sys
import time
import tkinter as tk
from pathlib import Path
from tkinter import font as tkfont
from tkinter import messagebox, ttk

from PIL import ImageTk

# Sharp text on scaled displays (125 %, 150 %…): must happen before Tk starts
if sys.platform == "win32":
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
        # Own taskbar entry (and icon) instead of the generic python.exe one
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("PolyaPipeline.GUI")
    except Exception:
        pass


# ── Import the pipeline ───────────────────────────────────────────────────────
try:
    import derivatives as der
except ImportError:
    root = tk.Tk()
    root.withdraw()
    messagebox.showerror(
        "Import error",
        "derivatives.py could not be found.\n"
        "Make sure this file is in the same folder.",
    )
    sys.exit(1)

USER_GUIDE = Path(__file__).resolve().parent / "userguide.txt"


# ── Palette: paper and ink ────────────────────────────────────────────────────
# Text colours are ≥ 4.5:1 on the background they are used on.
PAPER    = "#ffffff"   # content surfaces: molecule, table, log, fields
PANEL    = "#f1f3f6"   # form column and table headings
INK      = "#1a1f27"   # main text
INK2     = "#454e5c"   # labels, secondary text
INK3     = "#5b6472"   # hints
RULE     = "#c5ccd6"   # borders and separators
DISABLED = "#e3e7ec"   # background of disabled fields
ACCENT   = "#1c4e8c"   # primary action, focus, selection, position labels
ACCENT_H = "#163f73"   # primary action pressed / hover
SELECT   = "#dbe5f2"   # selected rows
OK       = "#1d6b35"
WARN     = "#a3261f"
STALE    = "#8a5a00"   # results that no longer match the inputs

UI   = "Segoe UI"
MONO = "Consolas"


def hex_to_rgb(color):
    """'#rrggbb' → (r, g, b) in 0-1, as RDKit expects."""
    return tuple(int(color[i:i + 2], 16) / 255 for i in (1, 3, 5))


def open_path(path):
    """Opens a file or folder with the system's default program."""
    if sys.platform == "win32":
        os.startfile(path)
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


# ── The pipeline runs in its own process, so a run can be cancelled ──────────
class _QueueWriter:
    """sys.stdout of the child process: every print goes to the GUI."""
    def __init__(self, q):
        self._q = q

    def write(self, text):
        if text:
            self._q.put(("text", text))

    def flush(self):
        pass


def pipeline_process(q, kwargs, output_dir):
    """Child process: runs der.run_pipeline and reports through the queue."""
    der.OUTPUT_DIR = Path(output_dir)
    sys.stdout = _QueueWriter(q)
    info = {}
    try:
        rows = der.run_pipeline(**kwargs, info=info)
        q.put(("done", rows, info))
    except Exception as e:
        q.put(("error", str(e), info))


def app_icon(size):
    """The window icon: a substituted aromatic ring (white) on the accent colour."""
    import math
    from PIL import Image, ImageDraw
    k = 4                                   # drawn 4× and scaled down for smooth edges
    s = size * k
    im = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    d.rounded_rectangle((0, 0, s - 1, s - 1), radius=int(s * 0.22), fill=ACCENT)
    cx = cy = s / 2
    r = s * 0.30
    pts = [(cx + r * math.cos(math.radians(90 + 60 * i)),
            cy - r * math.sin(math.radians(90 + 60 * i))) for i in range(6)]
    width = max(int(s * 0.07), k)
    d.line(pts + pts[:1], fill=PAPER, width=width, joint="curve")
    ri = r * 0.55                            # aromatic circle
    d.ellipse((cx - ri, cy - ri, cx + ri, cy + ri), outline=PAPER, width=max(width // 2, k))
    dot = s * 0.085                          # the substituted position
    x, y = pts[0]
    d.ellipse((x - dot, y - dot, x + dot, y + dot), fill="#9cc0ee", outline=PAPER,
              width=max(width // 3, k))
    return im.resize((size, size), Image.LANCZOS)


def read_guide(path):
    """userguide.txt → [(kind, text)], kind in title / section / sub / code / body / blank.

    The file is a '#'-commented text: banner lines (=====) frame the section titles,
    short capitalised lines are sub-headings, aligned or indented examples stay monospaced.
    """
    lines = [l.rstrip("\n") for l in path.read_text(encoding="utf-8").splitlines()]
    text = [l[2:] if l.startswith("# ") else l.lstrip("#") for l in lines]
    blocks, prev_blank = [], True
    for k, t in enumerate(text):
        s = t.strip()
        if set(s) == {"="}:
            continue
        framed = (k > 0 and set(text[k - 1].strip()) == {"="}
                  and k + 1 < len(text) and set(text[k + 1].strip()) == {"="})
        if framed:
            blocks.append(("title" if not blocks else "section", s))
        elif not s:
            if blocks and blocks[-1][0] != "blank":
                blocks.append(("blank", ""))
        elif (t.startswith("    ") and not s.startswith(("·", "-"))) or "   " in s \
                or s.startswith("--"):
            blocks.append(("code", t.rstrip()))
        elif prev_blank and len(s) < 40 and s[0].isupper() and not s.endswith((":", ".")) \
                and "." not in s and not s.startswith("-"):
            blocks.append(("sub", s))
        else:
            blocks.append(("body", t.rstrip()))
        prev_blank = not s
    return blocks


class GuideWindow(tk.Toplevel):
    """The user guide inside the program: sections on the left, text on the right, find."""

    def __init__(self, app):
        super().__init__(app)
        self.app = app
        self.title("Pólya Pipeline — User guide")
        self.configure(bg=PANEL)
        px = app.px
        w, h = int(app.winfo_width() * 0.82), int(app.winfo_height() * 0.88)
        self.geometry(f"{w}x{h}+{app.winfo_rootx() + (app.winfo_width() - w) // 2}"
                      f"+{app.winfo_rooty() + px(20)}")
        self.minsize(px(640), px(420))
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        # Find bar
        bar = tk.Frame(self, bg=PANEL, padx=px(16), pady=px(10))
        bar.grid(row=0, column=0, columnspan=3, sticky="ew")
        tk.Label(bar, text="User guide", bg=PANEL, fg=INK, font=app.f_head).pack(side="left")
        self.find_var = tk.StringVar()
        self.find_status = tk.Label(bar, bg=PANEL, fg=INK3, font=app.f_label)
        self.find_status.pack(side="right", padx=(px(8), 0))
        find = app._entry(bar, self.find_var, width=26)
        find.unbind("<Return>")
        find.bind("<Return>", lambda e: self._find())
        find.pack(side="right", ipady=px(3))
        tk.Label(bar, text="Find  (Ctrl+F)", bg=PANEL, fg=INK2,
                 font=app.f_label).pack(side="right", padx=(0, px(8)))
        self.find_entry = find

        # Sections
        self.toc = tk.Listbox(self, bg=PANEL, fg=INK, font=app.f_body, relief="flat", bd=0,
                              highlightthickness=0, activestyle="none",
                              selectbackground=SELECT, selectforeground=INK)
        self.toc.grid(row=1, column=0, sticky="ns", padx=(px(12), 0), pady=(0, px(12)))
        self.toc.bind("<<ListboxSelect>>", self._jump)

        # Text
        self.text = tk.Text(self, bg=PAPER, fg=INK, font=app.f_body, wrap="word",
                            relief="flat", bd=0, padx=px(28), pady=px(20),
                            highlightthickness=0, selectbackground=SELECT,
                            selectforeground=INK, spacing1=px(1), spacing3=px(1))
        self.text.grid(row=1, column=1, sticky="nsew", padx=(px(12), 0), pady=(0, px(12)))
        sb = ttk.Scrollbar(self, orient="vertical", command=self.text.yview)
        sb.grid(row=1, column=2, sticky="ns", pady=(0, px(12)), padx=(0, px(12)))
        self.text.configure(yscrollcommand=sb.set)
        t = self.text
        t.tag_configure("title", font=(UI, 17, "bold"), spacing3=px(10))
        t.tag_configure("section", font=(UI, 13, "bold"), spacing1=px(18), spacing3=px(6))
        t.tag_configure("sub", font=(UI, 10, "bold"), foreground=INK, spacing1=px(8),
                        spacing3=px(2))
        t.tag_configure("code", font=app.f_log, foreground=INK, background=PANEL,
                        lmargin1=px(12), lmargin2=px(48))   # a wrapped line stays indented
        t.tag_configure("body", foreground=INK)
        t.tag_configure("found", background="#ffe08a")

        self._sections = []
        for kind, line in read_guide(USER_GUIDE):
            if kind == "section":
                mark = f"sec{len(self._sections)}"
                t.mark_set(mark, "end-1c")
                t.mark_gravity(mark, "left")
                self._sections.append(mark)
                self.toc.insert("end", line)
            t.insert("end", line + "\n", kind)
        t.configure(state="disabled")
        self.toc.configure(width=max((len(self.toc.get(i)) for i in range(self.toc.size())),
                                     default=20) + 2)

        self.bind("<Escape>", lambda e: self.destroy())
        self.bind("<Control-f>", lambda e: (self.find_entry.focus_set(), "break"))
        self.bind("<F1>", lambda e: self.lift())
        self.text.focus_set()

    def _jump(self, _=None):
        sel = self.toc.curselection()
        if sel:
            self.text.yview(self._sections[sel[0]])

    def _find(self):
        word = self.find_var.get().strip()
        t = self.text
        t.tag_remove("found", "1.0", "end")
        if not word:
            self.find_status.configure(text="")
            return
        # Enter again with the same word → next match (wrapping to the top)
        start = t.index("insert +1c") if self._last == word else "1.0"
        pos = t.search(word, start, nocase=True, stopindex="end") or \
            t.search(word, "1.0", nocase=True, stopindex="end")
        self._last = word
        if not pos:
            self.find_status.configure(text="No match")
            return
        t.tag_add("found", pos, f"{pos}+{len(word)}c")
        t.mark_set("insert", pos)
        t.see(pos)
        self.find_status.configure(text="Enter = next match")

    _last = None


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"Pólya Pipeline {der.__version__}")
        self.configure(bg=PANEL)
        self._icons = [ImageTk.PhotoImage(app_icon(n)) for n in (16, 24, 32, 48, 64, 256)]
        self.iconphoto(True, *self._icons)          # also used by the guide window
        self._scale = self.winfo_fpixels("1i") / 96      # 1.25 at 125 %
        self.form_w = self.px(300)                       # width of the form fields
        self._proc = None
        self._queue = None
        self._t0 = 0.0
        self._result_inputs = None    # inputs of the results on screen
        self._preview_job = None
        self._preview_photo = None
        self._numbered = None
        self._row_positions = {}      # table item → "Positions" text of the engine row
        self._result_parent = None    # parent SMILES and substituents of the run on screen
        self._result_subst = []
        self._drawn_smiles = None
        self._guide = None

        self._fonts()
        self._style()
        self._build()

        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        w, h = min(self.px(1180), sw - self.px(40)), min(self.px(800), sh - self.px(90))
        self.geometry(f"{w}x{h}+{max((sw - w) // 2, 0)}+{self.px(10)}")
        self.minsize(min(self.px(900), w), min(self.px(600), h))

        self.bind("<Escape>", lambda e: self._escape())
        self.bind("<F1>", lambda e: self._open_guide())
        self.protocol("WM_DELETE_WINDOW", self._close)
        self._inputs_changed()
        self._fit_form()

    def px(self, n):
        """Logical pixels (designed at 96 dpi) → screen pixels."""
        return int(round(n * self._scale))

    # ── Fonts and ttk styles ──────────────────────────────────────────────────
    def _fonts(self):
        self.f_title = tkfont.Font(family=UI, size=15, weight="bold")
        self.f_group = tkfont.Font(family=UI, size=10, weight="bold")
        self.f_label = tkfont.Font(family=UI, size=9)
        self.f_hint  = tkfont.Font(family=UI, size=8)
        self.f_body  = tkfont.Font(family=UI, size=10)
        self.f_head  = tkfont.Font(family=UI, size=13, weight="bold")
        self.f_mono  = tkfont.Font(family=MONO, size=10)
        self.f_log   = tkfont.Font(family=MONO, size=9)

    def _style(self):
        s = ttk.Style(self)
        s.theme_use("clam")
        s.configure(".", background=PANEL, foreground=INK, font=self.f_body,
                    bordercolor=RULE, lightcolor=PANEL, darkcolor=PANEL,
                    troughcolor=PAPER, focuscolor=ACCENT)
        s.configure("TFrame", background=PANEL)
        s.configure("Paper.TFrame", background=PAPER)

        s.configure("TCheckbutton", background=PANEL, foreground=INK,
                    indicatorbackground=PAPER, indicatorforeground=PAPER,
                    indicatormargin=(0, 0, self.px(6), 0))
        s.map("TCheckbutton", background=[("active", PANEL)])
        # Check mark instead of the theme's ✕ (which reads as "remove")
        self._check_imgs = self._checkbox_images()
        s.element_create("Check.indicator", "image", self._check_imgs["off"],
                         ("selected", self._check_imgs["on"]),
                         ("active", "!selected", self._check_imgs["hover"]))
        s.layout("TCheckbutton", [("Checkbutton.padding", {"sticky": "nswe", "children": [
            ("Check.indicator", {"side": "left", "sticky": ""}),
            ("Checkbutton.focus", {"side": "left", "sticky": "w", "children": [
                ("Checkbutton.label", {"sticky": "nswe"})]})]})])

        s.configure("Primary.TButton", background=ACCENT, foreground=PAPER,
                    font=(UI, 10, "bold"), padding=(self.px(16), self.px(8)),
                    borderwidth=0, focuscolor=PAPER)
        s.map("Primary.TButton",
              background=[("pressed", ACCENT_H), ("active", ACCENT_H),
                          ("disabled", DISABLED)],
              foreground=[("disabled", INK3)])
        s.configure("TButton", background=PAPER, foreground=INK,
                    padding=(self.px(10), self.px(6)), borderwidth=1,
                    bordercolor=RULE, focuscolor=ACCENT)
        s.map("TButton", background=[("pressed", SELECT), ("active", SELECT)],
              foreground=[("disabled", INK3)])
        s.configure("Link.TButton", background=PANEL, foreground=ACCENT,
                    padding=(0, self.px(2)), borderwidth=0, font=(UI, 9))
        s.map("Link.TButton", background=[("active", PANEL), ("pressed", PANEL)],
              foreground=[("active", ACCENT_H)])
        s.configure("PaperLink.TButton", background=PAPER, foreground=ACCENT,
                    padding=(0, 0), borderwidth=0, font=(UI, 9))
        s.map("PaperLink.TButton", background=[("active", PAPER), ("pressed", PAPER)],
              foreground=[("active", ACCENT_H)])

        for orient in ("Vertical", "Horizontal"):
            s.configure(f"{orient}.TScrollbar", background=PANEL, troughcolor=PAPER,
                        bordercolor=PAPER, arrowcolor=INK2, lightcolor=PANEL,
                        darkcolor=PANEL, gripcount=0)
            s.map(f"{orient}.TScrollbar", background=[("active", RULE)])

        s.configure("TNotebook", background=PAPER, borderwidth=0,
                    tabmargins=(0, 0, 0, 0))
        s.configure("TNotebook.Tab", background=PANEL, foreground=INK2,
                    padding=(self.px(14), self.px(5)), bordercolor=RULE)
        s.map("TNotebook.Tab", background=[("selected", PAPER)],
              foreground=[("selected", INK)])

        row_h = int(self.f_body.metrics("linespace") * 1.55)
        s.configure("Treeview", background=PAPER, fieldbackground=PAPER,
                    foreground=INK, rowheight=row_h, borderwidth=0, font=self.f_body)
        s.map("Treeview", background=[("selected", SELECT)],
              foreground=[("selected", INK)])
        s.configure("Treeview.Heading", background=PANEL, foreground=INK,
                    font=(UI, 9, "bold"), relief="flat",
                    padding=(self.px(6), self.px(4)))
        s.map("Treeview.Heading", background=[("active", SELECT)])

        s.configure("Horizontal.TProgressbar", background=ACCENT, troughcolor=PANEL,
                    bordercolor=PANEL, lightcolor=ACCENT, darkcolor=ACCENT)

    def _checkbox_images(self):
        """Box and check-mark images at screen resolution (drawn 4× and scaled down)."""
        from PIL import Image, ImageDraw
        box, gap, k = self.px(16), self.px(8), 4
        imgs = {}
        for name, fill, border, tick in (("off", PAPER, INK3, False),
                                         ("hover", PAPER, ACCENT, False),
                                         ("on", ACCENT, ACCENT, True)):
            im = Image.new("RGBA", ((box + gap) * k, box * k), (0, 0, 0, 0))
            d = ImageDraw.Draw(im)
            d.rounded_rectangle((k // 2, k // 2, box * k - k // 2, box * k - k // 2),
                                radius=3 * k, fill=fill, outline=border, width=max(k, k * box // 14))
            if tick:
                s = box * k
                d.line([(0.24 * s, 0.52 * s), (0.42 * s, 0.70 * s), (0.76 * s, 0.32 * s)],
                       fill=PAPER, width=int(0.13 * s), joint="curve")
            imgs[name] = ImageTk.PhotoImage(im.resize((box + gap, box), Image.LANCZOS))
        return imgs

    # ── Widget helpers ────────────────────────────────────────────────────────
    def _entry(self, parent, var, mono=False, width=None):
        e = tk.Entry(parent, textvariable=var, font=self.f_mono if mono else self.f_body,
                     bg=PAPER, fg=INK, insertbackground=INK, relief="flat", bd=0,
                     highlightthickness=1, highlightbackground=RULE,
                     highlightcolor=ACCENT, disabledbackground=DISABLED,
                     disabledforeground=INK3, selectbackground=SELECT,
                     selectforeground=INK)
        if width:
            e.configure(width=width)
        e.bind("<Return>", lambda ev: self._run())
        return e

    def _label(self, parent, text, row, pady=(6, 1)):
        tk.Label(parent, text=text, bg=PANEL, fg=INK2, font=self.f_label,
                 anchor="w").grid(row=row, column=0, columnspan=2, sticky="ew", pady=pady)

    def _hint(self, parent, row, text=""):
        h = tk.Label(parent, text=text, bg=PANEL, fg=INK3, font=self.f_hint,
                     anchor="w", justify="left", wraplength=self.form_w)
        h.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(2, 0))
        h.default = text
        return h

    def _group(self, parent, text, row):
        tk.Label(parent, text=text, bg=PANEL, fg=INK, font=self.f_group,
                 anchor="w").grid(row=row, column=0, columnspan=2, sticky="ew",
                                  pady=(self.px(12), 0))

    # ── Layout ────────────────────────────────────────────────────────────────
    def _build(self):
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        self._build_form()
        tk.Frame(self, bg=RULE, width=1).grid(row=0, column=1, sticky="nsw")
        self._build_results()

    def _build_form(self):
        # The fields scroll when the window is short; Run and the status stay visible
        outer = tk.Frame(self, bg=PANEL)
        outer.grid(row=0, column=0, sticky="ns")
        outer.rowconfigure(0, weight=1)
        col_w = self.form_w + 2 * self.px(20)
        self.form_canvas = canvas = tk.Canvas(outer, bg=PANEL, highlightthickness=0, bd=0,
                                              width=col_w)
        canvas.grid(row=0, column=0, sticky="ns")
        form_sb = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=form_sb.set)
        f = tk.Frame(canvas, bg=PANEL, padx=self.px(20), pady=self.px(16))
        canvas.create_window(0, 0, window=f, anchor="nw", width=col_w)
        f.columnconfigure(0, weight=1)

        def fit(_=None, form=f):
            canvas.configure(scrollregion=(0, 0, col_w, form.winfo_reqheight()))
            if form.winfo_reqheight() > canvas.winfo_height() > 1:
                form_sb.place(in_=canvas, relx=1.0, rely=0, relheight=1.0, anchor="ne")
            else:
                form_sb.place_forget()
                canvas.yview_moveto(0)
        f.bind("<Configure>", fit)
        canvas.bind("<Configure>", fit)
        self._fit_form = lambda: self.after_idle(fit)
        wheel = lambda e: canvas.yview_scroll(int(-e.delta / 120), "units") \
            if form_sb.winfo_ismapped() else None
        canvas.bind("<Enter>", lambda e: canvas.bind_all("<MouseWheel>", wheel))
        canvas.bind("<Leave>", lambda e: canvas.unbind_all("<MouseWheel>"))

        footer = tk.Frame(outer, bg=PANEL, padx=self.px(20), pady=self.px(14))
        footer.grid(row=1, column=0, columnspan=2, sticky="ew")
        footer.columnconfigure(0, weight=1)
        tk.Frame(outer, bg=RULE, height=1).grid(row=1, column=0, columnspan=2, sticky="new")
        r = 0

        tk.Label(f, text="Pólya Pipeline", bg=PANEL, fg=INK, font=self.f_title,
                 anchor="w").grid(row=r, column=0, columnspan=2, sticky="ew"); r += 1
        tk.Label(f, text=f"Aromatic derivatives by substitution · v{der.__version__}",
                 bg=PANEL, fg=INK2,
                 font=self.f_label, anchor="w").grid(row=r, column=0, columnspan=2,
                                                     sticky="ew"); r += 1

        # ── Molecule ──────────────────────────────────────────────────────────
        self._group(f, "Molecule", r); r += 1
        self._label(f, "SMILES of the parent molecule", r); r += 1
        self.smiles_var = tk.StringVar(value="c1ccccc1")
        self.smiles_entry = self._entry(f, self.smiles_var, mono=True)
        self.smiles_entry.grid(row=r, column=0, columnspan=2, sticky="ew",
                               ipady=self.px(4)); r += 1
        self.smiles_hint = self._hint(
            f, r, "Tip: write [cH:1], [cH:2]… to choose the positions yourself."); r += 1

        self._label(f, "Name for the files", r); r += 1
        self.name_var = tk.StringVar()
        self._entry(f, self.name_var).grid(row=r, column=0, columnspan=2, sticky="ew",
                                          ipady=self.px(4)); r += 1
        self._hint(f, r, "Empty = molecular formula."); r += 1

        # ── Substitution ──────────────────────────────────────────────────────
        self._group(f, "Substitution", r); r += 1
        self._label(f, "Atoms that replace a ring CH", r); r += 1
        self.subst_var = tk.StringVar(value="N")
        self._entry(f, self.subst_var, mono=True).grid(
            row=r, column=0, columnspan=2, sticky="ew", ipady=self.px(4)); r += 1
        self.subst_hint = self._hint(f, r, "Comma-separated, e.g. N · P · N,P · [O+] · [N+]"); r += 1

        self._label(f, "Max. substituents per derivative", r); r += 1
        self.maxsubst_var = tk.StringVar(value="")
        self._entry(f, self.maxsubst_var, mono=True, width=6).grid(
            row=r, column=0, sticky="w", ipady=self.px(4)); r += 1
        self.max_hint = self._hint(f, r, "Empty = no limit. Use it for large molecules."); r += 1

        # ── Output files ──────────────────────────────────────────────────────
        self._group(f, "Output files", r); r += 1
        self.xyz_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(f, text=".xyz geometries", variable=self.xyz_var).grid(
            row=r, column=0, columnspan=2, sticky="w", pady=(self.px(6), 0)); r += 1
        self.gjf_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, text="Gaussian 16 inputs (.gjf)", variable=self.gjf_var,
                        command=self._toggle_route).grid(
            row=r, column=0, columnspan=2, sticky="w", pady=(self.px(4), 0)); r += 1
        self.route_var = tk.StringVar(value=der.DEFAULT_ROUTE)
        self.route_entry = self._entry(f, self.route_var, mono=True)
        self.route_entry.configure(state="disabled")
        self.route_entry.grid(row=r, column=0, columnspan=2, sticky="ew",
                              padx=(self.px(24), 0), pady=(self.px(4), 0),
                              ipady=self.px(3)); r += 1
        self.route_hint = self._hint(f, r, "Route line of every .gjf; must start with #.")
        self.route_hint.grid_configure(padx=(self.px(24), 0)); r += 1

        # ── Advanced (hidden by default) ──────────────────────────────────────
        self.adv_btn = ttk.Button(f, text="Show advanced options", style="Link.TButton",
                                  command=self._toggle_advanced, takefocus=True)
        self.adv_btn.grid(row=r, column=0, columnspan=2, sticky="w",
                          pady=(self.px(12), 0)); r += 1
        adv = self.adv_frame = tk.Frame(f, bg=PANEL)
        adv.columnconfigure(0, weight=1)
        self._adv_row = r; r += 1
        self._label(adv, "Symmetry tolerance (Å)", 0, pady=(2, 2))
        self.tol_var = tk.StringVar(value="0.15")
        self._entry(adv, self.tol_var, mono=True, width=6).grid(
            row=1, column=0, sticky="w", ipady=self.px(4))
        self.tol_hint = self._hint(
            adv, 2, "How far an atom may be from its symmetric position in the 3D check. "
                    "Raise it (e.g. 0.20) if that check disagrees.")
        self.charges_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(adv, text="Distinguish formal charges",
                        variable=self.charges_var).grid(
            row=3, column=0, columnspan=2, sticky="w", pady=(self.px(10), 0))
        self._hint(adv, 4, "Treats e.g. N and [N+] as different atoms when finding the symmetry.")

        # ── Actions (fixed footer) ────────────────────────────────────────────
        f = footer
        r = 0
        self.run_btn = ttk.Button(f, text="Run", style="Primary.TButton",
                                  command=self._run)
        self.run_btn.grid(row=r, column=0, columnspan=2, sticky="ew",
                          pady=(0, self.px(8))); r += 1
        links = tk.Frame(f, bg=PANEL)
        links.grid(row=r, column=0, columnspan=2, sticky="ew"); r += 1
        links.columnconfigure((0, 1), weight=1, uniform="b")
        ttk.Button(links, text="All results", command=self._open_output).grid(
            row=0, column=0, sticky="ew", padx=(0, self.px(4)))
        ttk.Button(links, text="User guide  (F1)", command=self._open_guide).grid(
            row=0, column=1, sticky="ew", padx=(self.px(4), 0))
        self.status_var = tk.StringVar()
        self.status = tk.Label(f, textvariable=self.status_var, bg=PANEL, fg=INK2,
                               font=self.f_label, anchor="w", justify="left",
                               wraplength=self.form_w)
        self.status.grid(row=r, column=0, columnspan=2, sticky="ew", pady=(self.px(10), 0))

        for var in (self.smiles_var, self.name_var, self.subst_var, self.maxsubst_var,
                    self.tol_var, self.route_var, self.xyz_var, self.gjf_var,
                    self.charges_var):
            var.trace_add("write", self._inputs_changed)

    def _build_results(self):
        right = tk.Frame(self, bg=PAPER, padx=self.px(24), pady=self.px(16))
        right.grid(row=0, column=1, sticky="nsew", padx=(1, 0))
        right.columnconfigure(0, weight=1)
        right.rowconfigure(3, weight=1)

        # Molecule, with the substitution positions numbered as in the table
        self.preview_box = tk.Frame(right, bg=PAPER, height=self.px(230))
        self.preview_box.grid(row=0, column=0, sticky="ew")
        self.preview_box.grid_propagate(False)
        self.preview_box.pack_propagate(False)
        self.preview = tk.Label(self.preview_box, bg=PAPER, fg=INK3, font=self.f_body,
                                wraplength=self.px(520), justify="center")
        self.preview.pack(fill="both", expand=True)
        self.preview_box.bind("<Configure>", lambda e: self._schedule_preview(150))
        # The drawing takes ~28 % of the height, so the table keeps room on small windows
        right.bind("<Configure>", lambda e: self.preview_box.configure(
            height=max(self.px(140), min(self.px(230), int(e.height * 0.28)))))

        cap = tk.Frame(right, bg=PAPER)
        cap.grid(row=1, column=0, sticky="ew", pady=(self.px(2), 0))
        self.caption = tk.Label(cap, bg=PAPER, fg=INK2, font=self.f_label, anchor="w")
        self.caption.pack(side="left")
        self.use_numbered_btn = ttk.Button(cap, text="Use numbered SMILES",
                                           style="PaperLink.TButton",
                                           command=self._use_numbered)
        self.show_parent_btn = ttk.Button(cap, text="Show the parent molecule  (Esc)",
                                          style="PaperLink.TButton",
                                          command=self._show_parent)

        # Summary of the run: the answer comes first, the log is the appendix
        tk.Frame(right, bg=RULE, height=1).grid(row=2, column=0, sticky="new",
                                                pady=(self.px(12), 0))
        body = tk.Frame(right, bg=PAPER)
        body.grid(row=3, column=0, sticky="nsew", pady=(self.px(12), 0))
        body.columnconfigure(0, weight=1)
        body.rowconfigure(2, weight=1)

        top = tk.Frame(body, bg=PAPER)
        top.grid(row=0, column=0, sticky="ew")
        top.columnconfigure(0, weight=1)
        self.headline = tk.Label(top, bg=PAPER, fg=INK, font=self.f_head, anchor="w")
        self.headline.grid(row=0, column=0, sticky="ew")
        self.detail = tk.Label(top, bg=PAPER, fg=INK2, font=self.f_body, anchor="w",
                               justify="left")
        self.detail.grid(row=1, column=0, sticky="ew", pady=(self.px(2), 0))
        def wrap_summary(e):
            # room left by the two buttons, whether or not they are shown yet
            buttons = self.open_files_btn.winfo_reqwidth() + self.copy_btn.winfo_reqwidth()
            width = max(e.width - buttons - self.px(28), self.px(220))
            self.headline.configure(wraplength=width)
            self.detail.configure(wraplength=width)
        top.bind("<Configure>", wrap_summary)
        self.headline.configure(justify="left")
        self.progress = ttk.Progressbar(top, mode="indeterminate", length=self.px(160))
        btns = tk.Frame(top, bg=PAPER)
        btns.grid(row=0, column=1, rowspan=2, sticky="ne", padx=(self.px(12), 0))
        self.open_files_btn = ttk.Button(btns, text="Open this run's files",
                                         command=self._open_run_files)
        self.copy_btn = ttk.Button(btns, text="Copy SMILES", command=self._copy_smiles)
        self._run_dir = None

        self.tabs = ttk.Notebook(body)
        self.tabs.grid(row=2, column=0, sticky="nsew", pady=(self.px(12), 0))
        self.table_note = tk.Label(body, bg=PAPER, fg=INK3, font=self.f_hint, anchor="w",
                                   justify="left")
        self.table_note.grid(row=3, column=0, sticky="ew", pady=(self.px(6), 0))
        body.bind("<Configure>", lambda e: self.table_note.configure(wraplength=e.width))

        # Before the first run: examples that fill in the form (instead of an empty table)
        self.examples = tk.Frame(body, bg=PAPER)
        tk.Label(self.examples, text="Try an example:", bg=PAPER, fg=INK2,
                 font=self.f_label).pack(side="left", padx=(0, self.px(10)))
        for label, smi, subst, name in (("Pyridine + N", "n1ccccc1", "N", "pyridine"),
                                        ("Naphthalene + N, P", "c1ccc2ccccc2c1", "N,P",
                                         "naphthalene"),
                                        ("Toluene + N", "Cc1ccccc1", "N", "toluene")):
            ttk.Button(self.examples, text=label,
                       command=lambda s=smi, u=subst, n=name: self._example(s, u, n)).pack(
                side="left", padx=(0, self.px(8)))

        # Derivatives table
        tab1 = ttk.Frame(self.tabs, style="Paper.TFrame")
        tab1.columnconfigure(0, weight=1)
        tab1.rowconfigure(0, weight=1)
        self.tree = ttk.Treeview(tab1, show="headings", selectmode="extended")
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb = ttk.Scrollbar(tab1, orient="vertical", command=self.tree.yview)
        vsb.grid(row=0, column=1, sticky="ns")
        hsb = ttk.Scrollbar(tab1, orient="horizontal", command=self.tree.xview)
        hsb.grid(row=1, column=0, sticky="ew")
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.tag_configure("invalid", foreground=WARN)
        self.tree.tag_configure("stale", foreground=INK3)
        self.tree.bind("<Control-c>", lambda e: self._copy_smiles())
        self.tree.bind("<<TreeviewSelect>>", self._on_row_selected)
        self.tabs.add(tab1, text="Derivatives")

        # Full log of the run
        tab2 = ttk.Frame(self.tabs, style="Paper.TFrame")
        tab2.columnconfigure(0, weight=1)
        tab2.rowconfigure(0, weight=1)
        self.log = tk.Text(tab2, bg=PAPER, fg=INK, font=self.f_log, wrap="none",
                           relief="flat", bd=0, padx=self.px(10), pady=self.px(8),
                           selectbackground=SELECT, selectforeground=INK,
                           state="disabled", highlightthickness=0)
        self.log.grid(row=0, column=0, sticky="nsew")
        lv = ttk.Scrollbar(tab2, orient="vertical", command=self.log.yview)
        lv.grid(row=0, column=1, sticky="ns")
        lh = ttk.Scrollbar(tab2, orient="horizontal", command=self.log.xview)
        lh.grid(row=1, column=0, sticky="ew")
        self.log.configure(yscrollcommand=lv.set, xscrollcommand=lh.set)
        self.log.tag_configure("ok", foreground=OK)
        self.log.tag_configure("warn", foreground=WARN)
        self.tabs.add(tab2, text="Log")

        self._show_empty_state()

    # ── Results area states ───────────────────────────────────────────────────
    def _set_summary(self, headline, detail, detail_fg=INK2):
        self.headline.configure(text=headline)
        self.detail.configure(text=detail, fg=detail_fg)

    def _show_empty_state(self):
        self._set_summary(
            "No results yet",
            "Enter a SMILES and the atoms that replace a ring CH, then press Run. "
            "The number of derivatives, their list and the files will appear here. "
            "To process many molecules at once, use the command line (see the user guide).")
        self.tree.configure(columns=())
        self._show_result_buttons(False)
        self.tabs.grid_remove()
        self.table_note.configure(text="")
        self.examples.grid(row=2, column=0, sticky="nw", pady=(self.px(16), 0))

    def _show_tabs(self):
        self.examples.grid_remove()
        self.tabs.grid()

    def _example(self, smiles, subst, name):
        self.smiles_var.set(smiles)
        self.subst_var.set(subst)
        self.name_var.set(name)
        self.maxsubst_var.set("")
        self.status_var.set("Example filled in. Press Enter to run.")
        self.run_btn.focus_set()

    def _show_result_buttons(self, show, copy=True):
        for b in (self.open_files_btn, self.copy_btn):
            b.pack_forget()
        if show:
            if self._run_dir:
                self.open_files_btn.pack(side="left", padx=(0, self.px(6)))
            if copy:
                self.copy_btn.pack(side="left")

    # ── Preview ───────────────────────────────────────────────────────────────
    def _schedule_preview(self, delay=350):
        if self._preview_job is not None:
            self.after_cancel(self._preview_job)
        self._preview_job = self.after(delay, self._draw)

    def _draw(self):
        """Draws the selected derivative, or the parent molecule with its positions numbered."""
        self._preview_job = None
        self.use_numbered_btn.pack_forget()
        self.show_parent_btn.pack_forget()
        selected = self.tree.focus() or next(iter(self.tree.selection()), "")
        if selected and selected in self.tree.selection() and selected in self._row_positions:
            self._draw_derivative(selected)
            return
        smi = self.smiles_var.get().strip()
        problem = der.smiles_problem(smi)
        if problem:
            self._preview_photo = None
            self.preview.configure(image="", fg=INK3,
                                   text="No drawing: the reason is shown under the SMILES."
                                   if smi else "The molecule will be drawn here.")
            self.caption.configure(text="")
            return
        w = max(self.preview_box.winfo_width(), self.px(200))
        h = max(self.preview_box.winfo_height(), self.px(120))
        img = der.smiles_image(smi, size=(w, h),
                               label_colour=hex_to_rgb(ACCENT))
        self._preview_photo = ImageTk.PhotoImage(img)   # keep a reference
        self.preview.configure(image=self._preview_photo, text="")

        mol = der.Chem.MolFromSmiles(smi)
        self._numbered = der.number_ring_positions(mol)
        positions, scheme = der.substitution_positions(mol)
        n = len(positions)
        what = f"{n} ring positions" if self._numbered else f"{n} substitution positions"
        how = {"iupac": "IUPAC numbering",
               "user": "numbered by you in the SMILES",
               "program": "program numbering, not IUPAC (this ring system is not in "
                          "the program's IUPAC table)"}[scheme]
        kept = "; CH outside the rings are kept" if self._numbered else ""
        self.caption.configure(text=f"{what}, {how}{kept}.   ")
        if self._numbered:
            self.use_numbered_btn.pack(side="left")

    def _preview_size(self):
        return (max(self.preview_box.winfo_width(), self.px(200)),
                max(self.preview_box.winfo_height(), self.px(120)))

    def _draw_derivative(self, item):
        """The selected row's derivative, in the parent's layout and numbering."""
        num = self.tree.set(item, "#")
        positions = self._row_positions[item]
        where = positions.replace(f"{self._result_subst[0]}: ", "") \
            if len(self._result_subst) == 1 else positions
        try:
            img, smi = der.derivative_image(self._result_parent, positions, self._result_subst,
                                            self._preview_size(),
                                            label_colour=hex_to_rgb(ACCENT))
        except ValueError as e:
            self._preview_photo = None
            self.preview.configure(image="", fg=WARN, text=f"Derivative {num}: {e}")
            self.caption.configure(text=f"Derivative {num}, substituted at {where}.   ")
        else:
            self._preview_photo = ImageTk.PhotoImage(img)
            self.preview.configure(image=self._preview_photo, text="")
            what = "the parent molecule" if positions == "—" else \
                f"substituted at {where} (shaded)"
            self.caption.configure(text=f"Derivative {num}: {what} · {smi}   ")
        self.show_parent_btn.pack(side="left")

    def _on_row_selected(self, _=None):
        self._schedule_preview(30)

    def _show_parent(self):
        self.tree.selection_remove(*self.tree.selection())
        self._schedule_preview(0)

    def _escape(self):
        if self._proc is not None:
            self._cancel()
        elif self.tree.selection():
            self._show_parent()

    def _use_numbered(self):
        if self._numbered:
            self.smiles_var.set(self._numbered[1])

    # ── Live validation ───────────────────────────────────────────────────────
    def _inputs(self):
        return (self.smiles_var.get().strip(), self.name_var.get().strip(),
                self.subst_var.get().strip(), self.maxsubst_var.get().strip(),
                self.tol_var.get().strip(), self.xyz_var.get(), self.gjf_var.get(),
                self.route_var.get().strip(), self.charges_var.get())

    def _field_errors(self):
        """{hint label: message} for every field that blocks a run."""
        errors = {}
        smi, _, subst, ms, tol, _, gjf, route, _ = self._inputs()
        problem = der.smiles_problem(smi)
        if problem:
            errors[self.smiles_hint] = problem
        atoms = [s.strip() for s in subst.split(",") if s.strip()]
        if not atoms:
            errors[self.subst_hint] = "Enter at least one atom, e.g. N."
        for a in atoms:
            het = der.read_heteroatom(a)
            if het is None:
                errors[self.subst_hint] = f"'{a}' is not a one-atom SMILES."
                break
            if het[0] in (1, 9, 17, 35, 53) and het[1] == 0:
                errors[self.subst_hint] = (f"{a} makes only one bond, so it cannot replace "
                                           "a ring CH. Use atoms such as N, P or [O+].")
                break
        if ms and not (ms.isdigit() and int(ms) > 0):
            errors[self.max_hint] = "Must be a whole number of 1 or more (or empty)."
        try:
            if float(tol) <= 0:
                raise ValueError
        except ValueError:
            errors[self.tol_hint] = "Must be a positive number of ångström, e.g. 0.15."
        if gjf and not route.startswith("#"):
            errors[self.route_hint] = "The route line must start with #, e.g. #b3lyp/6-31G(d) opt."
        return errors

    def _inputs_changed(self, *_):
        smi = self.smiles_var.get().strip()
        if smi != self._drawn_smiles:     # a new parent: show it, not an old derivative
            self._drawn_smiles = smi
            if self.tree.selection():
                self.tree.selection_remove(*self.tree.selection())
        self._schedule_preview()
        errors = self._field_errors()
        for hint in (self.smiles_hint, self.subst_hint, self.max_hint,
                     self.tol_hint, self.route_hint):
            hint.configure(text=errors.get(hint, hint.default),
                           fg=WARN if hint in errors else INK3)
        if not errors.keys() & {self.smiles_hint, self.subst_hint, self.max_hint}:
            self._show_job_size()
        if self.tol_hint in errors and not self.adv_frame.winfo_ismapped():
            self._toggle_advanced()
        running = self._proc is not None
        if not running:
            self.run_btn.configure(state="disabled" if errors else "normal")
            self.status_var.set("Fix the field marked in red to run." if errors
                                else "Ready. Press Enter to run.")
        self._mark_stale()

    def _show_job_size(self):
        """Number of substitution patterns the run will check, before it starts."""
        smi, _, subst, ms, *_ = self._inputs()
        n = len(der.substitution_positions(der.Chem.MolFromSmiles(smi))[0])
        nvars = len([s for s in subst.split(",") if s.strip()])
        count = der.count_patterns(n, nvars, int(ms) if ms else None)
        limit = "" if ms else "Empty = no limit. "
        if count > der.LARGE_ENUMERATION:
            self.max_hint.configure(
                fg=STALE, text=f"{limit}This run checks {count:,} patterns and may take "
                               "minutes. Set a maximum to make it faster.")
        else:
            self.max_hint.configure(fg=INK3, text=f"{limit}This run checks {count:,} patterns.")

    def _mark_stale(self):
        if self._result_inputs is None or self._proc is not None:
            return
        stale = self._inputs() != self._result_inputs
        for item in self.tree.get_children():
            tags = [t for t in self.tree.item(item, "tags") if t != "stale"]
            self.tree.item(item, tags=tags + (["stale"] if stale else []))
        if stale:
            self.detail.configure(text="Inputs changed since this run. Press Run to update "
                                       "the results.", fg=STALE)
        else:
            self.detail.configure(text=self._detail_text, fg=self._detail_fg)

    def _toggle_route(self):
        self.route_entry.configure(state="normal" if self.gjf_var.get() else "disabled")

    def _toggle_advanced(self):
        if self.adv_frame.winfo_ismapped():
            self.adv_frame.grid_remove()
            self.adv_btn.configure(text="Show advanced options")
        else:
            self.adv_frame.grid(row=self._adv_row, column=0, columnspan=2, sticky="ew",
                                pady=(self.px(6), 0))
            self.adv_btn.configure(text="Hide advanced options")
        self._fit_form()

    # ── Run / cancel ──────────────────────────────────────────────────────────
    def _run(self):
        if self._proc is not None:
            return
        errors = self._field_errors()
        if errors:
            self.status_var.set("Fix the field marked in red to run.")
            return
        smi, name, subst, ms, tol, xyz, gjf, route, charges = self._inputs()
        max_subst = int(ms) if ms else None

        # Earlier results with the same name and substituents are replaced: ask first
        if name:
            file_name = name
        else:
            from rdkit.Chem import rdMolDescriptors
            file_name = rdMolDescriptors.CalcMolFormula(der.Chem.MolFromSmiles(smi))
        tag = der.safe_filename("-".join(s.strip() for s in subst.split(",") if s.strip()))
        suffix = f"_max{max_subst}" if max_subst is not None else ""
        csv = der.OUTPUT_DIR / "tables" / f"{der.safe_filename(file_name)}_{tag}{suffix}.csv"
        if csv.exists() and not messagebox.askyesno(
                "Replace earlier results?",
                f"Results for “{file_name}” with {subst} already exist in generated_files.\n"
                "Running again replaces that table and its .xyz/.gjf files.\n\nContinue?",
                icon="warning", parent=self):
            return

        kwargs = dict(parent_smiles=smi, subst_text=subst, name=name or None,
                      max_subst=max_subst, tol=float(tol), generate_xyz_files=xyz,
                      respect_charges=charges, generate_gaussian_files=gjf,
                      gaussian_route=route)
        self._started_inputs = self._inputs()
        self._clear_results()
        self._show_tabs()
        self.tabs.select(1)             # follow the log while it runs
        ctx = mp.get_context("spawn")
        self._queue = ctx.Queue()
        self._proc = ctx.Process(target=pipeline_process, daemon=True,
                                 args=(self._queue, kwargs, str(der.OUTPUT_DIR)))
        self._proc.start()
        self._t0 = time.perf_counter()
        self.run_btn.configure(text="Cancel  (Esc)", command=self._cancel, state="normal")
        self._set_summary("Running…", "Counting and enumerating the derivatives. "
                                      "Press Esc to cancel.")
        self.progress.grid(row=0, column=1, sticky="e")
        self.progress.start(12)
        self.status_var.set("Running…")
        self.after(60, self._poll)

    def _clear_results(self):
        self._result_inputs = None
        self._run_dir = None
        self._show_result_buttons(False)
        self.tree.delete(*self.tree.get_children())
        self.tree.configure(columns=())
        self._row_positions = {}
        self._result_parent = None
        self._schedule_preview(0)
        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        self.log.configure(state="disabled")

    def _poll(self):
        if self._proc is None:
            return
        try:
            while True:
                msg = self._queue.get_nowait()
                if msg[0] == "text":
                    self._log(msg[1])
                elif msg[0] == "done":
                    self._finish(rows=msg[1], info=msg[2])
                    return
                else:
                    self._finish(error=msg[1], info=msg[2])
                    return
        except queue.Empty:
            pass
        if not self._proc.is_alive() and self._queue.empty():
            self._finish(error="The run ended unexpectedly.", info={})
            return
        self.after(60, self._poll)

    def _cancel(self):
        if self._proc is None:
            return
        self._proc.terminate()
        self._end_run()
        self._set_summary("Cancelled", "No results. Files written before cancelling are "
                                       "kept in generated_files.")
        self.status_var.set("Cancelled.")

    def _end_run(self):
        self._proc = None
        self.progress.stop()
        self.progress.grid_remove()
        self.run_btn.configure(text="Run", command=self._run)
        self._inputs_changed()

    def _finish(self, rows=None, info=None, error=None):
        elapsed = time.perf_counter() - self._t0
        proc = self._proc
        self._end_run()
        if proc is not None:
            proc.join(timeout=1)
        if error:
            self._set_summary("The run stopped", error, WARN)
            self.status_var.set("Stopped with an error. See the Log tab for details.")
            self.tabs.select(1)
            return
        self._result_parent = self._started_inputs[0]
        self._result_subst = info.get("substituents", [])
        self._show_table(rows, info)
        self._detail_text, self._detail_fg = self._summary_detail(info)
        self._set_summary(self._summary_headline(info), self._detail_text, self._detail_fg)
        self._run_dir = info.get("gjf_dir") or info.get("xyz_dir") or (
            str(Path(info["csv_path"]).parent) if info.get("csv_path") else None)
        self._show_result_buttons(True)
        self._result_inputs = self._started_inputs
        self.status_var.set(f"Done in {elapsed:.1f} s.")
        self.tabs.select(0)
        self._mark_stale()

    # ── Summary text ──────────────────────────────────────────────────────────
    def _summary_headline(self, info):
        n = info["distinct"]
        if info["duplicates"] or info["invalid"]:
            return f"{n} distinct compound{'s' * (n != 1)}"
        return f"{n} derivative{'s' * (n != 1)}, all different compounds"

    def _summary_detail(self, info):
        parts, fg = [], INK2
        if info["duplicates"]:
            parts.append(f"Pólya counts {info['patterns']} substitution patterns; "
                         f"{info['duplicates']} repeat a compound and were merged.")
        elif info.get("max_subst") is not None:
            parts.append(f"Up to {info['max_subst']} substituents "
                         f"(Pólya total for any number: {info['polya_total']}).")
        else:
            parts.append(f"Pólya count: {info['polya_total']}.")
        if info["invalid"]:
            parts.append(f"{info['invalid']} would have an impossible valence: listed in red, "
                         "no files written.")
            fg = WARN
        check = info.get("check_3d")
        if check == "match":
            parts.append("Symmetry checked two ways (3D geometry and graph agree).")
        elif check == "less":
            parts.append(f"The 3D geometry is less symmetric than the graph "
                         f"(|G| {info.get('group_order_3d')} vs {info['group_order']}); "
                         "the graph symmetry is used.")
        elif check == "mismatch":
            parts.append("3D and graph symmetry disagree; the graph symmetry is used. "
                         "Try another tolerance.")
            fg = WARN
        else:
            parts.append("3D check skipped: some positions carry more than one H.")
        if info.get("xyz_dir") or info.get("gjf_dir"):
            kinds = " and ".join(k for k, d in ((".xyz", info.get("xyz_dir")),
                                                (".gjf", info.get("gjf_dir"))) if d)
            parts.append(f"{info['files_written']} {kinds} files written.")
            if info.get("file_failures"):
                parts.append(f"{info['file_failures']} geometries failed (see the Log).")
                fg = WARN
        return " ".join(parts), fg

    # ── Table and log ─────────────────────────────────────────────────────────
    def _show_table(self, rows, info):
        subst = info.get("substituents", [])
        single = len(subst) == 1        # one atom type: its count column and "N:" add nothing
        cols = ["#", "Substituents", *([] if single else subst),
                "Equivalent patterns", "Positions", "SMILES"]
        self.tree.configure(columns=cols)
        widths = {"#": 50, "Substituents": 105, "Equivalent patterns": 140,
                  "Positions": 140, "SMILES": 260}
        for c in cols:
            self.tree.heading(c, text=c, command=lambda c=c: self._sort(c, False))
            self.tree.column(c, width=self.px(widths.get(c, 56)),
                             anchor="w" if c in ("Positions", "SMILES") else "center",
                             stretch=c == "SMILES")
        self._row_positions = {}
        for row in rows:
            num, n_subst, *counts, positions, patterns, smi = row
            shown = positions.replace(f"{subst[0]}: ", "") if single else positions
            if single:
                counts = []
            tags = ["invalid"] if "invalid" in smi else []
            item = self.tree.insert("", "end", tags=tags,
                                    values=(num, n_subst, *counts, patterns, shown, smi))
            self._row_positions[item] = positions
        numbering = {"iupac": "IUPAC numbering", "user": "your numbering from the SMILES",
                     "program": "program numbering (not IUPAC)"}[info.get("numbering", "program")]
        self.table_note.configure(
            text=f"Positions use {numbering}, as on the drawing. Equivalent patterns: how many "
                 "ways of placing the substituents give this same derivative by symmetry.")

    def _sort(self, col, reverse):
        def key(item):
            v = self.tree.set(item, col)
            try:
                return (0, float(v))
            except ValueError:
                return (1, v)
        items = sorted(self.tree.get_children(), key=key, reverse=reverse)
        for i, item in enumerate(items):
            self.tree.move(item, "", i)
        self.tree.heading(col, command=lambda: self._sort(col, not reverse))

    def _copy_smiles(self):
        items = self.tree.selection() or self.tree.get_children()
        smiles = [self.tree.set(i, "SMILES") for i in items]
        smiles = [s for s in smiles if s and "invalid" not in s]
        if smiles:
            self.clipboard_clear()
            self.clipboard_append("\n".join(smiles))
            what = "selected" if self.tree.selection() else "all"
            self.status_var.set(f"Copied {len(smiles)} SMILES ({what}).")

    def _log(self, text):
        self.log.configure(state="normal")
        for line in text.splitlines(keepends=True):
            tag = "ok" if "✓" in line else "warn" if "⚠" in line else ()
            self.log.insert("end", line, tag)
        self.log.configure(state="disabled")

    # ── Files and help ────────────────────────────────────────────────────────
    def _open(self, path):
        try:
            open_path(path)
        except Exception as e:
            messagebox.showerror("Could not open", f"{path}\n\n{e}", parent=self)

    def _open_output(self):
        der.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        self._open(der.OUTPUT_DIR)

    def _open_run_files(self):
        if self._run_dir:
            self._open(self._run_dir)

    def _open_guide(self):
        """Opens the user guide inside the program (one window, brought to the front)."""
        if self._guide is not None and self._guide.winfo_exists():
            self._guide.deiconify()
            self._guide.lift()
            self._guide.focus_force()
            return
        try:
            self._guide = GuideWindow(self)
        except OSError as e:
            messagebox.showerror("User guide not found", f"{USER_GUIDE}\n\n{e}", parent=self)

    def _close(self):
        if self._proc is not None:
            self._proc.terminate()
        self.destroy()


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    mp.freeze_support()
    app = App()
    app.mainloop()
