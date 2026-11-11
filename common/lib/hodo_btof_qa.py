"""Draw CBtof0 (btof0) QA pages from decoded Hodo ROOT for hdphc books."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from . import hodo_btof_monitor, hodo_qa_pdf


#_____________________________________________________________________________
def btof0_pdf_path(run: int, *, kaon: bool = False) -> Path:
    suffix = hodo_qa_pdf.particle_suffix(kaon)
    return hodo_qa_pdf.img_dir(run) / f"run{int(run):05d}_hodo_btof0_{suffix}.pdf"


#_____________________________________________________________________________
def _particle_tag(kaon: bool) -> str:
    return "_K" if kaon else "_Pi"


#_____________________________________________________________________________
def _hist_candidates(kaon: bool) -> list[tuple[str, str]]:
    """(key, title) pairs; prefer particle-suffixed, fall back to unsuffixed."""
    tag = _particle_tag(kaon)
    names = [
        ("CBtof0", "CBtof0 (btof0)"),
        ("CBtof0_vs_deT0Seg", "CBtof0 vs dE(T0)"),
        ("CBtof0_vs_deBtof0Seg", "CBtof0 vs dE(Btof0 seg)"),
    ]
    out = []
    for base, title in names:
        out.append((f"{base}{tag}", title))
        out.append((base, title))
    return out


#_____________________________________________________________________________
def _pick_hist(keys: set[str], kaon: bool) -> list[tuple[str, str]]:
    """Unique hist keys present in the file, preserving preferred order."""
    seen = set()
    picked = []
    for key, title in _hist_candidates(kaon):
        # uproot keys may include ";1"
        match = None
        if key in keys:
            match = key
        else:
            for k in keys:
                if k.split(";")[0] == key:
                    match = k.split(";")[0]
                    break
        if match and match not in seen:
            seen.add(match)
            picked.append((match, title))
    return picked


#_____________________________________________________________________________
def draw_btof0_pdf(
    run: int,
    *,
    kaon: bool = False,
    root_path: Path | None = None,
    out_path: Path | None = None,
) -> tuple[bool, str, Optional[Path]]:
    """Write a small PDF of CBtof0 (+ related 2D). Returns (ok, msg, path)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import uproot

    root = Path(root_path) if root_path is not None else hodo_btof_monitor.hodo_root(run)
    out = Path(out_path) if out_path is not None else btof0_pdf_path(run, kaon=kaon)
    if not root.is_file():
        return False, f"Hodo ROOT missing: {root}", None

    try:
        with uproot.open(root) as rf:
            keys = {k.split(";")[0] for k in rf.keys()}
            picked = _pick_hist(keys, kaon)
            if not picked:
                return False, f"no CBtof0 hist in {root.name}", None
            out.parent.mkdir(parents=True, exist_ok=True)
            with PdfPages(out) as pdf:
                for name, title in picked:
                    obj = rf[name]
                    fig, ax = plt.subplots(figsize=(8.0, 5.5))
                    values = obj.values(flow=False)
                    if values.ndim == 1:
                        edges = obj.axes[0].edges()
                        centers = 0.5 * (edges[:-1] + edges[1:])
                        ax.step(centers, values, where="mid", color="C0")
                        ax.fill_between(
                            centers, values, step="mid", alpha=0.25, color="C0",
                        )
                        ax.set_xlabel("btof0 [ns]")
                        ax.set_ylabel("counts")
                    elif values.ndim == 2:
                        xedges = obj.axes[0].edges()
                        yedges = obj.axes[1].edges()
                        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                        im = ax.imshow(
                            values.T,
                            origin="lower",
                            aspect="auto",
                            extent=extent,
                            interpolation="nearest",
                        )
                        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                        ax.set_xlabel("dE [mip]")
                        ax.set_ylabel("btof0 [ns]")
                    else:
                        plt.close(fig)
                        continue
                    ax.set_title(
                        f"run{int(run):05d} {title} "
                        f"({hodo_qa_pdf.particle_suffix(kaon)})"
                    )
                    fig.tight_layout()
                    pdf.savefig(fig)
                    plt.close(fig)
    except Exception as exc:
        return False, f"btof0 PDF failed: {type(exc).__name__}: {exc}", None

    return True, f"wrote {out}", out


#_____________________________________________________________________________
def ensure_btof0_pdf(
    run: int,
    *,
    kaon: bool = False,
    root_path: Path | None = None,
    force: bool = False,
) -> tuple[bool, str, Optional[Path]]:
    """Return existing btof0 PDF or draw a new one."""
    out = btof0_pdf_path(run, kaon=kaon)
    if out.is_file() and not force:
        root = (
            Path(root_path) if root_path is not None
            else hodo_btof_monitor.hodo_root(run)
        )
        if root.is_file() and out.stat().st_mtime >= root.stat().st_mtime:
            return True, f"reuse {out.name}", out
    return draw_btof0_pdf(
        run, kaon=kaon, root_path=root_path, out_path=out,
    )
