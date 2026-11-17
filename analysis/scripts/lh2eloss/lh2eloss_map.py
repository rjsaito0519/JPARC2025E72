"""Python twin of the analyzer's Lh2ElossMan: map file I/O, Classify and Loss (vectorised).

Same rules as src/Lh2ElossMan.cc of the analyzer (JPARC2025E72):
  region (vertex, then ray):
    VESSEL     r < VesselRadius and |y| < VesselHalfHeight
    GAP        else r < HolderInnerRadius and y > HolderBottomY, ray does not enter the vessel (L = 0)
    GAP_CROSS  as GAP, ray enters the vessel (L = chord through the liquid, may be 0)
    OUTSIDE    anything else (no loss)
  L / exit   chord of the ray through the liquid cylinder and the surface it leaves through (SIDE / CAP)
  holder     where the ray leaves the holder cylinder: WINDOW or G10 (bottom plate and open top count as G10)
  cos        incidence cosine there (bottom plate: |u_y|, open top: 0)
  loss       trilinear in (p, L, cos), clamped to the axis ends, >= 0

Integer codes are the C++ enum values and the numbers written in the map file.
Grids are keyed by (particle, region, exit, holder); a grid is dict(p=, L=, c=, loss=array(np, nL, nc)).

Map file: geometry 'Key value' lines (any order, all required), then per grid the label lines
'particle N' / 'region N' / 'exit N' / 'holder N', the p / L / cos axis lines and the rows
'p L loss(cos_0) ... loss(cos_n)' (p outer, L inner). Every particle in the file has all REQUIRED grids.
"""
import numpy as np

PROTON, PION, KAON = 0, 1, 2
VESSEL, GAP, GAP_CROSS, OUTSIDE = 0, 1, 2, 3
SIDE, CAP = 0, 1
WINDOW, G10 = 0, 1

PARTICLE_NAMES = ("proton", "pion", "kaon")
REGION_NAMES = ("vessel", "gap", "gap_cross")
EXIT_NAMES = ("side", "cap")
HOLDER_NAMES = ("window", "g10")

# (region, exit, holder) every particle must have; GAP has no liquid chord, so only SIDE
REQUIRED = tuple((rg, ex, ho) for rg in (VESSEL, GAP, GAP_CROSS)
                 for ex in ((SIDE,) if rg == GAP else (SIDE, CAP)) for ho in (WINDOW, G10))

GEO_KEYS = ("TargetRadius", "TargetHalfHeight", "TargetCenterX", "TargetCenterY", "TargetCenterZ",
            "VesselRadius", "VesselHalfHeight", "HolderInnerRadius", "HolderRadius",
            "HolderWindowHalfPhi", "HolderWindowHalfHeight", "HolderBottomY")
GEO_UNIT = {k: ("deg" if k == "HolderWindowHalfPhi" else "mm") for k in GEO_KEYS}
EPS = 1e-12


class MapError(ValueError):
    pass


def grid_name(key):
    pa, rg, ex, ho = key
    return f"{PARTICLE_NAMES[pa]} {REGION_NAMES[rg]} {EXIT_NAMES[ex]} {HOLDER_NAMES[ho]}"


# ---------------------------------------------------------------- geometry

def ray_in_cylinder(x, u, R, ylo, yhi):
    """Part t >= 0 of the rays x + t u inside the cylinder r_xz < R, ylo < y < yhi.

    x, u: (N, 3) relative to the target centre, u unit. Returns (hit, t0, t1, side_ends)."""
    a = u[:, 0] ** 2 + u[:, 2] ** 2
    c = x[:, 0] ** 2 + x[:, 2] ** 2 - R * R
    b = x[:, 0] * u[:, 0] + x[:, 2] * u[:, 2]
    disc = b * b - a * c
    tilted = a > EPS
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.sqrt(np.maximum(disc, 0))
        a1 = np.where(tilted, a, 1)
        r_lo = np.where(tilted, (-b - s) / a1, -np.inf)
        r_hi = np.where(tilted, (-b + s) / a1, np.inf)
        r_ok = np.where(tilted, disc > 0, c < 0)
        uy = u[:, 1]
        moving = np.abs(uy) > EPS
        uy1 = np.where(moving, uy, 1)
        ta, tb = (ylo - x[:, 1]) / uy1, (yhi - x[:, 1]) / uy1
        y_lo = np.where(moving, np.minimum(ta, tb), -np.inf)
        y_hi = np.where(moving, np.maximum(ta, tb), np.inf)
        y_ok = np.where(moving, True, (ylo < x[:, 1]) & (x[:, 1] < yhi))
    t0 = np.maximum(np.maximum(r_lo, y_lo), 0.0)
    t1 = np.minimum(r_hi, y_hi)
    return r_ok & y_ok & (t1 > t0 + EPS), t0, t1, r_hi <= y_hi


def classify(geo, v, d):
    """Classify vertices v (N, 3) [mm] with directions d (N, 3).

    Returns dict(region, in_liquid, L, exit, holder, cos_alpha); holder -1 and cos NaN for OUTSIDE."""
    g = geo
    v = np.atleast_2d(np.asarray(v, float))
    d = np.atleast_2d(np.asarray(d, float))
    n = np.linalg.norm(d, axis=1)
    u = d / np.where(n > 0, n, 1)[:, None]
    x = v - np.array([g["TargetCenterX"], g["TargetCenterY"], g["TargetCenterZ"]])
    r = np.hypot(x[:, 0], x[:, 2])
    in_vessel = (r < g["VesselRadius"]) & (np.abs(x[:, 1]) < g["VesselHalfHeight"])
    in_holder = (r < g["HolderInnerRadius"]) & (x[:, 1] > g["HolderBottomY"])
    enters, *_ = ray_in_cylinder(x, u, g["VesselRadius"], -g["VesselHalfHeight"], g["VesselHalfHeight"])
    region = np.where(in_vessel, VESSEL, np.where(in_holder, np.where(enters, GAP_CROSS, GAP), OUTSIDE))
    region = np.where(n > 0, region, OUTSIDE)
    inside = region != OUTSIDE
    in_liquid = inside & (r < g["TargetRadius"]) & (np.abs(x[:, 1]) < g["TargetHalfHeight"])

    hit, t0, t1, side = ray_in_cylinder(x, u, g["TargetRadius"], -g["TargetHalfHeight"], g["TargetHalfHeight"])
    chord = inside & (region != GAP) & hit
    L = np.where(chord, t1 - t0, 0.0)
    exit_ = np.where(chord & ~side, CAP, SIDE)

    hit, _, t1, side = ray_in_cylinder(x, u, g["HolderRadius"], g["HolderBottomY"], np.inf)
    up = ~hit | ~np.isfinite(t1)            # straight up the axis: leaves through the open top
    bottom = ~up & ~side
    h = x + np.where(up | bottom, 0.0, t1)[:, None] * u
    phi = np.degrees(np.arctan2(h[:, 0], h[:, 2]))
    hp = g["HolderWindowHalfPhi"]
    window = ((np.abs(phi) < hp) | (np.abs(phi) > 180 - hp)) & (np.abs(h[:, 1]) < g["HolderWindowHalfHeight"])
    holder = np.where(up | bottom | ~window, G10, WINDOW)
    cosa = np.where(up, 0.0, np.where(bottom, np.abs(u[:, 1]),
                                      np.clip((h[:, 0] * u[:, 0] + h[:, 2] * u[:, 2]) / g["HolderRadius"], 0, 1)))
    return dict(region=region, in_liquid=in_liquid, L=L, exit=exit_,
                holder=np.where(inside, holder, -1), cos_alpha=np.where(inside, cosa, np.nan))


# ---------------------------------------------------------------- lookup

def _axis(nodes, x):
    e = np.asarray(nodes, float)
    if len(e) == 1:
        z = np.zeros(len(x), int)
        return z, z, np.zeros(len(x))
    i = np.clip(np.searchsorted(e, x, side="right") - 1, 0, len(e) - 2)
    return i, i + 1, np.clip((x - e[i]) / (e[i + 1] - e[i]), 0, 1)


def grid_eval(gr, p, L, c):
    """Trilinear value of one grid, clamped to the axis ends, >= 0."""
    (p0, p1, wp), (l0, l1, wl), (c0, c1, wc) = _axis(gr["p"], p), _axis(gr["L"], L), _axis(gr["c"], c)
    T = gr["loss"]
    s = 0.0
    for ip, a in ((p0, 1 - wp), (p1, wp)):
        for il, b in ((l0, 1 - wl), (l1, wl)):
            for ic, w in ((c0, 1 - wc), (c1, wc)):
                s = s + a * b * w * T[ip, il, ic]
    return np.maximum(s, 0.0)


def loss(grids, trk, p, particle):
    """Loss [MeV/c] = |p|_vertex - |p|_TPC for the tracks of classify(); p is |p| at the TPC.

    0 for OUTSIDE; NaN if the particle has no grids or p <= 0."""
    p = np.asarray(p, float)
    out = np.zeros(len(p))
    for rg, ex, ho in REQUIRED:
        s = (trk["region"] == rg) & (trk["exit"] == ex) & (trk["holder"] == ho)
        if s.any():
            gr = grids.get((particle, rg, ex, ho))
            out[s] = np.nan if gr is None else grid_eval(gr, p[s], trk["L"][s], trk["cos_alpha"][s])
    out[~(p > 0)] = np.nan
    return out


class Map:
    """Map file + geometry: Map(path).classify(v, d), .loss(trk, p, particle)."""

    def __init__(self, path):
        self.path = path
        self.grids, self.geo = read_map(path)

    def classify(self, v, d):
        return classify(self.geo, v, d)

    def loss(self, trk, p, particle):
        return loss(self.grids, trk, p, particle)


# ---------------------------------------------------------------- file I/O

def write_map(path, grids, geo, header):
    """Write grids {(particle, region, exit, holder): grid} and geometry; header: comment lines (no '#')."""
    check_complete(grids)
    with open(path, "w") as f:
        for h in header:
            f.write(f"# {h}\n" if h else "#\n")
        f.write("\n")
        for k in GEO_KEYS:
            f.write(f"{k:<24}{geo[k]:<8g}# {GEO_UNIT[k]}\n")
        for key in sorted(grids):
            pa, rg, ex, ho = key
            gr = grids[key]
            f.write(f"\nparticle {pa}  # {PARTICLE_NAMES[pa]}\nregion   {rg}  # {REGION_NAMES[rg]}\n"
                    f"exit     {ex}  # {EXIT_NAMES[ex]}\nholder   {ho}  # {HOLDER_NAMES[ho]}\n")
            for ax, lab in (("p", "p [MeV/c]"), ("L", "L [mm]"), ("c", "cos")):
                f.write(f"# {lab}\n" + " ".join(f"{x:g}" for x in gr[ax]) + "\n")
            f.write(f"# p L loss(cos_0) ... loss(cos_{len(gr['c']) - 1})\n")
            for ip, pv in enumerate(gr["p"]):
                for il, lv in enumerate(gr["L"]):
                    f.write(f"{pv:g} {lv:g} " + " ".join(f"{x:.6f}" for x in gr["loss"][ip, il]) + "\n")


def check_complete(grids):
    for pa in sorted({k[0] for k in grids}):
        miss = [(pa,) + k for k in REQUIRED if (pa,) + k not in grids]
        if miss:
            raise MapError(f"{PARTICLE_NAMES[pa]} grids incomplete, missing: "
                           + ", ".join(grid_name(k) for k in miss))


def _numbers(tokens, where):
    try:
        return [float(t) for t in tokens]
    except ValueError:
        raise MapError(f"{where}: non-numeric token in {' '.join(tokens)!r}") from None


def read_map(path):
    """Strict reader with the C++ rules. Returns (grids, geo); raises MapError with the line number."""
    lines = []
    for n, raw in enumerate(open(path), 1):
        tok = raw.split("#")[0].split()
        if tok:
            lines.append((n, tok))
    geo, grids, i = {}, {}, 0

    def label(name, size):
        n, tok = lines[i]
        if tok[0] != name or len(tok) != 2:
            raise MapError(f"line {n}: expected '{name} <value>'")
        x = _numbers(tok[1:], f"line {n}")[0]
        if x != int(x) or not 0 <= x < size:
            raise MapError(f"line {n}: {name} {tok[1]} out of range 0-{size - 1}")
        return int(x)

    while i < len(lines):
        n, tok = lines[i]
        if tok[0] != "particle":
            if tok[0] not in GEO_KEYS:
                raise MapError(f"line {n}: unknown key {tok[0]!r}")
            if tok[0] in geo or len(tok) != 2:
                raise MapError(f"line {n}: duplicate key or not one value: {tok[0]}")
            geo[tok[0]] = _numbers(tok[1:], f"line {n}")[0]
            i += 1
            continue
        if i + 7 > len(lines):
            raise MapError(f"line {n}: truncated grid block")
        key = []
        for name, size in (("particle", 3), ("region", 3), ("exit", 2), ("holder", 2)):
            key.append(label(name, size))
            i += 1
        key = tuple(key)
        if key[1] == GAP and key[2] != SIDE:
            raise MapError(f"line {n}: gap grids use exit 0")
        if key in grids:
            raise MapError(f"line {n}: duplicate grid {grid_name(key)}")
        ax = {}
        for name in ("p", "L", "c"):
            m, tok = lines[i]
            ax[name] = np.array(_numbers(tok, f"line {m}"))
            if not np.all(np.isfinite(ax[name])) or np.any(np.diff(ax[name]) <= 0):
                raise MapError(f"line {m}: {name} axis not finite and strictly increasing")
            i += 1
        T = np.empty((len(ax["p"]), len(ax["L"]), len(ax["c"])))
        for ip, pv in enumerate(ax["p"]):
            for il, lv in enumerate(ax["L"]):
                if i >= len(lines):
                    raise MapError(f"{grid_name(key)}: rows end early")
                m, tok = lines[i]
                t = _numbers(tok, f"line {m}")
                if len(t) != 2 + len(ax["c"]) or t[0] != pv or t[1] != lv:
                    raise MapError(f"line {m}: expected row 'p L loss x{len(ax['c'])}' for p {pv:g} L {lv:g}")
                T[ip, il] = t[2:]
                i += 1
        if not np.all(np.isfinite(T)) or np.any(T < 0):
            raise MapError(f"{grid_name(key)}: loss not finite and >= 0")
        grids[key] = dict(p=ax["p"], L=ax["L"], c=ax["c"], loss=T)
    missing = [k for k in GEO_KEYS if k not in geo]
    if missing:
        raise MapError(f"{path}: missing keys {' '.join(missing)}")
    if not (geo["TargetRadius"] <= geo["VesselRadius"] < geo["HolderInnerRadius"] < geo["HolderRadius"]
            and geo["TargetHalfHeight"] <= geo["VesselHalfHeight"]):
        raise MapError(f"{path}: geometry not nested (target <= vessel < holder)")
    check_complete(grids)
    return grids, geo
