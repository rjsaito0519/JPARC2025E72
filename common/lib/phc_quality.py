"""Fixed-population PHC quality from native TOF profiles and decoded hit checks."""
from __future__ import annotations

import math
import re
from statistics import mean

DEFAULT_POLICY = {
    "flow_fraction": .001,
    "flow_count_floor": 1.,
    "tail_fraction": .001,
    "tail_count_floor": 1.,
    "center_floor_ns": .02,
    "center_width_multiplier": 3.,
}


CORE_POLICY = {**DEFAULT_POLICY, "flow_fraction": .005, "tail_fraction": .005,
               "center_floor_ns": .05, "warning_fraction": .001,
               "minimum_width_gain": .01, "minimum_center_gain_ns": .05,
               "maximum_width_loss": .02}


def policy(values=None, *, version=2):
    defaults = CORE_POLICY if version == 3 else DEFAULT_POLICY
    result = dict(defaults if values is None else values)
    if set(result) != set(defaults) or any(
            isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or value < 0 for value in result.values()):
        raise ValueError("invalid PHC quality policy")
    return result


def active(segment):
    try:
        return (math.isfinite(float(segment["width"])) and float(segment["width"]) >= 0
                and math.isfinite(float(segment["n_entries"])) and float(segment["n_entries"]) > 0)
    except (TypeError, ValueError, KeyError):
        return False


def available(segment):
    try:
        return math.isfinite(float(segment["n_total"])) and float(segment["n_total"]) > 0
    except (TypeError, ValueError, KeyError):
        return False


def profile_metrics(counts, edges):
    """Uniform density within native bins; preserve both TOF flow bins."""
    if not isinstance(counts, list) or not isinstance(edges, list) or len(counts) != len(edges) + 1:
        raise ValueError("TOF profile shape unavailable")
    if (not all(math.isfinite(value) and value >= 0 for value in counts)
            or not all(math.isfinite(value) for value in edges)
            or any(right <= left for left, right in zip(edges, edges[1:]))):
        raise ValueError("invalid TOF profile")
    total = sum(counts)
    if total <= 0:
        return {"n_total": 0., "median": None, "half68": None, "tail_fraction": None, "tof_flow": 0.}
    def quantile(q):
        target = q * total
        if target <= counts[0] or target >= total - counts[-1]:
            return None
        cumulative = counts[0]
        for weight, lo, hi in zip(counts[1:-1], edges, edges[1:]):
            if weight > 0 and cumulative + weight >= target:
                return lo + (hi - lo) * (target - cumulative) / weight
            cumulative += weight
        return None
    q16, median, q84 = (quantile(q) for q in (.16, .5, .84))
    tails = counts[0] + counts[-1]
    for weight, lo, hi in zip(counts[1:-1], edges, edges[1:]):
        inside = max(0., min(hi, 1.5) - max(lo, -1.5)) / (hi - lo)
        tails += weight * (1. - inside)
    return {"n_total": total, "median": median,
            "half68": (q84 - q16) * .5 if q16 is not None and q84 is not None else None,
            "tail_fraction": tails / total, "tof_flow": counts[0] + counts[-1]}



def peak_profile_metrics(counts, edges, window):
    """Rank the declared peak while retaining the full population and tails."""
    full = profile_metrics(counts, edges)
    lower, upper = window
    first = max(0, sum(edge <= lower for edge in edges) - 1)
    last = min(len(edges) - 2, sum(edge <= upper for edge in edges) - 1)
    selected = [0.] + [weight if first <= index <= last else 0.
                       for index, weight in enumerate(counts[1:-1])] + [0.]
    core = profile_metrics(selected, edges)
    full.update(median=core["median"], half68=core["half68"],
                peak_entries=core["n_total"], peak_window=list(window),
                peak_effective_edges=[edges[first], edges[last + 1]])
    return full


def aggregate(detector, names=None):
    rows = detector.get("segments") or []
    if names is not None:
        lookup = {row["histogram"]: row for row in rows}
        rows = [lookup[name] for name in names if name in lookup]
    rows = [row for row in rows if re.search(r"seg\d+U_", row["histogram"])]
    edges = detector.get("tof_edges")
    if not rows or not edges:
        raise ValueError("fixed U-side TOF profiles unavailable")
    counts = [0.] * (len(edges) + 1)
    for row in rows:
        values = row.get("tof_counts")
        profile_metrics(values, edges)
        if not math.isclose(sum(values), float(row["n_total"]), rel_tol=0., abs_tol=1.e-8):
            raise ValueError("TOF profile counts disagree with fill total")
        for index, value in enumerate(values):
            counts[index] += value
    window = detector.get("peak_window")
    return peak_profile_metrics(counts, edges, window) if window else profile_metrics(counts, edges)


def integrity_from_score(score):
    result = {"scoring_context": score.get("scoring_context"),
            "per_det": {det: {"time_validity": item.get("time_validity"),
                             "tof_edges": item.get("tof_edges"),
                             "segments": [{"histogram": row["histogram"], "n_total": row.get("n_total"),
                                           "adc_flow": row.get("n_de_flow", 0.) + row.get("n_both_flow", 0.),
                                           "axes_sha": row.get("axes_sha"),
                                           "tof_counts": row.get("tof_counts")}
                                          for row in item.get("segments", [])]}
                        for det, item in score.get("per_det", {}).items()}}

    for detector, item in score.get("per_det", {}).items():
        if item.get("peak_window"):
            result["per_det"][detector]["peak_window"] = list(item["peak_window"])
    return result


def compare(journal):
    if not journal:
        return []
    first = journal[0]
    context = first.get("scoring_context") or {}
    version = context.get("version")
    settings = policy(context.get("quality_policy"), version=version)
    detectors = context.get("required_detectors") or []
    primary = context.get("primary_detectors") or []
    population = {det: {row["histogram"]: row for row in first.get("per_det", {}).get(det, {}).get("segments", [])
                        if (available(row) if version == 3 else active(row))} for det in detectors}
    initial = first.get("input_integrity")
    common_reasons = []
    if (context.get("beam_selection") or {}).get("ok") is False:
        common_reasons.append("declared particle/BAC selection requirement is not met")
    if version == 3 and not all((context.get("fixed_parameters") or {}).get(name) for name in ("HDPRM", "USER")):
        common_reasons.append("fixed ADC/TDC/USER parameters unavailable")
    if not detectors or not primary or any(det not in detectors for det in primary):
        common_reasons.append("primary detector definition unavailable")
    if not initial or initial.get("scoring_context") != context:
        common_reasons.append("initial decode integrity unavailable; restart/re-evaluate the series")
    for det in detectors:
        if not population[det] and det in primary:
            common_reasons.append(f"{det}: initial evaluation has no active segments")
    baseline, input_baseline = {}, {}
    for det in detectors:
        if not population[det] and det not in primary:
            continue
        try:
            baseline[det] = aggregate(first["per_det"][det], population[det])
            input_baseline[det] = aggregate(initial["per_det"][det], population[det])
        except (ValueError, TypeError, KeyError) as error:
            common_reasons.append(f"{det}: {error}")
    results = []
    for entry in journal:
        reasons = list(common_reasons)
        if not entry.get("ok"):
            reasons.append("evaluation incomplete")
        if entry.get("scoring_context") != context:
            reasons.append("particle/detector/scoring policy changed")
        details, quality, allowed, warnings = [], {}, [], []
        widths = []
        totals_changed = set()
        blocked_flows = set()
        for det in detectors:
            item = entry.get("per_det", {}).get(det, {})
            current_rows = item.get("segments") or []
            current = {row["histogram"]: row for row in current_rows}
            if len(current) != len(current_rows):
                reasons.append(f"{det}: duplicate histogram records")
            if item.get("peak_window") != first.get("per_det", {}).get(det, {}).get("peak_window"):
                reasons.append(f"{det}: declared peak window changed")
            if item.get("tof_edges") != first.get("per_det", {}).get(det, {}).get("tof_edges"):
                reasons.append(f"{det}: TOF profile edges changed")
            original = (initial or {}).get("per_det", {}).get(det, {})
            valid0, valid1 = original.get("time_validity"), item.get("time_validity")
            try:
                if not isinstance(valid0, dict) or not isinstance(valid1, dict):
                    raise ValueError
                if (not valid0.get("raw_hit_signature")
                        or valid0.get("raw_hit_signature") != valid1.get("raw_hit_signature")):
                    reasons.append(f"{det}: raw hit population or uncorrected timing changed")
                if valid0["n_cmt"] != valid1["n_cmt"] or valid0["n_events"] != valid1["n_events"]:
                    reasons.append(f"{det}: decoded hit population changed")
                if not set(map(tuple, valid1["nonfinite_keys"])).issubset(set(map(tuple, valid0["nonfinite_keys"]))):
                    reasons.append(f"{det}: newly non-finite corrected time")
            except (KeyError, TypeError, ValueError):
                reasons.append(f"{det}: corrected-time validity unavailable")
            for row in original.get("segments", []):
                name = row["histogram"]
                value = current.get(name, {})
                if value.get("n_total") != row.get("n_total"):
                    reasons.append(f"{name}: fill total changed from initial decode")
                    totals_changed.add(name)
                if not row.get("axes_sha") or value.get("axes_sha") != row.get("axes_sha"):
                    reasons.append(f"{name}: axes changed from initial decode")
                if value and value.get("n_de_flow", 0.) + value.get("n_both_flow", 0.) != row.get("adc_flow"):
                    reasons.append(f"{name}: ADC-domain population changed")
            for name, previous in population[det].items():
                value = current.get(name) or {}
                detail = {"detector": det, "histogram": name, "reference_entries": previous.get("n_entries"),
                          "n_entries": value.get("n_entries"), "reference_total": previous.get("n_total"),
                          "n_total": value.get("n_total"), "width": value.get("width")}
                if value.get("n_total") and previous.get("n_total"):
                    detail.update(reference_flow_fraction=1. - previous["n_entries"] / previous["n_total"],
                                  flow_fraction=1. - value["n_entries"] / value["n_total"],
                                  entries_delta=value["n_entries"] - previous["n_entries"],
                                  total_delta=value["n_total"] - previous["n_total"])
                for key in ("n_tof_flow", "n_de_flow", "n_both_flow"):
                    detail["reference_" + key] = previous.get(key)
                    detail[key] = value.get(key)
                details.append(detail)
                if not (available(value) if version == 3 else active(value)):
                    reasons.append(f"{name}: missing/empty/invalid fixed segment")
                    continue
                if not previous.get("axes_sha") or value.get("axes_sha") != previous.get("axes_sha"):
                    reasons.append(f"{name}: histogram axes changed")
                if value.get("n_total") != previous.get("n_total"):
                    reasons.append(f"{name}: fill total changed")
                    totals_changed.add(name)
                profile = value.get("tof_counts")
                if (not isinstance(profile, list) or len(profile) < 3
                        or not math.isclose(sum(profile), float(value["n_total"]), rel_tol=0., abs_tol=1.e-8)
                        or not math.isclose(sum(profile[1:-1]), value["n_entries"] + value.get("n_de_flow", 0.), rel_tol=0., abs_tol=1.e-8)
                        or not math.isclose(profile[0] + profile[-1], value.get("n_tof_flow", 0.) + value.get("n_both_flow", 0.), rel_tol=0., abs_tol=1.e-8)):
                    reasons.append(f"{name}: TOF profile/fill counters inconsistent")
                before_y = previous.get("n_tof_flow", 0.) + previous.get("n_both_flow", 0.)
                after_y = value.get("n_tof_flow", 0.) + value.get("n_both_flow", 0.)
                delta = after_y - before_y
                tolerance = max(settings["flow_count_floor"], settings["flow_fraction"] * previous["n_total"])
                if version == 3:
                    tolerance = max(tolerance, settings["flow_fraction"] *
                                    baseline.get(det, {}).get("n_total", previous["n_total"]))
                if delta > tolerance + 1.e-9:
                    reasons.append(f"{name}: TOF flow increase exceeds allowance")
                    blocked_flows.add(name)
                elif delta > 0:
                    allowed.append({"histogram": name, "delta": delta, "allowance": tolerance})
                    if version == 3 and delta > max(settings["flow_count_floor"], settings["warning_fraction"] * previous["n_total"]):
                        warnings.append(f"{name}: small TOF flow increase; raw fills unchanged")
            if not population[det] and det not in primary:
                quality[det] = {"status": "no_evaluable_population_in_first_round"}
                continue
            try:
                current_global = aggregate(item, population[det])
                base = baseline[det]
                allowances = {}
                for label, reference in (("first updated round", base),
                                         ("before first update", input_baseline[det])):
                    total = reference["n_total"]
                    flow_limit = max(settings["flow_count_floor"], settings["flow_fraction"] * total)
                    if "peak_entries" in reference:
                        peak_loss = reference["peak_entries"] - current_global.get("peak_entries", 0.)
                        if peak_loss > flow_limit + 1.e-9:
                            reasons.append(f"{det}: main-peak population loss exceeds allowance ({label})")
                        elif version == 3 and peak_loss > max(settings["flow_count_floor"], settings["warning_fraction"] * total):
                            warnings.append(f"{det}: small main-peak population decrease ({label})")
                    if current_global["tof_flow"] - reference["tof_flow"] > flow_limit + 1.e-9:
                        reasons.append(f"{det}: aggregate TOF flow increase exceeds allowance ({label})")
                    tail_limit = max(settings["tail_fraction"], settings["tail_count_floor"] / total) if total else 1.
                    if (reference["tail_fraction"] is not None and
                            current_global["tail_fraction"] - reference["tail_fraction"] > tail_limit + 1.e-12):
                        reasons.append(f"{det}: global TOF tail fraction worsened ({label})")
                    if version == 3:
                        if current_global["tof_flow"] - reference["tof_flow"] > max(settings["flow_count_floor"], settings["warning_fraction"] * total):
                            warnings.append(f"{det}: small TOF flow increase ({label})")
                        if reference["tail_fraction"] is not None and current_global["tail_fraction"] - reference["tail_fraction"] > max(settings["warning_fraction"], settings["tail_count_floor"] / total if total else 1.):
                            warnings.append(f"{det}: small tail increase ({label})")
                    center_limit = settings["center_floor_ns"]
                    if reference["median"] is not None:
                        if current_global["median"] is None:
                            reasons.append(f"{det}: global TOF center unavailable ({label})")
                        else:
                            count = max(1, min(total, valid0["n_events"], valid1["n_events"]))
                            width = max(reference["half68"] or 0., current_global["half68"] or 0.)
                            center_limit = max(center_limit, settings["center_width_multiplier"] * width / math.sqrt(count))
                            if abs(current_global["median"]) > abs(reference["median"]) + center_limit:
                                reasons.append(f"{det}: global TOF center moved away from zero ({label})")
                    allowances[label] = {"flow_count_allowance": flow_limit,
                                         "tail_fraction_allowance": tail_limit,
                                         "center_allowance_ns": center_limit}
                quality[det] = {"reference": base, "input": input_baseline[det],
                                "current": current_global, "allowances": allowances}
                if det in primary:
                    if current_global["half68"] is None:
                        reasons.append(f"{det}: primary TOF width unavailable")
                    else:
                        widths.append(current_global["half68"])
            except (ValueError, KeyError, TypeError) as error:
                reasons.append(f"{det}: TOF profile quality unavailable ({error})")
        width = mean(widths) if len(widths) == len(primary) and primary else None
        reasons = list(dict.fromkeys(reasons))
        results.append({"round": int(entry["round"]), "comparable": not reasons and width is not None,
                        "width": width, "reasons": reasons, "warnings": list(dict.fromkeys(warnings)),
                        "segments": details, "quality": quality,
                        "allowed_flow_changes": allowed, "reference_round": int(first["round"]),
                        "missing_segments": sum(not (available(row) if version == 3 else active(row)) for row in details),
                        "statistics_decreased": len(totals_changed), "flow_increased": len(blocked_flows),
                        "reference_min_entries": min((row["reference_entries"] for row in details), default=None),
                        "required_segments": sum(map(len, population.values()))})
    return results

def summary_lines(comparison):
    lines = []
    def number(value):
        return f"{value:+.3f}" if value is not None and math.isfinite(value) else "不明"
    for detector, details in comparison.get("quality", {}).items():
        current = details.get("current")
        if not current:
            lines.append(f"{detector}: 評価に使えるデータなし")
            continue
        width = f"{current['half68']:.3f}" if current.get("half68") is not None else "不明"
        tail = f"{100 * current['tail_fraction']:.2f}%" if current.get("tail_fraction") is not None else "不明"
        lines.append(f"{detector}: 中心 {number(current.get('median'))} ns / 幅 {width} ns / 裾 {tail}")
        if "peak_entries" in current:
            lines.append(f"  主ピーク内 {current['peak_entries']:.0f} / 全体 {current['n_total']:.0f} hit（All）")
    return lines


def select_candidate(comparisons, context):
    """Prefer the initial correction when a refit offers no material peak gain."""
    eligible = [row for row in comparisons if row["comparable"]]
    if not eligible:
        return None
    if context.get("version") != 3:
        return min(eligible, key=lambda row: (row["width"], row["round"]))
    initial = next((row for row in eligible if row["round"] == 0), None)
    if initial is None:
        return min(eligible, key=lambda row: (row["width"], row["round"]))
    settings = policy(context.get("quality_policy"), version=3)
    candidates = []
    primary = context["primary_detectors"]
    for row in eligible:
        if row["round"] == 0 or row["width"] > initial["width"] * (1 + settings["maximum_width_loss"]):
            continue
        center_gain = mean(abs(initial["quality"][det]["current"]["median"]) -
                           abs(row["quality"][det]["current"]["median"]) for det in primary)
        if (row["width"] <= initial["width"] * (1 - settings["minimum_width_gain"])
                or center_gain >= settings["minimum_center_gain_ns"]):
            candidates.append(row)
    return min(candidates, key=lambda row: (row["width"], row["round"])) if candidates else initial


def stable_input(score):
    """Keep only quantities that a PHC correction is not allowed to change."""
    audit = integrity_from_score(score)
    for detector in audit["per_det"].values():
        validity = detector.get("time_validity") or {}
        detector["time_validity"] = {key: validity.get(key)
                                     for key in ("raw_hit_signature", "n_cmt", "n_events")}
        for row in detector["segments"]:
            row.pop("tof_counts", None)
        detector["segments"].sort(key=lambda row: row["histogram"])
    return audit
