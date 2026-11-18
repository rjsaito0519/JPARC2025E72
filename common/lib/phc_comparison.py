"""Compare PHC widths on the first evaluation's fixed histogram population."""
from __future__ import annotations

import math
from statistics import median


#_____________________________________________________________________________
def _active(segment):
    try:
        return (math.isfinite(float(segment["width"]))
                and float(segment["width"]) >= 0
                and math.isfinite(float(segment["n_entries"]))
                and float(segment["n_entries"]) > 0)
    except (KeyError, TypeError, ValueError):
        return False


#_____________________________________________________________________________
def compare_candidates(journal):
    """Return fixed-population scores and the evidence for each candidate."""
    if not journal:
        return []
    if (journal[0].get("scoring_context") or {}).get("version") == 2:
        from . import phc_quality
        return phc_quality.compare(journal)
    reference = journal[0]
    context = reference.get("scoring_context") or {}
    detectors = context.get("required_detectors") or []
    population = {}
    baseline_reasons = []
    if context.get("version") != 1 or not detectors:
        baseline_reasons.append("initial evaluation lacks scoring context; restart/re-evaluate the series")
    for detector in detectors:
        segments = (reference.get("per_det", {}).get(detector) or {}).get("segments") or []
        population[detector] = {segment["histogram"]: segment for segment in segments if _active(segment)}
        if not population[detector]:
            baseline_reasons.append(f"{detector}: initial evaluation has no active segments")
    results = []
    for entry in journal:
        reasons = list(baseline_reasons)
        if not entry.get("ok"):
            reasons.append("evaluation incomplete")
        if entry.get("scoring_context") != reference.get("scoring_context"):
            reasons.append("particle/detector scoring context changed")
        details = []
        detector_widths = []
        for detector, expected in population.items():
            current_rows = (entry.get("per_det", {}).get(detector) or {}).get("segments") or []
            current = {segment["histogram"]: segment for segment in current_rows}
            if len(current) != len(current_rows):
                reasons.append(f"{detector}: duplicate histogram records")
            widths = []
            for name, baseline in expected.items():
                segment = current.get(name) or {}
                detail = {"detector": detector, "histogram": name,
                          "reference_entries": baseline.get("n_entries"),
                          "n_entries": segment.get("n_entries"),
                          "reference_total": baseline.get("n_total"),
                          "n_total": segment.get("n_total"),
                          "width": segment.get("width")}
                details.append(detail)
                if not _active(segment):
                    reasons.append(f"{name}: missing/empty/invalid fixed segment")
                    continue
                if not baseline.get("axes_sha") or baseline.get("axes_sha") != segment.get("axes_sha"):
                    reasons.append(f"{name}: histogram axes unknown/changed")
                widths.append(float(segment["width"]))
                try:
                    base_total = float(baseline["n_total"])
                    total = float(segment["n_total"])
                    base_entries = float(baseline["n_entries"])
                    entries = float(segment["n_entries"])
                    if not all(math.isfinite(value) for value in (base_total, total, base_entries, entries)):
                        raise ValueError("nonfinite counts")
                    if base_total < base_entries or total < entries:
                        raise ValueError("invalid counts")
                    base_flow = 1 - base_entries / base_total
                    flow = 1 - entries / total
                    if entries < base_entries or total < base_total:
                        reasons.append(f"{name}: statistics decreased")
                    if flow > base_flow:
                        reasons.append(f"{name}: out-of-range fraction increased")
                    for flow_key in ("n_tof_flow", "n_de_flow", "n_both_flow"):
                        if flow_key not in baseline:
                            continue
                        base_count = float(baseline[flow_key])
                        count = float(segment[flow_key])
                        if not math.isfinite(base_count) or not math.isfinite(count) or min(base_count, count) < 0:
                            raise ValueError("invalid axis flow")
                        detail["reference_" + flow_key] = base_count
                        detail[flow_key] = count
                        if count / total > base_count / base_total:
                            reasons.append(f"{name}: {flow_key} out-of-range fraction increased")
                    detail.update(entries_delta=entries-base_entries,
                                  reference_flow_fraction=base_flow,
                                  flow_fraction=flow, flow_fraction_delta=flow-base_flow,
                                  total_delta=total-base_total)
                except (KeyError, TypeError, ValueError, ZeroDivisionError):
                    reasons.append(f"{name}: flow/statistics unavailable")
            if widths:
                detector_widths.append(median(widths))
        fixed_width = median(detector_widths) if detector_widths else None
        results.append({"round": int(entry["round"]), "comparable": not reasons and fixed_width is not None,
                        "width": fixed_width, "reasons": reasons, "segments": details,
                        "reference_round": int(reference["round"]),
                        "missing_segments": sum(not _active(segment) for segment in details),
                        "statistics_decreased": len({reason.split(":")[0] for reason in reasons
                                                    if "statistics decreased" in reason}),
                        "flow_increased": len({reason.split(":")[0] for reason in reasons
                                              if "out-of-range fraction increased" in reason}),
                        "reference_min_entries": min((segment["reference_entries"] for segment in details),
                                                     default=None),
                        "required_segments": sum(len(segments) for segments in population.values())})
    return results
