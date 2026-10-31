"""Core of the methods-number registry: draft normalization, phrase/token matching, formatting, row assembly."""
import re, time
import numpy as np

WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}
SUPERSCRIPT_MINUS = "−"


def normalize_draft(raw):
    """Text as extracted from the PDF/docx -> one line: drop page-number lines, collapse whitespace, and write superscript
    exponents the way the document's plain text reads (10−7: U+2212, no caret, no space; punctuation closed up after it)."""
    lines = [l for l in raw.replace("\r", "\n").split("\n") if not re.fullmatch(r"\s*\d{1,2}\s*", l)]
    t = re.sub(r"\s+", " ", " ".join(lines))
    t = re.sub(r"10 ([−-]\d+)", r"10\1", t)
    t = re.sub(r"(10[−-]\d+) ([;.,])", r"\1\2", t)
    return t.replace(" )", ")").replace("( ", "(")


def token_count(s, tok, numeric):
    if not numeric:
        return s.count(tok), [m.start() for m in re.finditer(re.escape(tok), s)]
    pat = re.compile(r"(?<![\d.,^−])" + re.escape(tok) + r"(?![\d]|[.,]\d)")
    ms = list(pat.finditer(s))
    return len(ms), [m.start() for m in ms]


def fmt(v, how):
    """Format a recomputed value the way the draft writes it."""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return None
    kind, _, arg = how.partition(":")
    if kind == "int": return f"{int(round(v)):,}"
    if kind == "intplain": return f"{int(round(v))}"
    if kind == "word": return WORDS.get(int(v), str(v))
    if kind == "Word":  # a count that opens a sentence ("Three clusters")
        return WORDS.get(int(v), str(v)).capitalize()
    if kind == "dec": return f"{float(v):.{int(arg)}f}"
    if kind == "pct": return f"{100 * float(v):.{int(arg)}f}%"            # v is a fraction
    if kind == "pctnum": return f"{float(v):.{int(arg)}f}%"              # v is already a percentage
    if kind == "pctnum_plain": return f"{float(v):.{int(arg)}f}"
    # Exponents are written as the document's plain text shows superscripts: 10−7 (U+2212, no caret). Only negative
    # exponents are supported; a positive one (10^5 -> "105") would be ambiguous without the caret.
    if kind == "sci":
        m, e = f"{float(v):.{int(arg)}e}".split("e")
        assert int(e) < 0, f"positive exponent {e} cannot be written without a caret"
        return f"{m} × 10{SUPERSCRIPT_MINUS}{abs(int(e))}"
    if kind == "log10sci":                                               # v = log10(x); arg = significant figures
        m, e = mant_exp(float(v), int(arg))
        assert e < 0, f"positive exponent {e} cannot be written without a caret"
        return f"{m:.{int(arg) - 1}f} × 10{SUPERSCRIPT_MINUS}{abs(e)}"
    if kind == "pow10":                                                  # 1e-7 -> 10−7
        e = int(round(np.log10(float(v))))
        assert e < 0, f"positive exponent {e} cannot be written without a caret"
        return f"10{SUPERSCRIPT_MINUS}{abs(e)}"
    if kind == "str": return str(v)
    raise ValueError(how)


def mant_exp(log10x, sig):
    """log10(x) -> (mantissa rounded to `sig` significant figures, exponent), without forming x (no underflow)."""
    e = int(np.floor(log10x)); m = round(10 ** (log10x - e), sig - 1)
    if m >= 10: m, e = round(m / 10, sig - 1), e + 1
    return m, e


def raw_value(v, how):
    """value_raw: the unrounded value. Small floats keep 6 significant figures (rounding to 6 decimals would zero them);
    log10-valued rows give the value itself in e-notation, computed from the log."""
    if v is None: return ""
    if how.startswith("log10sci"):
        m, e = mant_exp(float(v), 6)
        return f"{m:.5f}e{e}"
    if isinstance(v, (float, np.floating)):
        v = float(v)
        return format(v, ".6g") if 0 < abs(v) < 1e-4 else repr(round(v, 6))
    return str(v)


class Row:
    __slots__ = ("key", "phrase", "shown", "unit", "kind", "step", "desc", "how", "nth", "fn")

    def __init__(self, key, phrase, shown, how, fn, unit="", kind="count", step="", desc="", nth=0):
        self.key, self.phrase, self.shown, self.how, self.fn = key, phrase, shown, how, fn
        self.unit, self.kind, self.step, self.desc, self.nth = unit, kind, step, desc, nth


def evaluate(rows, ctx, draft):
    """rows: list of Row; ctx: object passed to every fn; draft: normalized draft text. fn(ctx) -> (value, source_file[, note])."""
    out, now = [], time.strftime("%Y-%m-%d %H:%M")
    seen = set()
    for r in rows:
        assert r.key not in seen, f"duplicate key {r.key}"; seen.add(r.key)
        numeric = r.how.split(":")[0] not in ("str",)
        rec = dict(key=r.key, value=None, unit=r.unit, description=r.desc, source_file="", source_step=r.step, doc_phrase=r.phrase, updated=now,
                   draft_value=r.shown, display_value=None, status="", kind=r.kind, note="")
        nphr = draft.count(r.phrase)
        ntok, _ = token_count(r.phrase, r.shown, numeric)
        try:
            res = r.fn(ctx)
            val, src = res[0], res[1]
            note = res[2] if len(res) > 2 else ""
            rec.update(value=val, source_file=src, note=note)
            rec["display_value"] = fmt(val, r.how) if val is not None else None
            rec["value_raw"] = raw_value(val, r.how)
        except Exception as e:
            rec.update(note=f"not recomputed: {type(e).__name__}: {e}"[:300])
        if nphr != 1:
            rec["status"] = f"phrase_found_{nphr}x"
        elif ntok <= r.nth:
            rec["status"] = "value_not_in_phrase"
        elif rec["display_value"] is None:
            rec["status"] = "no_source"
        else:
            rec["status"] = "match" if rec["display_value"] == r.shown else "MISMATCH"
        out.append(rec)
    return out


def apply_to_text(draft, rows_out):
    """Return the draft with every MISMATCH replaced by its display value (token-exact, inside its phrase). Not used to edit any file."""
    t = draft
    for r in rows_out:
        if r["status"] != "MISMATCH": continue
        i = t.find(r["doc_phrase"])
        if i < 0: continue
        numeric = r.get("_numeric", True)
        n, pos = token_count(r["doc_phrase"], r["draft_value"], numeric)
        k = pos[0] if pos else None
        if k is None: continue
        new_phrase = r["doc_phrase"][:k] + r["display_value"] + r["doc_phrase"][k + len(r["draft_value"]):]
        t = t[:i] + new_phrase + t[i + len(r["doc_phrase"]):]
    return t
