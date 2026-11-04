#!/usr/bin/env python3
"""Build refs.bib from Crossref, arXiv and DataCite (Zenodo) records saved in tools/records/.

Every title, author list, year, venue, volume, issue, pages, DOI and arXiv
identifier is copied from a saved record.  If a record is missing it is fetched
once (api.crossref.org/works/DOI or export.arxiv.org/api/query?id_list=ID) and
the raw response is saved before use; existing records are never overwritten.
Only mechanical transformations are applied:
  * HTML entities are unescaped, inline HTML tags removed, and non-ASCII
    letters converted to LaTeX accents;
  * family (given) names entirely in capitals are converted to title case;
  * an author record with no given name is set as one braced name;
  * titles are wrapped in braces so that their capitalisation is kept;
  * page ranges use an en dash;
  * year = year of print publication if present, otherwise of first issue
    (online, posting, or as a last resort record creation, in that order).
The Evidence Press release is copied from its DataCite record (title, creators, year, DOI).
The script stops with an error
if a record is missing and cannot be fetched, or a character cannot be converted.

Usage: python3 -I tools/make_refs.py refs.bib tools/records
"""
import html
import json
import os
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET

CROSSREF = {
    "10.1109/TNNLS.2025.3551275": "Li2025",
    "10.1215/ijm/1446819297": "Lyons2014",
    "10.1016/j.jfa.2025.110975": "PS2025",
    "10.1007/s10959-008-0151-0": "BomanLindskog2009",
    "10.1007/978-1-4419-6055-9": "Helgason2011",
    "10.1017/CBO9781107359949.008": "McDiarmid1989",
    "10.1007/b13794": "Tsybakov2009",
}
ARXIV = {
    "2211.10066": "Bonet2023",
    "2403.06560": "BDC2025",
    "2106.03306": "Chami2021",
    "2003.05783": "Nadjahi2020",
    "2311.15874": "KT2026",
}
SOFTWARE = ""
DATACITE = {"10.5281/zenodo.23295188": "EPnoninteger2026"}

ACCENTS = {
    "č": r"{\v{c}}", "Č": r"{\v{C}}", "ş": r"{\c{s}}", "Ş": r"{\c{S}}",
    "á": r"{\'a}", "Á": r"{\'A}", "é": r"{\'e}", "É": r"{\'E}", "í": r"{\'i}",
    "ó": r"{\'o}", "ú": r"{\'u}", "ñ": r"{\~n}", "Ñ": r"{\~N}", "ï": r'{\"i}',
    "ü": r'{\"u}', "ö": r'{\"o}', "ä": r'{\"a}', "è": r"{\`e}", "ç": r"{\c c}",
    "Ö": r'{\"O}', "Ü": r'{\"U}', "ß": r"{\ss}",
    "‐": "-", "‑": "-", "–": "--", "—": "---",
    "“": "``", "”": "''", "‘": "`", "’": "'", " ": "~",
}
UA = {"User-Agent": "make_refs (mailto:anonymous@example.org)"}


def tex(s):
    s = re.sub(r"<[^>]+>", "", html.unescape(s))
    s = " ".join(s.split())
    out = []
    for ch in s:
        if ch in ACCENTS:
            out.append(ACCENTS[ch])
        elif ord(ch) < 128:
            out.append({"&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_"}.get(ch, ch))
        else:
            raise SystemExit(f"unconvertible character {ch!r} in {s!r}")
    return "".join(out)


def family(name):
    name = html.unescape(name)
    if name.isupper():
        name = name.capitalize()
    return tex(name)


def given(name):
    name = html.unescape(name)
    if name.isupper():
        name = name.title()
    return tex(name)


def fname(ident):
    return "".join(c if c.isalnum() or c in ".-" else "_" for c in ident)


def fetch(url, path):
    if os.path.exists(path):
        return
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read()
    with open(path, "xb") as fh:
        fh.write(data)


def year(m):
    for k in ("published-print", "issued", "published-online", "posted", "created"):
        dp = m.get(k, {}).get("date-parts")
        if dp and dp[0] and dp[0][0]:
            return dp[0][0]
    raise SystemExit("no year")


def fmt(kind, key, f):
    body = ",\n".join(f"  {k} = {v}" if v.startswith("{{") else f"  {k} = {{{v}}}"
                      for k, v in f.items())
    return f"@{kind}{{{key},\n{body}\n}}\n"


def crossref_entry(doi, key, m):
    authors = " and ".join(
        f"{family(a['family'])}, {given(a['given'])}" if a.get("given")
        else "{" + family(a["family"]) + "}" for a in m["author"])
    f = {"author": authors, "title": "{{" + tex(m["title"][0]) + "}}",
         "year": str(year(m))}
    typ = m["type"]
    container = tex(m["container-title"][0]) if m.get("container-title") else ""
    if typ == "journal-article":
        kind = "article"
        f["journal"] = container
        if m.get("volume"):
            f["volume"] = m["volume"]
        if m.get("issue"):
            f["number"] = m["issue"]
        if m.get("page"):
            f["pages"] = m["page"].replace("-", "--")
        elif m.get("article-number"):
            f["pages"] = m["article-number"]
    elif typ in ("proceedings-article", "book-chapter"):
        kind = "inproceedings" if typ == "proceedings-article" else "incollection"
        f["booktitle"] = container
        if m.get("volume"):
            f["volume"] = m["volume"]
        if m.get("page"):
            f["pages"] = m["page"].replace("-", "--")
        if m.get("publisher"):
            f["publisher"] = tex(m["publisher"])
    elif typ in ("book", "monograph"):
        kind = "book"
        if container:
            f["series"] = container
        if m.get("publisher"):
            f["publisher"] = tex(m["publisher"])
    else:
        raise SystemExit(f"unhandled type {typ}")
    f["doi"] = m["DOI"]
    return fmt(kind, key, f)


def arxiv_entry(aid, key, raw):
    ns = {"a": "http://www.w3.org/2005/Atom", "x": "http://arxiv.org/schemas/atom"}
    e = ET.fromstring(raw).find("a:entry", ns)
    if e is None or aid not in e.find("a:id", ns).text:
        raise SystemExit(f"arXiv record mismatch for {aid}")
    names = [a.find("a:name", ns).text.strip() for a in e.findall("a:author", ns)]

    def bibname(n):
        parts = n.split()
        return f"{tex(parts[-1])}, {tex(' '.join(parts[:-1]))}"
    f = {"author": " and ".join(bibname(n) for n in names),
         "title": "{{" + tex(e.find("a:title", ns).text) + "}}",
         "year": e.find("a:published", ns).text[:4],
         "eprint": aid, "archivePrefix": "arXiv",
         "howpublished": f"arXiv:{aid}"}
    c = e.find("x:comment", ns)
    if c is not None and c.text:
        f["note"] = "arXiv comment: " + tex(c.text)
    return fmt("misc", key, f)


def main():
    out, rec = sys.argv[1], sys.argv[2]
    os.makedirs(rec, exist_ok=True)
    parts = ["% Generated by tools/make_refs.py from saved Crossref/arXiv/DataCite records. Do not edit by hand.\n"]
    for doi, key in CROSSREF.items():
        path = os.path.join(rec, "crossref_" + fname(doi) + ".json")
        fetch("https://api.crossref.org/works/" + doi, path)
        with open(path) as fh:
            m = json.load(fh)["message"]
        if m["DOI"].lower() != doi.lower():
            raise SystemExit(f"DOI mismatch {doi} vs {m['DOI']}")
        parts.append(crossref_entry(doi, key, m))
    for aid, key in ARXIV.items():
        path = os.path.join(rec, "arxiv_" + fname(aid) + ".xml")
        fetch("https://export.arxiv.org/api/query?id_list=" + aid, path)
        with open(path, "rb") as fh:
            parts.append(arxiv_entry(aid, key, fh.read()))
    for doi, key in DATACITE.items():
        path = os.path.join(rec, "datacite_" + fname(doi) + ".json")
        fetch("https://api.datacite.org/dois/" + doi, path)
        with open(path) as fh:
            at = json.load(fh)["data"]["attributes"]
        if at["doi"].lower() != doi.lower():
            raise SystemExit(f"DOI mismatch {doi}")
        names = " and ".join("{" + tex(c["name"]) + "}" for c in at["creators"])
        parts.append(fmt("misc", key, dict([("author", names), ("title", "{" + tex(at["titles"][0]["title"]) + "}"),
                                      ("year", str(at["publicationYear"])), ("publisher", tex(at["publisher"] if isinstance(at["publisher"], str) else at["publisher"]["name"])),
                                      ("note", "Version " + tex(at.get("version") or "") + ", Evidence Press research release"),
                                      ("doi", at["doi"])])))
    with open(out, "w") as fh:
        fh.write("\n".join(parts))
    print(f"wrote {len(CROSSREF) + len(ARXIV) + len(DATACITE)} entries to {out}")


if __name__ == "__main__":
    main()
