"""HTML report of a WHISCY prediction, with a Mol* view (MolViewSpec).

The report refers to the files by relative names, so the directory can be
served as is, e.g. with `python -m http.server -d <directory>`.
"""

import os
import shutil
import warnings
from html import escape
from string import Template

with warnings.catch_warnings():
    # molviewspec 1.8 passes a deprecated keyword to pydantic's Field.
    warnings.filterwarnings("ignore", "Using extra keyword arguments on `Field`")
    from molviewspec import (
        MVSJ,
        ComponentExpression,
        GlobalMetadata,
        States,
        create_builder,
    )

MOLSTAR_VERSION = "5.13.1"
# Ends of the score scale, light enough for black text on them.
BLUE, RED = (67, 147, 195), (214, 96, 77)

TEMPLATE = Template(
    """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>WHISCY Report</title>
<script src="https://cdn.jsdelivr.net/npm/molstar@${molstar}/build/mvs-stories/mvs-stories.js"></script>
<link rel="stylesheet" type="text/css"
  href="https://cdn.jsdelivr.net/npm/molstar@${molstar}/build/mvs-stories/mvs-stories.css">
<style>
  body { margin: 0; font-family: sans-serif; }
  header, section { padding: 0 1rem; }
  .container { display: flex; height: 600px; }
  #viewer { flex: 2 1 0; min-width: 0; position: relative; }
  #controls { flex: 1 1 0; min-width: 300px; padding: 1rem;
              border-left: 1px solid #ccc; background: #f6f5f3; }
  table { border-collapse: collapse; margin: 0.5rem 0 1rem; }
  th, td { border: 1px solid #a0a0a0; padding: 4px 8px; text-align: right; }
  tr.active { font-weight: bold; }
  pre { white-space: pre-wrap; }
  table.sortable th { cursor: pointer; }
  #file-hint { background: #fff3cd; padding: 0.5rem; }
  @media (max-width: 900px) {
    .container { flex-direction: column; height: auto; }
    #viewer { height: 50vh; }
    #controls { border-left: none; border-top: 1px solid #ccc; }
  }
</style>
</head>
<body>
<header>
<h2>WHISCY Report</h2>
<details><summary>Options</summary><ul>${options}</ul></details>
<p>Download: ${downloads}</p>
<p id="file-hint" hidden>The browser does not let a page opened from disk load
the files of the viewer. Serve this directory with
<code>python -m http.server -d &lt;directory&gt;</code> and open
<a href="http://localhost:8000/report.html">http://localhost:8000/report.html</a>.</p>
</header>
<div class="container">
  <div id="viewer"><mvs-stories-viewer /></div>
  <div id="controls"><mvs-stories-snapshot-markdown /></div>
</div>
<section>
<h3>HADDOCK restraints</h3>
<p><b>Active residues</b> (high WHISCY score), as in active.txt:</p>
<pre>${active}</pre>
<p><b>Passive residues</b> (surface residues next to an active residue), as in
passive.txt:</p>
<pre>${passive}</pre>
<p>You should not trust WHISCY predictions blindly, but always try to relate
them to information that you already have. WHISCY works very well for proteins
with a single conserved protein binding site. It should not be used for
non-conserved interactions, such as antibody-antigens. If you know that your
protein has multiple binding sites, you should verify if a predicted site
could be the one you are interested in.</p>
<h3>WHISCY scores</h3>
<p><i>Click a column header to sort the table by that column; click again to
reverse the order.</i></p>
${table}</section>
<script>
// Sort a table by a column on a click on its header, again to reverse.
document.querySelectorAll('table.sortable th').forEach(th => {
  th.title = 'Click to sort';
  th.addEventListener('click', () => {
    const table = th.closest('table');
    const col = th.cellIndex;
    const rows = [...table.rows].slice(1);
    const key = r => {
      const text = r.cells[col].textContent;
      const x = parseFloat(text);
      return isNaN(x) ? text : x;
    };
    const dir = th.dataset.dir === 'desc' ? 1 : -1;
    th.dataset.dir = dir === 1 ? 'asc' : 'desc';
    rows.sort((a, b) => {
      const x = key(a), y = key(b);
      return dir * (x < y ? -1 : x > y ? 1 : 0);
    });
    rows.forEach(r => table.tBodies[0].appendChild(r));
  });
});
if (location.protocol === 'file:') {
  document.getElementById('file-hint').hidden = false;
} else {
  mvsStories.loadFromURL('state.mvsj', { format: 'mvsj' });
}
</script>
</body>
</html>
"""  # noqa: E501
)


def _shade(t, rgb):
    """Return the hex color at fraction t from white to rgb."""
    return "#{:02x}{:02x}{:02x}".format(*[round(255 + t * (c - 255)) for c in rgb])


def _score_color(score, smin, smax):
    """Blue (lowest score) to white (0) to red (highest score)."""
    if score >= 0:
        return _shade(score / smax if smax > 0 else 0, RED)
    return _shade(score / smin, BLUE)


def _snapshot(pdb, key, title, description, colors):
    """Return a Mol* snapshot of the structure, colored by residue number."""
    builder = create_builder()
    surface = (
        builder.download(url=pdb)
        .parse(format="pdb")
        .model_structure()
        .component()
        .focus()
        .representation(type="surface")
        .color(color="lightgrey")
    )
    for resi, color in colors:
        surface.color(color=color, selector=ComponentExpression(auth_seq_id=resi))
    return builder.get_snapshot(
        key=key,
        title=title,
        description=description,
        description_format="markdown",
    )


def generate_report(directory, pdb, files, scores, active, passive, options):
    """Write report.html and state.mvsj in directory.

    pdb is the input structure (copied to directory, so the report can load
    it), files the other files to link to, scores a dict of WHISCY scores per
    residue (e.g. {"V75": 0.16}), active and passive the HADDOCK residue ids
    and options the options shown in the report. Returns the path of
    report.html.
    """
    pdb_name = os.path.basename(pdb)
    target = os.path.join(directory, pdb_name)
    if not (os.path.exists(target) and os.path.samefile(pdb, target)):
        shutil.copyfile(pdb, target)

    smin, smax = min(scores.values()), max(scores.values())
    colors = [
        (int(res[1:]), _score_color(score, smin, smax)) for res, score in scores.items()
    ]
    haddock = [(int(res[1:]), "red") for res in active]
    haddock += [(int(res[1:]), "green") for res in passive]
    snapshots = [
        _snapshot(
            pdb_name,
            "scores",
            "WHISCY scores",
            "Surface residues colored by their WHISCY score, from blue "
            "({:.3f}) to white (0) to red ({:.3f}); residues without a score "
            "in grey.".format(smin, smax),
            colors,
        ),
        _snapshot(
            pdb_name,
            "haddock",
            "HADDOCK residues",
            "Active residues ({:d}) in red, passive residues ({:d}) in "
            "green.".format(len(active), len(passive)),
            haddock,
        ),
    ]
    state = MVSJ(
        data=States(
            snapshots=snapshots,
            metadata=GlobalMetadata(description="WHISCY prediction"),
        )
    )
    with open(os.path.join(directory, "state.mvsj"), "w", encoding="utf-8") as f:
        f.write(state.dumps(indent=2))

    rows = ["<table class='sortable'><tr><th>Rank</th><th>Residue</th>"]
    rows.append("<th>Score</th><th>HADDOCK</th></tr>")
    ranked = sorted(scores.items(), key=lambda item: -item[1])
    for rank, (res, score) in enumerate(ranked, 1):
        kind = "active" if res in active else "passive" if res in passive else ""
        rows.append(
            "<tr class='{}'><td>{:d}</td><td>{}</td>"
            "<td style='background: {}'>{:.5f}</td><td>{}</td></tr>".format(
                kind, rank, res, _score_color(score, smin, smax), score, kind
            )
        )
    rows.append("</table>")

    def residues(ids):
        return ",".join(str(n) for n in sorted(int(res[1:]) for res in ids)) or "None"

    links = [pdb_name] + [os.path.relpath(f, directory) for f in files]
    html = TEMPLATE.substitute(
        molstar=MOLSTAR_VERSION,
        options="".join(
            "<li>{}: {}</li>".format(escape(k), escape(str(v)))
            for k, v in options.items()
        ),
        downloads=" ".join('<a href="{0}">{0}</a>'.format(escape(f)) for f in links),
        active=residues(active),
        passive=residues(passive),
        table="".join(rows),
    )
    fname = os.path.join(directory, "report.html")
    with open(fname, "w", encoding="utf-8") as f:
        f.write(html)
    return fname
