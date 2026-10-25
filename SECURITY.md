# Security and data handling

## Where your sequences go

Nowhere. The application binds to localhost, every analysis runs in your own
process, and no sequence, structure or report is transmitted anywhere. There is
no telemetry, no error reporting and no usage counter.

Reference and model downloads reach the network during setup or first use;
they do not carry your submitted sequence or structure:

- `setup/build_reference_data.py` downloads Oxford's Thera-SAbDab CSV.
- `setup/fetch_deepsp.py` downloads DeepSP's weights and training labels.
- `setup/fetch_constant_regions.py` queries UniProt for constant regions.
- Model checkpoints (ImmuneBuilder, AntiBERTy) are fetched on first use.
- `scripts/setup_stability.py` downloads pinned ThermoMPNN source, model
  checkpoints and license notices. Stability inference itself runs locally.

If you are screening unpublished molecules, this matters, so it is worth
stating plainly: run the setup steps once, and after that the tool needs no
network at all. You can verify that with your own firewall.

## Browser access

Use `http://127.0.0.1:8000` (or the port printed by the launcher).
The server accepts loopback Host names only (`127.0.0.1`, `localhost`, or
`[::1]`) to prevent DNS rebinding through an unrelated website. Browser
requests carrying an Origin must match the server's scheme, host, and port;
opaque origins and cross-site browser requests are rejected. Opening the home
page through a top-level browser navigation is allowed, so a link or launcher
can open Atlas; this exception does not apply to API access. Local CLI
requests without an Origin remain supported. These protections do not add
user authentication or make a network deployment safe.

## What `.gitignore` is doing

Analysis output is excluded from version control by default — reports, PDBs,
PDFs, result directories — because that output contains the sequences you put
in. If you fork this for your own molecules, keep those rules. The most likely
way to disclose an unpublished sequence with this tool is to commit a report.

`tools/make_oss_export.py` is the gate for anything going public. It refuses to
emit a tree containing restricted sequences or molecule names, and it checks
that the files the app needs at runtime are present. Run it before any public
push.

## Reporting a vulnerability

Open a GitHub security advisory on this repository ("Security" tab → "Report a
vulnerability"), which keeps the report private until it is resolved. Please do
not open a public issue for anything that would let someone read or exfiltrate
another user's sequences.

Expect an acknowledgement within a week. This is a research tool maintained by
one person, so a fix may take longer than that, but you will be told where it
stands.

## What is explicitly out of scope

This tool is designed to be run by the person whose data it is, on their own
machine. It has no authentication, no authorization and no multi-user model,
because it has no users other than the person who started it. Do not expose the
port to a network you do not control, and do not run it as a shared service: if
you do, anyone who can reach the port can read every analysis on it. That is
not a vulnerability in the tool, it is a deployment the tool was never built
for, and `run.py` binds to localhost to make the intended shape obvious.
