"""
roman_cutouts.cli
~~~~~~~~~~~~~~~~~
Command-line interface for the Roman WFI cutout pipeline.

Entry point: ``roman-make-cutout`` (after ``pip install -e .``)
or ``python -m roman_cutouts``.
"""
import argparse
import logging

from .pipeline import run_cutout_pipeline


def build_parser():
    """Return the argument parser (separated so tests can reuse it)."""
    parser = argparse.ArgumentParser(
        prog="roman-make-cutout",
        description=(
            "Query MAST for Roman WFI exposures and extract sky cutout cubes. "
            "Observation hierarchy flags (--execution-plan … --visit) narrow the "
            "query to a single field; omit them to retrieve all fields for the "
            "program/detector/filter combination."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── Target position ────────────────────────────────────────────────────────
    pos = parser.add_argument_group("Target position")
    pos.add_argument("--ra", type=float, nargs="+", required=True,
                     help="Right Ascension(s) of the cutout centre(s) (deg, ICRS).")
    pos.add_argument("--dec", type=float, nargs="+", required=True,
                     help="Declination(s) of the cutout centre(s) (deg, ICRS); "
                          "same count as --ra.")

    # ── Data selection ─────────────────────────────────────────────────────────
    data = parser.add_argument_group("Data selection")
    data.add_argument("--program-id", type=int, required=True,
                      help="Roman program (proposal) ID.")
    data.add_argument("--detector", "--sca", dest="detector", type=int, required=True,
                      help="WFI detector number 1–18 (also accepted as --sca).")
    data.add_argument("--filter", type=str, default="F146",
                      help="WFI optical element.")
    data.add_argument("--file-type", type=str, default="l2",
                      choices=["l1", "l2", "l4"],
                      help="Data level: l1 (uncal), l2 (cal), l4 (catalog).")

    # ── Observation hierarchy ──────────────────────────────────────────────────
    hier = parser.add_argument_group(
        "Observation hierarchy",
        "Narrow the MAST query to a specific sky field. "
        "File names encode: r{PPPPP}{CC}{AAA}{SSS}{OOO}{VVV}_{eeee}_… "
        "where PPPPP=program, CC=execution_plan, AAA=pass, SSS=segment, "
        "OOO=observation (primary field ID), VVV=visit (dither).",
    )
    hier.add_argument("--execution-plan", type=int, default=None,
                      help="Execution plan within the program (CC).")
    hier.add_argument("--pass-num", type=int, default=None,
                      help="Pass / epoch number (AAA).")
    hier.add_argument("--segment", type=int, default=None,
                      help="Segment within the pass (SSS).")
    hier.add_argument("--observation", type=int, default=None,
                      help="Observation number — primary sky-field ID (OOO).")
    hier.add_argument("--visit", type=int, default=None,
                      help="Visit (dither position) within the observation (VVV).")

    # ── Cutout and output ──────────────────────────────────────────────────────
    out = parser.add_argument_group("Cutout and output")
    out.add_argument("--cutout-size", type=str, nargs="+", default=["25"],
                     metavar="SIZE",
                     help="Cutout size in pixels: N (square) or NY,NX. One value is "
                          "shared by all targets; otherwise give one per target.")
    out.add_argument("--max-files", type=int, default=None,
                     help="Process at most this many files (useful for testing).")
    out.add_argument("--output-path", type=str, default=None,
                     help="Directory for the output ASDF file(s). Omit to skip saving.")
    out.add_argument("--verbose", action="store_true",
                     help="Print per-file pixel-coordinate messages.")

    par = parser.add_argument_group("Parallelization")
    par.add_argument("--n-workers", type=int, default=8,
                     help="Number of parallel workers (1 = serial).")
    par.add_argument("--backend", choices=["thread", "process"], default="thread",
                     help="Parallel backend.")

    return parser


def _parse_sizes(tokens):
    """Convert ``['25', '40,30']`` to ``[25, (40, 30)]``; a single entry is unwrapped."""
    sizes = []
    for tok in tokens:
        parts = [int(p) for p in tok.split(",")]
        if len(parts) == 1:
            sizes.append(parts[0])
        elif len(parts) == 2:
            sizes.append((parts[0], parts[1]))
        else:
            raise ValueError(f"Invalid --cutout-size value {tok!r}; use N or NY,NX.")
    return sizes[0] if len(sizes) == 1 else sizes


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    )

    args = build_parser().parse_args()

    run_cutout_pipeline(
        ra=args.ra,
        dec=args.dec,
        program_id=args.program_id,
        detector=args.detector,
        filter_name=args.filter,
        cutout_size=_parse_sizes(args.cutout_size),
        file_type=args.file_type,
        execution_plan=args.execution_plan,
        pass_num=args.pass_num,
        segment=args.segment,
        observation=args.observation,
        visit=args.visit,
        max_files=args.max_files,
        output_path=args.output_path,
        verbose=args.verbose,
        n_workers=args.n_workers,
        backend=args.backend,
    )
