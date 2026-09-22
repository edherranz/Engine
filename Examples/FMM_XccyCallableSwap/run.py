#!/usr/bin/env python
"""Cancellable cross-currency swap template: stock ORE CrossAssetModel comparator runs, driven
through the fork's ORE Python package (no executable)."""

import os
import sys
import time

from ORE import OREApp, Parameters

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = [
    ("stock ORE CrossAssetModel (LGM x2 + FX, MC): swap + cancellation right", "Input/ore_cam.xml"),
    ("consistency check: the right with the USD floating leg telescoped away", "Input/ore_cam_reduced.xml"),
]


def run(ore_xml):
    params = Parameters()
    params.fromFile(ore_xml)
    app = OREApp(params, False)
    t0 = time.time()
    app.run()
    errors = app.getErrors()
    if errors:
        raise RuntimeError("\n".join(errors))
    npv = app.getReport("npv")
    col = {npv.header(i): i for i in range(npv.columns())}
    ids = npv.dataAsString(col["TradeId"])
    values = npv.dataAsReal(col["NPV(Base)"])
    for i, v in zip(ids, values):
        print("  %-45s %14.2f USD" % (i, v))
    print("  run time %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    os.chdir(HERE)
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for title, ore_xml in RUNS:
        if only and only not in ore_xml:
            continue
        print("+ " + title)
        run(ore_xml)
