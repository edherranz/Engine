#!/usr/bin/env python

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from ore_examples_helper import OreExample  # noqa

oreex = OreExample(sys.argv[1] if len(sys.argv) > 1 else False)

print("+-----------------------------------------------------+")
print("| Callable bond: LGM Grid drops the first-event coupon |")
print("+-----------------------------------------------------+")

oreex.print_headline("LGM callable bond engine, event-time grid (Grid)")
oreex.run("Input/ore_grid.xml")
oreex.print_headline("LGM callable bond engine, finite differences (FD)")
oreex.run("Input/ore_fd.xml")
oreex.print_headline("With the stock EUR credit curve: Grid, FD and MC (the MC engine needs a credit curve)")
oreex.run("Input/ore_grid_credit.xml")
oreex.run("Input/ore_fd_credit.xml")
oreex.run("Input/ore_mc_credit.xml")
