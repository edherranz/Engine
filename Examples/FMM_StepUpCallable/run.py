#!/usr/bin/env python

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from ore_examples_helper import OreExample  # noqa

oreex = OreExample(sys.argv[1] if len(sys.argv) > 1 else False)

print("+-----------------------------------------------------+")
print("| FMM step-up callable note (USD-SOFR, 2025-02-10)      |")
print("+-----------------------------------------------------+")

oreex.print_headline("Run ORE: FMM / LSM")
oreex.run("Input/ore.xml")
oreex.print_headline("Run ORE: LGM comparator on the same trade")
oreex.run("Input/ore_lgm.xml")
