#!/usr/bin/env python

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from ore_examples_helper import OreExample  # noqa

oreex = OreExample(sys.argv[1] if len(sys.argv) > 1 else False)

print("+-----------------------------------------------------+")
print("| FMM Bermudan swaption (USD-SOFR, 2025-02-10)        |")
print("+-----------------------------------------------------+")

oreex.print_headline("Run ORE: FMM / LSM (Longstaff-Schwartz with dual bound)")
oreex.run("Input/ore.xml")
oreex.print_headline("Run ORE: LGM / Grid on the same trade (comparison)")
oreex.run("Input/ore_lgm.xml")
oreex.print_headline("Run ORE: FMM / LSM with the coterminal ATM calibration strategy")
oreex.run("Input/ore_atm.xml")
oreex.print_headline("Run ORE: FMM / LSM calibrated jointly to the SOFR optionlets and the coterminal swaptions")
oreex.run("Input/ore_joint.xml")
