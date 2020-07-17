"""
Unfortunately, it is hard to design unit tests with matplotlib. 
However, this file can be run to have an idea that everything seems fine.
"""

from ThermalConductivity import Analysis

FILENAME = "Data-TS-0.0T-BOT-SCOC_1903A-2019-04-09.dat"

conductivity = Analysis.Conductivity(FILENAME)

conductivity.Plot("kxx", "-o")
