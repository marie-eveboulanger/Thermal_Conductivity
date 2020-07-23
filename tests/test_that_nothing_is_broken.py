"""
Unfortunately, it is hard to design unit tests with matplotlib. 
However, this file can be run to have an idea that everything seems fine.
"""

from ThermalConductivity import Analysis
import os

DIR = os.path.dirname(os.path.realpath(__file__))
FILENAME = "Data-TS-0.0T-BOT-SCOC_1903A-2019-04-09.dat"
FILENAME = os.path.join(DIR,FILENAME)

conductivity = Analysis.Conductivity(FILENAME)

conductivity.Plot("kxx", "-o")
