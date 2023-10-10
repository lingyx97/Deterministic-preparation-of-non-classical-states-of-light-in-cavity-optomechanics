import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent/"common"))


from generateOps import *
# from qutipPowerByNumpy import powerByNP as PNP
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import numpy as np
import qutip as qtp
import scipy as scp
import matplotlib._color_data as mcd
# from reconstructEvoSO import recEvoPYSO
from lindblad_numeric import *
from linear_grayscale_colormap import *
from lindblad_trotter import *

# from parseEvoString import NumericStringParser

os.chdir(str(Path(__file__).parent/"tempQU"))

qtp.settings.auto_tidyup = False