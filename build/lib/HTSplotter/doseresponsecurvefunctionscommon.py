import numpy as np
import os
from grupping import Groupping
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
import math
import scipy.interpolate as inter
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit


class DoseResponseFunctions:
    def __init__(self, analysistype):
        # self.analysistype = analysistype

        self.grup = None
        self.response_range = None
        self.response_rangestd = None



