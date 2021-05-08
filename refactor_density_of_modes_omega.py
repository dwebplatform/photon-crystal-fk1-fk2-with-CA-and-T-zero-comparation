import matplotlib.pyplot as plt
import numpy as np

import math
from scipy import interpolate, linalg as la
import scipy.signal
from scipy.interpolate import splev, splrep, CubicSpline, UnivariateSpline, interp1d
import pydash as pyd

import sys
sys.path.append('./utils')
from density_modes_utils import createOmegaList,getTransferIndexList,getInterpolatedForDensity, getDensityFromOmega
from matrix_utils import getSuperMatrix,getLeftMatrix, getRightMatrix, getPMatrix, getRefractiveFromList, getReflection, getMatrixAbsDirect, getThreeCalcMatrix
from reader_utils import mutateTxtData, readTransmittanceForTEZeroDegrees, readExperimentalData,readExperimentalDataFromThreeCols, readExperimentalDataForRefractive
from layer_utils import getResultMatrix, matrixForPeriod,matrixForLayer
from collection_utils import  converDictToArray
sys.path.append('./layer_manager')
from layer_with_omega import LayerOmega
 
 
CA_data, PVK_data,Glass_data = (readExperimentalData('./photon_crystal_data/CA.txt'), readExperimentalData('./photon_crystal_data/PVK.txt'),readExperimentalDataFromThreeCols('./photon_crystal_data/Стекло.txt'))

[CA_data,PVK_data,_] = mutateTxtData(CA_data,PVK_data,Glass_data) 

def getSuperMatrix(leftMatrix, glassLayer, rightMatrix):
    [leftMatrixAbsDirect, directDAbsMatrix, rightMatrixAbsDirect] = [getMatrixAbsDirect(
        leftMatrix), glassLayer.getDMatrix(), getMatrixAbsDirect(rightMatrix)]
    return getThreeCalcMatrix(leftMatrixAbsDirect, glassLayer.getDMatrix(), rightMatrixAbsDirect)


def getPlotData(leftMatrix):
    r = getReflection(leftMatrix)
    T = 1 / (leftMatrix[0][0])
    R = leftMatrix[1][0]/leftMatrix[0][0]
    return [r, T, R]
def derivative(f, a, h=0.001):
    return (f(a + h) - f(a - h))/(2*h)


# ns = 1.52
"""длина волны, продабуированная"""
lamdaList = list(range(280, 800, 1))
#
plotArray = []
"""параметры начального луча """
n0 = 1
alfa0 = (math.pi/180)*0

nData = []
yTestData = []
nCAData =[]
omegaList = createOmegaList(280, 800)
c = 300000000

n1, n2, ns, h1, h2, d =(1.475,1.648,1, 50,150,10)
""" массив tx, ty """
txData, tyData = getTransferIndexList(omegaList,n1, n2, ns, h1, h2, d)

(xInterpolated,yInterpolated) = getInterpolatedForDensity(txData, tyData,omegaList)

derivativeXData = []
derivativeYData = []

xInterpolatedData = []

densityOfModesPlot = getDensityFromOmega(200_0,omegaList,xInterpolated, yInterpolated)
plt.plot(omegaList, np.absolute(densityOfModesPlot), 'y')
plt.show()
