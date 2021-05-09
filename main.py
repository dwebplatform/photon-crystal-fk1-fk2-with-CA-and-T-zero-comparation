import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import interpolate, linalg as la
from scipy.interpolate import CubicSpline, UnivariateSpline, interp1d
import pydash as pyd

import sys
sys.path.append('./utils')
from matrix_utils import  getPMatrix, getRefractiveFromList, getReflection, getMatrixAbsDirect, getThreeCalcMatrix
from reader_utils import mutateTxtData, readExperimentalData, readExperimentalDataFromThreeCols,readExperimentalDataFromThreeColsInArray
from layer_utils import getResultMatrix, matrixForPeriod

sys.path.append('./layer_manager')
from layer import Layer

CA_data = readExperimentalData('./photon_crystal_data/CA.txt')
PVK_data = readExperimentalData('./photon_crystal_data/PVK.txt')
Glass_data = readExperimentalDataFromThreeCols('./photon_crystal_data/Стекло.txt')
Glass_dataArray = readExperimentalDataFromThreeColsInArray('./photon_crystal_data/Стекло.txt')

# Glass_dataArray
[CA_data,PVK_data,Glass_data] = mutateTxtData(CA_data,PVK_data,Glass_data) 
def getSuperMatrix(leftMatrix, glassLayer, rightMatrix):
    [leftMatrixAbsDirect, directDAbsMatrix, rightMatrixAbsDirect] = [getMatrixAbsDirect(
        leftMatrix), getMatrixAbsDirect(glassLayer.getDMatrix()), getMatrixAbsDirect(rightMatrix)]
    return getThreeCalcMatrix(leftMatrixAbsDirect, directDAbsMatrix, rightMatrixAbsDirect)

def getPlotData(resultMatrix, superResultMatrix):
    r = getReflection(resultMatrix)
    T = 1 / (superResultMatrix[0][0])
    R = superResultMatrix[1][0]/superResultMatrix[0][0]
    return [r, T, R]

# ns = 1.52
"""длина волны, протабуированная"""
lamdaList = list(range(280, 800, 1))
#
plotArray = []
"""параметры начального луча """
n0 = 1
alfa0 = (math.pi/180)*0
# массив tx, ty
txData = []
tyData = []

yNGlass = []
yTestData = []
nPvkData =[]
nCaData =[]

for lamda in lamdaList:
    experimentalDataList = [ [
        PVK_data, lamda, 1.683], [CA_data, lamda, 1.475]]
    """ns стекло, n1 PVK, n2 CA"""
    [n1, n2] = map(getRefractiveFromList, experimentalDataList)
    if lamda == 431:
      n2 = complex(1.478704894405802,-0.0002410654919080258)
    nCaData.append(n2)
    lamdaAndNs = pyd.find(Glass_dataArray,lambda x: x[0] == lamda)
    if not lamdaAndNs is None:
      [_,ns] = lamdaAndNs
    yNGlass.append(ns.imag)
    """ширина слоя h1, h2, и число периодов"""
    [h1, h2, d] = [50, 150, 10]
    """первый и второй слои"""
    [layerFirst, layerSecond] = [Layer(lamda, n1, h1), Layer(lamda, n2, h2)]
    # [layerFirst, layerSecond] = [Layer(lamda, 1.683, h1), Layer(lamda, 1.475, h2)]
    """матрица всех периодов"""
    resultMatrix = getResultMatrix([layerFirst, layerSecond], d)
    """последний слой стекло"""
    PSMatrix = getPMatrix(layerFirst.alfa0, ns)
    """вычисляем матрицу левую"""
    leftMatrix = np.dot(resultMatrix, PSMatrix)
    glassLayer = Layer(lamda, ns, 1000000)
    """правая матрица"""
    rightMatrix = np.dot(glassLayer.getPReverseMatrix(),
                         getPMatrix(layerFirst.alfa0, layerFirst.n0))
    
    superResultMatrix = getSuperMatrix(leftMatrix, glassLayer, rightMatrix)
    """вычисляем коэффициент отражения"""
    [r, T, R] = getPlotData(resultMatrix, superResultMatrix)
    #t = x(l) + i*y(l)
    txData.append(T.real)
    tyData.append(T.imag)
    yTestData.append(R.real)
    plotArray.append(R.real)
# !расскоментировать:
# plt.plot(plotArray, color='red')
# plt.show()


def derivative(f, a, h=0.001):
    return (f(a + h) - f(a - h))/(2*h)

# интерполировать x по lamda y по lamda создать массив [[x,y],[x,y],[x,y]]

xInterpolated = CubicSpline(lamdaList, txData, bc_type='natural')
yInterpolated = CubicSpline(lamdaList, tyData, bc_type='natural')
densityOfModesPlot = []

for lamda in range(280,800):
    [x,y] = [xInterpolated(lamda),yInterpolated(lamda)]
# for lamda in range(280, 800):
#         [x, y] = [xInterpolated(lamda), yInterpolated(lamda)]
#         [xDerive, yDerive] = [derivative(
#             xInterpolated, lamda), derivative(yInterpolated, lamda)]
#         densityOfModes = (lamda*lamda*(yDerive *
#                                       x - xDerive * y)) / (x * x + y * y)
#         densityOfModesPlot.append(densityOfModes)

# densityOfModesLoop(densityOfModesPlot)
# for lamda in range(280, 800):
    # [x, y] = [xInterpolated(lamda), yInterpolated(lamda)]
    # [xDerive, yDerive] = [derivative(
        # xInterpolated, lamda), derivative(yInterpolated, lamda)]
    # densityOfModes = (lamda*lamda*(yDerive *
                                  #  x - xDerive * y)) / (x * x + y * y)
    # densityOfModesPlot.append(densityOfModes)

# x_new = lamdaList
# y_new = cubicInterpolation(x_new)
# for el in y_new:
#     yDerivative = derivative(cubicInterpolation, el)
# plt.plot(lamdaList,yNGlass,'y')

# plt.plot(lamdaList,densityOfModesPlot,'y')
plt.plot(lamdaList, nCaData, 'y')
plt.show()
