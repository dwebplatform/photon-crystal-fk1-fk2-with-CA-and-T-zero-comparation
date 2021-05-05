import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import interpolate, linalg as la
from scipy.interpolate import CubicSpline, UnivariateSpline, interp1d
import pydash as pyd

import sys
sys.path.append('./utils')
from matrix_utils import getSuperMatrix,getLeftMatrix, getRightMatrix, getPMatrix, getRefractiveFromList, getReflection, getMatrixAbsDirect, getThreeCalcMatrix
from reader_utils import mutateTxtData, readTransmittanceForTEZeroDegrees, readExperimentalData,readExperimentalDataFromThreeCols, readExperimentalDataForRefractive
from layer_utils import getResultMatrix, matrixForPeriod,matrixForLayer
from collection_utils import  converDictToArray
sys.path.append('./layer_manager')
from layer import Layer
 
 
CA_data,PVK_data = (readExperimentalData('./photon_crystal_data/CA.txt'), readExperimentalData('./photon_crystal_data/PVK.txt'))
Glass_data = readExperimentalDataFromThreeCols('./photon_crystal_data/Стекло.txt')


[CA_data,PVK_data,_] = mutateTxtData(CA_data,PVK_data,Glass_data) 

"""массивы значений показателя преломления для CA, PVK, и стекла"""
CADataList, PVKDataList,glassDataList = (converDictToArray(CA_data),converDictToArray(PVK_data),readExperimentalDataForRefractive('./photon_crystal_data/Стекло.txt'))

def getPlotData(allLayersMatrix, superResultMatrix):
    r = getReflection(allLayersMatrix)
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
nCaData = []
 

"""ns стекло, n1 PVK, n2 CA"""
for (lamda, n1, n2, ns) in zip(lamdaList, PVKDataList, CADataList, glassDataList):
    if lamda == 431:
      n2 = complex(1.478704894405802,-0.0002410654919080258)
    nCaData.append(n1)
    """ширина слоя h1, h2, и число периодов"""
    [h1, h2, d] = [45.5, 119,10]
    """первый , второй слои и слой стекла """
    [layerFirst, layerSecond, glassLayer] = [Layer(lamda, n1, h1), Layer(lamda, n2, h2),Layer(lamda, ns, 1000000)]
    # [layerFirst, layerSecond] = [Layer(lamda, 1.683, h1), Layer(lamda, 1.475, h2)]
    """allLayersMatrix - матрица всех периодов PSMatrix - последний слой стекло  """ 
    allLayersMatrix, PSMatrix = (getResultMatrix([layerFirst, layerSecond], d), getPMatrix(layerFirst.alfa0, ns))
    
    """левая и правая матрица"""
    leftMatrix, rightMatrix = (getLeftMatrix(allLayersMatrix, PSMatrix), getRightMatrix(glassLayer,layerFirst))
    superResultMatrix = getSuperMatrix(leftMatrix, glassLayer, rightMatrix)
    """вычисляем коэффициент отражения"""
    [r, T, R] = getPlotData(allLayersMatrix, superResultMatrix)
    #t = x(l) + i*y(l)
    txData.append(T.real)
    print(T)
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

    # with open(fileName) as f:
    #     experimentalData = [line.split() for line in f]
    #     lamdaWithIndex = []
    #     for threeCol in experimentalData:


zeroTData = readTransmittanceForTEZeroDegrees('TE_пропускание_длинноволновый_0_град.txt')

plt.plot(lamdaList, zeroTData, linestyle='--')

plt.show()
