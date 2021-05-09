import numpy as np
import math

from plot_utils import getPlotData
from layer_utils import getResultMatrix, getMatrixPower
from matrix_utils import getPMatrix, getReflection, getThreeCalcMatrix, getLeftMatrix, getRightMatrix, getSuperMatrix
import scipy.signal
from scipy.interpolate import CubicSpline, UnivariateSpline, interp1d
import sys
sys.path.append('./layer_manager')
from layer_with_omega import LayerOmega
"""c - скорость света """
c = 300000000
def createOmegaList(lamdamin, lamdamax):
  omegaList = []
  minOmega,maxOmega = (2*math.pi*c/lamdamax,2*math.pi*c/lamdamin)
  deltaOmega = (maxOmega - minOmega)/1000
  count = minOmega
  while(count<maxOmega):
    count+=deltaOmega
    omegaList.append(count)
  return omegaList
""" """


def getTransferForMatrix(matrix):
    r = getReflection(matrix)
    T = 1 / (matrix[0][0])
    R = matrix[1][0]/matrix[0][0]
    return [r, T, R]


def getTransferForIndexListResonator(omegaList,n1, n2, ns, h1, h2, d, alfa=None):
  txData,tyData = ([],[]) 
  for omega in omegaList:
    [glassLayer] = [LayerOmega(omega, ns, 1000000,alfa)]

    """ левое зеркало""" 
      # ПВК 48 нм и СА 137 нм
    [hLeftPVKLayer, hLeftCALayer] =[45.5, 119]
    [leftLayerFirst,leftLayerSecond] = [LayerOmega(omega, n1, hLeftPVKLayer,alfa), LayerOmega(omega, n2, hLeftCALayer,alfa)]

    """ правое зеркало"""
      # ПВК 48 нм и СА 137 нм
    [hRightPVKLayer, hRightCALayer] =[48, 137]
    [rightLayerFirst,rightLayerSecond] = [LayerOmega(omega, n1, hRightPVKLayer,alfa), LayerOmega(omega, n2, hRightCALayer,alfa)]


    """allLayersMatrix - матрица всех периодов PSMatrix - последний слой стекло  """ 
      
    """потом слой СА 146 нм (т.е. 2 слоя СА на стыке, так получилось), 
    потом слой ПВК 48 нм, 
    потом слой СА 215 нм, далее"""
    middleMatrix = LayerOmega.getMiddleLayerForResonator([LayerOmega(omega,n2,146,alfa),LayerOmega(omega,n1,48,alfa),LayerOmega(omega,n2,215,alfa)])
    allLeftMatrix, allRightMatrix = (getMatrixPower([leftLayerFirst, leftLayerSecond], 7), getMatrixPower([rightLayerFirst, rightLayerSecond], 8))
    layersWithMiddle = getThreeCalcMatrix(allLeftMatrix,middleMatrix,allRightMatrix)

    allLayersMatrix = LayerOmega.getMatrixWithPInversed(leftLayerFirst,layersWithMiddle)
    PSMatrix = getPMatrix(glassLayer.alfa0, ns)
    """левая и правая матрица"""
    leftMatrix, rightMatrix = (getLeftMatrix(allLayersMatrix, PSMatrix), getRightMatrix(glassLayer,glassLayer))
    superResultMatrix = getSuperMatrix(leftMatrix, glassLayer, rightMatrix)
    """вычисляем коэффициент отражения"""
    [r, T, R] = getPlotData(allLayersMatrix, superResultMatrix)
    txData.append(T.real)
    tyData.append(T.imag)
  return (txData,tyData) 

def getTransferIndexList(omegaList,n1, n2, ns, h1, h2, d, alfa=None):
  txData,tyData = ([],[]) 
  for omega in omegaList:
    [layerFirst, layerSecond] = [LayerOmega(omega, n1, h1,alfa), LayerOmega(omega, n2, h2,alfa)]
    """матрица всех периодов"""
    resultMatrix = getResultMatrix([layerFirst, layerSecond], d)
    """последний слой стекло"""
    PSMatrix = getPMatrix(layerFirst.alfa0, ns)
    """вычисляем матрицу левую"""
    leftMatrix = np.dot(resultMatrix, PSMatrix)
    [r, T, R] = getTransferForMatrix(leftMatrix)
    txData.append(T.real)
    tyData.append(T.imag)
  return (txData,tyData)


def derivative(f, a, h=0.001):
    return (f(a + h) - f(a - h))/(2*h)

def getDensityFromOmegaResonator(D,omegaList,xInterpolated,yInterpolated):
  densityOfModesArray = []
  for omega in omegaList:
    [x, y] = [xInterpolated(omega), yInterpolated(omega)]
    [xDerive, yDerive] = [derivative(xInterpolated, omega), derivative(yInterpolated, omega)]
    densityOfModes = (1/D)*((yDerive *x - xDerive * y)) /(x * x + y * y)
    densityOfModesArray.append(densityOfModes) 
  return  densityOfModesArray

def getDensityFromOmega(D,omegaList,xInterpolated,yInterpolated):
  densityOfModesArray = []
  for omega in omegaList:
    [x, y] = [xInterpolated(omega), yInterpolated(omega)]
    [xDerive, yDerive] = [derivative(xInterpolated, omega), derivative(yInterpolated, omega)]
    densityOfModes = (1/D)*((yDerive *x - xDerive * y)) /(x * x + y * y)
    densityOfModesArray.append(densityOfModes) 
  return  densityOfModesArray
def getInterpolatedForDensity(txData,tyData, omegaList):
  new_tx_Data = scipy.signal.savgol_filter(txData, 37, 4)
  new_ty_Data = scipy.signal.savgol_filter(tyData, 37, 4)
  xInterpolated = CubicSpline(omegaList, new_tx_Data, bc_type='natural')
  yInterpolated = CubicSpline(omegaList, new_ty_Data, bc_type='natural')  
  return (xInterpolated,yInterpolated)
