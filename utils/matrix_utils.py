import numpy as np
import math


def getPMatrix(angle, n):
    return np.array(
        [[1, 1], [n*math.cos(angle), -n*math.cos(angle)]])


def getRefractiveFromList(refractiveData):
    [experimentalData, lamda, defaultRefractiveIndex] = refractiveData
    n = defaultRefractiveIndex
    if lamda in experimentalData:
        n = complex(experimentalData[lamda])
    return n


def getThreeCalcMatrix(a, b, c):
    res = np.dot(a, b)
    res = np.dot(res, c)
    return res
def getMatrixAbsDirect(matrix):
    matrix[0][0] = np.absolute( matrix[0][0]) *np.absolute( matrix[0][0]) 
    matrix[0][1] = np.absolute( matrix[0][1]) * np.absolute( matrix[0][1])
    matrix[1][0] = np.absolute( matrix[1][0])* np.absolute( matrix[1][0])
    matrix[1][1] = np.absolute( matrix[1][1])* np.absolute( matrix[1][1])
    return matrix


def getMatrixForLayer(pMatrix, dMatrix, pMatrixReveresed):
    res = np.dot(pMatrix, dMatrix)
    res = np.dot(res, pMatrixReveresed)
    return res

def getMatrixForLayerInstance(layer):
    return  getMatrixForLayer(layer.getPMatrix(),  layer.getDMatrix(), layer.getPReverseMatrix())
def getReflection(matrix):
    A = matrix[0][0]
    C = matrix[1][0]
    amplitude = C / A
    r = amplitude * amplitude.conjugate()
    return r



def getLeftMatrix(allLayersMatrix, PSMatrix):
    return np.dot(allLayersMatrix, PSMatrix) 
def getRightMatrix(glassLayer,layerFirst):
    return np.dot(glassLayer.getPReverseMatrix(), getPMatrix(layerFirst.alfa0, layerFirst.n0))



def getSuperMatrix(leftMatrix, glassLayer, rightMatrix):
    [leftMatrixAbsDirect, directDAbsMatrix, rightMatrixAbsDirect] = [getMatrixAbsDirect(
        leftMatrix), getMatrixAbsDirect(glassLayer.getDMatrix()), getMatrixAbsDirect(rightMatrix)]
    return getThreeCalcMatrix(leftMatrixAbsDirect, directDAbsMatrix, rightMatrixAbsDirect)
