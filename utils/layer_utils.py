import numpy as np
from numpy import linalg as LA
from matrix_utils import getPMatrix

def matrixForPeriod(layers):
    [M1Matrix, M2Matrix] = map(matrixForLayer, layers)
    MofTwoMatrix = np.dot(M1Matrix, M2Matrix)
    return MofTwoMatrix
def matrixForLayer(layer):
    result = np.dot(layer.getPMatrix(), layer.getDMatrix())
    result = np.dot(result, layer.getPReverseMatrix())
    return result
# layers =[layerFirst, layerSecond]
def getResultMatrix(layers, d):
    [layerFirst, layerSecond] = layers
    AllMMatrix = LA.matrix_power(matrixForPeriod(layers), d)
    resultMatrix = np.dot(LA.inv(getPMatrix(
        layerFirst.alfa0, layerFirst.n0)), AllMMatrix)
    return resultMatrix


# layers =[layerFirst, layerSecond]
def getMatrixPower(layers,d):
    [layerFirst, layerSecond] = layers
    AllMMatrix = LA.matrix_power(matrixForPeriod(layers), d)
    return AllMMatrix

def getResultMatrixWithCa(layers, lastCaLayer, d):
    [layerFirst, layerSecond] = layers
    AllMMatrix = LA.matrix_power(matrixForPeriod(layers), d)
    leftMultiplication =np.dot( LA.inv(getPMatrix(layerFirst.alfa0, layerFirst.n0)),matrixForLayer(lastCaLayer))
    resultMatrix = np.dot(leftMultiplication, AllMMatrix)
    return resultMatrix
