import math
import numpy as np
from numpy import linalg as LA
import sys
sys.path.append('../utils')
from matrix_utils import getPMatrix
from layer_utils import matrixForLayer
# n0*sin(alfa0) = ni *sin(alfai)

class LayerOmega:
    n0 = 1
    alfa0 = (math.pi/180)*0
    def __init__(self, omega, n, h,alfa=None):
        self.c = 3*100000000
        self.n = n
        self.h = h
        self.omega = omega
        if alfa is None:
            pass
        else:
            self.alfa0 = alfa
    def getRefractiveIndex(self):
        return self.n

    def getAngle(self):
        angleSin = (math.sin(self.alfa0) * self.n0) / self.n
        return math.asin(angleSin)

    def getKsi(self):
        Ksi = ((self.n * self.n) - self.n0 * self.n0 *
               math.sin(self.alfa0) * math.sin(self.alfa0))**(0.5)
        return Ksi

    def getK(self):
        return (self.omega/self.c) * self.getKsi()
    def getPMatrix(self):
        return np.array([[1, 1], [self.getKsi(), -self.getKsi()]])
    def getDMatrix(self):
        return np.array([[np.exp(1j * self.getK() * self.h), 0], [0, np.exp(-1j * self.getK() * self.h)]])
    def getPReverseMatrix(self):
        return LA.inv(self.getPMatrix())
    @staticmethod
    def getMatrixWithPInversed(layer,matrix):
      return np.dot(LA.inv(getPMatrix(
        layer.alfa0, layer.n0)), matrix)
    @staticmethod
    # middleLayers = [Layer(lamda,n2,146),Layer(lamda,n1,48),Layer(lamda,n2,215)]
    def getMiddleLayerForResonator(middleLayers):
        middleMatrix = None
        for middleLayer in middleLayers:
            if middleMatrix is None:
                middleMatrix = matrixForLayer(middleLayer)
            else:
                middleMatrix = np.dot(middleMatrix ,matrixForLayer(middleLayer))
        return middleMatrix