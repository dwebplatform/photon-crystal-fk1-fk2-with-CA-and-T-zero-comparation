from matrix_utils import  getReflection
def getPlotData(resultMatrix, superResultMatrix):
    r = getReflection(resultMatrix)
    T = 1 / (superResultMatrix[0][0])
    print(T)
    R = superResultMatrix[1][0]/superResultMatrix[0][0]
    return [r, T, R]