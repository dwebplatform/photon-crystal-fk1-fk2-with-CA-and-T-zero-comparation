
def getPlotData(resultMatrix, superResultMatrix):
    r = getReflection(resultMatrix)
    T = 1 / (superResultMatrix[0][0])
    R = superResultMatrix[1][0]/superResultMatrix[0][0]
    return [r, T, R]