
def readExperimentalData(fileName):
    with open(fileName) as f:
        data = [line.split() for line in f]
    arrayOfRefractiveIndexes = []
    for value in data:
        lamda = int(float(value[0])*100)
        n = value[1]
        if "0.*1j" in n:
            n = n.replace("0.*1j", "0j")
        if "*1j" in n:
            n = n.replace("*1j", "j")
        arrayOfRefractiveIndexes.append(complex(n))
    return arrayOfRefractiveIndexes

def readExperimentalDataFromThreeColsInArray(fileName):
    with open(fileName) as f:
        experimentalData = [line.split() for line in f]
        lamdaWithIndex = []
        for threeCol in experimentalData:
            [lamda, RealN, ImN] = threeCol
            [hungred, decimal] = lamda.split('.')
            if(decimal ==''): # 00
              lamda = int(hungred)*100
            elif(decimal=='1' or decimal=='2' or decimal=='3' or decimal=='4' or decimal=='5' or decimal=='6' or decimal=='7' or decimal=='8' or decimal=='9'):
              lamda = int(hungred)*100 + int(decimal)*10   
            else:
              lamda = 100*int(hungred) + int(decimal)
            ImN = float(ImN)/1000000
            lamdaWithIndex.append(complex(float(RealN), ImN))
        return lamdaWithIndex

  
