import pydash as pyd

def readExperimentalDataTest(fileName):
    with open(fileName) as f:
        data = [line.split() for line in f]
    objectData = {}
    for value in data:
        value[0] = int(float(value[0])*100)
        if "0.*1j" in value[1]:
            value[1] = value[1].replace("0.*1j", "0j")
        if "*1j" in value[1]:
            value[1] = value[1].replace("*1j", "j")
        objectData[value[0]] = value[1]
    return objectData

def readExperimentalData(fileName):
    with open(fileName) as f:
        data = [line.split() for line in f]
    objectData = {}
    for value in data:
        value[0] = int(float(value[0])*100)
        if "0.*1j" in value[1]:
            value[1] = value[1].replace("0.*1j", "0j")
        if "*1j" in value[1]:
            value[1] = value[1].replace("*1j", "j")
        objectData[value[0]] = value[1]
        # objectData.append(value[1])
    return objectData

def readExperimentalDataFromThreeCols(fileName):
    with open(fileName) as f:
        experimentalData = [line.split() for line in f]
        finalRes = {}
        for threeCol in experimentalData:
            [lamda, RealN, ImN] = threeCol
            lamda = int(float(lamda) * 100)
            ImN = float(ImN)/1000000
            finalRes[lamda] = complex(float(RealN), ImN)
    return finalRes



def readExperimentalDataFromThreeColsInArray(fileName):
    with open(fileName) as f:
        experimentalData = [line.split() for line in f]
        #  массив [[lamda1, n1],[lamda2, n2],... и тд]
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
            lamdaWithIndex.append([lamda,complex(float(RealN), ImN)])
        return lamdaWithIndex


def readTransmittanceForTEZeroDegrees(fileName):
    with open(fileName) as f:
        experimentalData = [line.split() for line in f]
        lamdaWithT = []
        for twoColumns in experimentalData:
            [lamda, T] = twoColumns
            lamdaWithT.append(float(T))
        return lamdaWithT
def readExperimentalDataForRefractive(fileName):
    with open(fileName) as f:
        experimentalData = [line.split() for line in f]
        refractiveIndexList = []
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
            refractiveIndexList.append(complex(float(RealN), ImN))
        return refractiveIndexList






def mutateTxtData(CA_data,PVK_data,Glass_data):
    CA_data[477] = 1.4738791595816867-0.00020608105693277592j
    CA_data[402] = 1.482839407001352-0.00027069787699413527j
    CA_data[406] = 1.482203674613889-0.00026614862684739556j
    CA_data[410] = 1.4815907222447684-0.00026176121947212277j
    CA_data[414] = 1.4809994306485148-0.00025752705922949323j
    CA_data[427] = 1.4792150807638462-0.0002447313291047784j
    CA_data[431] = 1.7413909681827986+0.0j
    CA_data[435] = 1.478211528088493-0.0002375164564511407j
    CA_data[439] = 1.4777342235403816-0.00023407862720977422j
    CA_data[452] = 1.4762855428905464-0.00022361247704396093j
    CA_data[456] = 1.4758689981871864-0.00022059251940484425j
    CA_data[460] = 1.4754651826208782-0.00021765956861090325j
    CA_data[464] = 1.4750735655830536-0.00021480983340122403j
    CA_data[481] = 1.4735340875657656-0.00020354746476096083j
    CA_data[485] = 1.4731988503637146-0.00020108041739213634j
    CA_data[489] = 1.472873066455529-0.0001986772585081511j
    CA_data[502] = 1.4718752901922678-0.0001912786868129822j
    CA_data[506] = 1.4715858428792674-0.0001891205090177918j
    CA_data[510] = 1.4713041296501064-0.0001870143734681904j

    PVK_data[402] = 1.764790867611399+0.0j
    PVK_data[406] = 1.7620832800876554+0.0j
    PVK_data[410] = 1.7580276181779861+0.0j
    PVK_data[414] = 1.7531363128154878+0.0j
    PVK_data[427] = 1.743605851666202+0.0j
    PVK_data[431] = 1.7413909681827986+0.0j
    PVK_data[435] = 1.73831557532483+0.0j
    PVK_data[439] = 1.7354843169376928+0.0j
    PVK_data[452] = 1.728232258045814+0.0j
    PVK_data[456] = 1.7257361619959835+0.0j
    PVK_data[460] = 1.7232944784091628+0.0j
    PVK_data[464] = 1.7214130179449232+0.0j
    PVK_data[477] = 1.7160897441561955+0.0j
    PVK_data[481] = 1.7143088604822019+0.0j
    PVK_data[485] = 1.7125932416629286+0.0j
    PVK_data[489] = 1.7111135391305181+0.0j
    PVK_data[502] = 1.7072385293843606+0.0j
    PVK_data[506] = 1.7060960037837483+0.0j
    PVK_data[510] = 1.7048511290735953+0.0j


    Glass_data[402] = complex(1.5144819344853788, - 0.49276820442878305)
    Glass_data[406] = complex(	1.5135919408582557, - 0.48533487792528196)
    Glass_data[410] = complex(1.5127121990401105, - 0.5245137895205729)
    Glass_data[414] = complex(1.511841923791581, - 0.567021114183855)
    Glass_data[427] = complex(1.509070088193674, - 0.694799291073842)
    Glass_data[431] = complex(1.5082320902673476, - 0.794056020323045)
    Glass_data[435] = complex(1.5073999808684662, - 0.777866028008483)
    Glass_data[439] = complex(1.5065731652849772, - 0.8462197926242889)
    Glass_data[452] = complex(1.5039158112576005, - 0.883384050670314)
    Glass_data[456] = complex(1.5031053682335564, - 0.865335819425766)
    Glass_data[460] = complex(1.5022974409809216, - 0.833989698900241)
    Glass_data[464] = complex(1.5014915552256625, - 0.833964421083963)
    Glass_data[477] = complex(1.4988810600783196, - 0.7824126427541189)
    Glass_data[481] = complex(1.4980788789452524, - 0.725779332352428)
    Glass_data[485] = complex(1.4972764649771644, - 0.770099070417813)
    Glass_data[489] = complex(1.496473418891503,	-0.779700980839764)
    Glass_data[502] = complex(1.4938544882900369, - 0.7298811080430601)
    Glass_data[506] = complex(	1.4930444926128519, - 0.724173469145571)
    Glass_data[510] = complex(1.4922318967479695, - 0.711814743521476)
    return [CA_data,PVK_data,Glass_data]