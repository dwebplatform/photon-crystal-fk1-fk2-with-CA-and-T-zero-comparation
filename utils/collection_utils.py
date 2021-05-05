
def converDictToArray(someDictionary):
    resultList =[]
    for item in sorted(someDictionary.keys()):
        resultList.append(complex(someDictionary[item]))    
    return resultList