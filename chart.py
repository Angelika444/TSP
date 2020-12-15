import matplotlib.pyplot as plt
import math as m
import statistics

algorithNo = 6

def result2():
    file = open("results2.txt", "r")
    instanceNo = int(file.readline())
    instanceNames = []
    optimal = []
    algorithmsNames = []
    results = []
    bestResults = []
    #worstResults = []
    meanResults = []
    deviationResults = []
    solutionNo = []
    steps = []
    times = []
    iterationTimes = []
    for i in range(instanceNo):
        instanceNames.append(file.readline()[:-5])
        optimal.append(float(file.readline()))
        algorithmsNames.append([])
        results.append([])
        bestResults.append([])
        #worstResults.append([])
        solutionNo.append([])
        steps.append([])
        times.append([])
        iterationTimes.append([])
        meanResults.append([])
        deviationResults.append([])
        for j in range(algorithNo):
            algorithmsNames[i].append(file.readline())
            results[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            bestResults[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            #worstResults[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            solutionNo[i].append(sum([int(x) for x in file.readline().split()]) / 10)
            steps[i].append(sum([int(x) for x in file.readline().split()]) / 10)
            times[i].append(float(file.readline()))
            iterationTimes[i].append(float(file.readline()))
        
        
            #oblicza srednia i odchylenie z rozwiazan
            mean = sum(results[i][j]) / 10
            deviation = 0
            for x in results[i][j]:
                deviation += (x - mean) * (x - mean)
            deviationResults[i].append(m.sqrt(deviation / 10))
            meanResults[i].append(mean)
    
    file.close()
    
    dataBest = []
    dataMean = []
    #dataWorst = []
    dataDeviation = []
    dataSolutionNo = []
    dataSteps = []
    dataTimes = []
    dataIterationTimes = []
    for i in range(algorithNo):
        dataBest.append([])
        dataMean.append([])
        dataDeviation.append([])
        #dataWorst.append([])
        dataSolutionNo.append([])
        dataSteps.append([])
        dataTimes.append([])
        dataIterationTimes.append([])
        for j in range(len(bestResults)):
            dataBest[i].append(bestResults[j][i][-1])
            dataMean[i].append(meanResults[j][i])
            dataDeviation[i].append(deviationResults[j][i])
            #dataWorst[i].append(worstResults[j][i][-1])
            dataSolutionNo[i].append(solutionNo[j][i])
            dataSteps[i].append(steps[j][i])
            dataTimes[i].append(times[j][i])
            dataIterationTimes[i].append(iterationTimes[j][i])
            
    
    plt.figure()
    for i in range(algorithNo):
        plt.plot(dataBest[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Best results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2best.pdf', bbox_inches='tight')
    #plt.show()
    

    plt.figure()
    x = [i for i in range(instanceNo)]
    for i in range(algorithNo):
        plt.errorbar(x, dataMean[i], yerr=dataDeviation[i], fmt='o', label = algorithmsNames[0][i]) 
    plt.xticks(range(instanceNo), instanceNames)
    plt.title('Mean results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.yscale('log')
    plt.savefig('charts/2mean.pdf', bbox_inches='tight')
    #plt.show()
    
    """plt.figure()
    for i in range(algorithNo):
        plt.plot(dataWorst[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Worst results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2worst.pdf', bbox_inches='tight')
    #plt.show()"""
    
    plt.figure()
    for i in range(algorithNo):
        plt.plot(dataSolutionNo[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Solution number')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2solutionNo.pdf', bbox_inches='tight')
    #plt.show()
    
    plt.figure()
    for i in range(3, 5):
        plt.plot(dataSteps[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Number of steps')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2steps.pdf', bbox_inches='tight')
    #plt.show()
    
    plt.figure()
    for i in range(algorithNo):
        plt.plot(dataTimes[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Mean time of 1 run')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2time.pdf', bbox_inches='tight')
    #plt.show()
    
    plt.figure()
    for i in range(algorithNo):
        plt.plot(dataIterationTimes[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Mean iteration time')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2iterationTime.pdf', bbox_inches='tight')
    #plt.show()
    
    return optimal, dataBest, dataMean, dataDeviation, dataTimes
    
    
def result3(optimal):
    file = open("results3.txt", "r")
    instanceNo = int(file.readline())
    iterationNo = int(file.readline())
    instanceNames = []
    algorithmsNames = ['G', 'S']
    startSolutionG = []
    startSolutionS = []
    finishSolutionG = []
    finishSolutionS = []
    for i in range(instanceNo):
        instanceNames.append(file.readline()[:-5])
        file.readline()
        startSolutionG.append([])
        finishSolutionG.append([])
        for j in range(iterationNo):
            line = file.readline().split()
            startSolutionG[i].append((float(line[0]) - optimal[i]) / optimal[i])
            finishSolutionG[i].append((float(line[1]) - optimal[i]) / optimal[i])
        
        file.readline()
        startSolutionS.append([])
        finishSolutionS.append([])
        for j in range(iterationNo):
            line = file.readline().split()
            startSolutionS[i].append((float(line[0]) - optimal[i]) / optimal[i])
            finishSolutionS[i].append((float(line[1]) - optimal[i]) / optimal[i])
    file.close()
    
    file = open("correlation.txt", "w")
    corS = []
    corG = []
    for i in range(instanceNo):
        pomS = []
        pomG = []
        for j in range(len(startSolutionS[i])):
            pomS.append(startSolutionS[i][j] * finishSolutionS[i][j])
            pomG.append(startSolutionG[i][j] * finishSolutionG[i][j])
        covS = (sum(pomS) - sum(startSolutionS[i]) * sum(finishSolutionS[i]) / len(pomS)) / len(pomS)
        covG = (sum(pomG) - sum(startSolutionG[i]) * sum(finishSolutionG[i]) / len(pomG)) / len(pomG)
        corS.append(covS / (statistics.stdev(startSolutionS[i]) - statistics.stdev(finishSolutionS[0])))
        corG.append(covG / (statistics.stdev(startSolutionG[i]) - statistics.stdev(finishSolutionG[0])))
    file.write('instancja G S \n')
    for i in range(instanceNo):
        file.write(instanceNames[i] + ' ' + str(corG[i]) + ' ' + str(corS[i]) + '\n')
    print(instanceNames)
    print(corS)
    print(corG)
    file.close()
    
    for i in range(instanceNo):
        plt.figure()
        plt.scatter(startSolutionG[i], finishSolutionG[i], label = 'G')
        plt.scatter(startSolutionS[i], finishSolutionS[i], label = 'S')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title(instanceNames[i])
        plt.xlabel('starting solution quality')
        plt.ylabel('finish solution quality')
        plt.savefig('charts/3' + instanceNames[i] + '.pdf', bbox_inches='tight')
    
    
def result4(optimal):
    file = open("results4.txt", "r")
    instanceNo = int(file.readline())
    iterationNo = int(file.readline())
    instanceNames = []
    algorithmsNames = ['G', 'S']
    bestSolutionG = []
    bestSolutionS = []
    meanSolutionG = []
    meanSolutionS = []
    for i in range(instanceNo):
        instanceNames.append(file.readline()[:-5])
        file.readline()
        bestSolutionG.append([])
        meanSolutionG.append([])
        for j in range(iterationNo):
            line = file.readline().split()
            bestSolutionG[i].append((float(line[0]) - optimal[i]) / optimal[i])
            meanSolutionG[i].append((float(line[1]) - optimal[i]) / optimal[i])
        
        file.readline()
        bestSolutionS.append([])
        meanSolutionS.append([])
        for j in range(iterationNo):
            line = file.readline().split()
            bestSolutionS[i].append((float(line[0]) - optimal[i]) / optimal[i])
            meanSolutionS[i].append((float(line[1]) - optimal[i]) / optimal[i])
    file.close()
    
    x = [i for i in range(iterationNo)]
    for i in range(instanceNo):
        plt.figure()
        plt.plot(x, bestSolutionG[i], label = 'G')
        plt.plot(x, bestSolutionS[i], label = 'S')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title('Best ' + instanceNames[i])
        plt.xlabel('iteration')
        plt.ylabel('solution quality')
        plt.savefig('charts/4best_' + instanceNames[i] + '.pdf', bbox_inches='tight')
        
        plt.figure()
        plt.plot(x, meanSolutionG[i], label = 'G')
        plt.plot(x, meanSolutionS[i], label = 'S')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title('Mean ' + instanceNames[i])
        plt.xlabel('iteration')
        plt.ylabel('solution quality')
        plt.savefig('charts/4mean_' + instanceNames[i] + '.pdf', bbox_inches='tight')
        
def readOpimal(instanceNames, dim):
    j = 0
    permutation = []
    for name in instanceNames:
        file = open('optimal/' + name + '.opt.tour', "r")
        permutation.append([])
        for i in range(dim[j]):
            permutation[j].append(int(file.readline()))
        file.close()
        j += 1
    return permutation

def similarity(optimal, algorithm, instanceNames, iterationNo):
    similar = []
    for i in range(len(instanceNames)):
        similar.append([])
        for k in range(len(algorithm[i])):
            counter = 0
            for j in range(len(optimal[i])):
                x = optimal[i][j]
                if j < len(optimal[i]) - 1:
                    y = optimal[i][j+1]
                else:
                    y = optimal[i][0]
                idx1 = algorithm[i][k].index(x)
                idx2 = algorithm[i][k].index(y)
                if abs(idx1 - idx2) == 1 or abs(idx1 - idx2) == len(optimal[i]):
                    counter += 1
            similar[i].append(counter / len(optimal[i]))
    return similar
        

def result5(optimal):
    file = open("results5.txt", "r")
    instanceNo = int(file.readline())
    iterationNo = int(file.readline())
    instanceNames = []
    algorithmsNames = ['G', 'S']
    solutionG = []
    solutionS = []
    permutationG = []
    permutationS = []
    dim = []
    l = 0
    for i in range(instanceNo):
        name = file.readline()[:-5]
        if name != 'pr1002' and name != 'rd100':
            instanceNames.append(name)
            file.readline()
            solutionG.append([])
            permutationG.append([])
            for j in range(iterationNo):
                solutionG[i-l].append((float(file.readline()) - optimal[i]) / optimal[i])
                permutationG[i-l].append([int(x) for x in file.readline().split()])
            
            file.readline()
            solutionS.append([])
            permutationS.append([])
            for j in range(iterationNo):
                solutionS[i-l].append((float(file.readline()) - optimal[i]) / optimal[i])
                permutationS[i-l].append([int(x) for x in file.readline().split()])
            dim.append(len(permutationS[i-l][0]))
        else:
            l += 1
            for j in range(iterationNo * 4 + 2):
                file.readline()
        
    file.close()
    permutationO = readOpimal(instanceNames, dim)
    similarG = similarity(permutationO, permutationG, instanceNames, iterationNo)
    similarS = similarity(permutationO, permutationS, instanceNames, iterationNo)
    #print(solutionG)
    #print(similarS)
    
    for i in range(instanceNo - l):
        plt.figure()
        plt.scatter(solutionG[i], similarG[i], label = 'G')
        plt.scatter(solutionS[i], similarS[i], label = 'S')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title(instanceNames[i])
        plt.xlabel('solution quality')
        plt.ylabel('solution similarity')
        plt.savefig('charts/5' + instanceNames[i] + '.pdf', bbox_inches='tight')
        
def startH(optimal, dataBest, dataMean, dataDeviation, dataTimes):
    print(dataDeviation)
    file = open("resultsH.txt", "r")
    instanceNo = int(file.readline())
    instanceNames = []
    algorithmsNames = []
    results = []
    bestResults = []
    meanResults = []
    deviationResults = []
    times = []
    for i in range(instanceNo):
        instanceNames.append(file.readline()[:-5])
        algorithmsNames.append([])
        results.append([])
        bestResults.append([])
        times.append([])
        meanResults.append([])
        deviationResults.append([])
        for j in range(2):
            algorithmsNames[i].append(file.readline())
            results[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            bestResults[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            times[i].append(float(file.readline()))
        
        
            #oblicza srednia i odchylenie z rozwiazan
            mean = sum(results[i][j]) / 10
            deviation = 0
            for x in results[i][j]:
                deviation += (x - mean) * (x - mean)
            deviationResults[i].append(m.sqrt(deviation / 10))
            meanResults[i].append(mean)
    
    file.close()
    
    dataHBest = []
    dataHMean = []
    dataHDeviation = []
    dataHTimes = []
    for i in range(2):
        dataHBest.append([])
        dataHMean.append([])
        dataHDeviation.append([])
        dataHTimes.append([])
        for j in range(len(bestResults)):
            dataHBest[i].append(bestResults[j][i][-1])
            dataHMean[i].append(meanResults[j][i])
            dataHDeviation[i].append(deviationResults[j][i])
            dataHTimes[i].append(times[j][i])
    
    print(dataHDeviation)
    plt.figure()
    for i in range(2):
        plt.plot(dataHBest[i], 'o', label = algorithmsNames[0][i][:-1] + '-H start') 
    for i in range(2):
        plt.plot(dataBest[i + 3], 'o', label = algorithmsNames[0][i][:-1] + '-R start')
    plt.title('Best results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    #plt.yscale('log')
    plt.savefig('charts/Hbest.pdf', bbox_inches='tight')
    
    plt.figure()
    x = [i for i in range(instanceNo)]
    for i in range(2):
        plt.errorbar(x, dataHMean[i], yerr=dataHDeviation[i], fmt='o', label = algorithmsNames[0][i][:-1] + '-H start')
    for i in range(2):
        plt.errorbar(x, dataMean[i + 3], yerr=dataDeviation[i + 3], fmt='o', label = algorithmsNames[0][i][:-1] + '-R start') 
    plt.xticks(range(instanceNo), instanceNames)
    plt.title('Mean results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    #plt.yscale('log')
    plt.savefig('charts/Hmean.pdf', bbox_inches='tight')
    
    plt.figure()
    for i in range(2):
        plt.plot(dataHTimes[i], 'o', label = algorithmsNames[0][i][:-1] + '-H start') 
    for i in range(2):
        plt.plot(dataTimes[i + 3], 'o', label = algorithmsNames[0][i][:-1] + '-R start') 
    plt.title('Mean time of 1 run')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/Htime.pdf', bbox_inches='tight')


def compareSwap(optimal, dataBest, dataMean, dataDeviation, dataTimes):
    file = open("resultsSwap.txt", "r")
    instanceNo = int(file.readline())
    instanceNames = []
    algorithmsNames = []
    results = []
    bestResults = []
    meanResults = []
    deviationResults = []
    times = []
    for i in range(instanceNo):
        instanceNames.append(file.readline()[:-5])
        algorithmsNames.append([])
        results.append([])
        bestResults.append([])
        times.append([])
        meanResults.append([])
        deviationResults.append([])
        for j in range(2):
            algorithmsNames[i].append(file.readline())
            results[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            bestResults[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            times[i].append(float(file.readline()))
        
        
            #oblicza srednia i odchylenie z rozwiazan
            mean = sum(results[i][j]) / 10
            deviation = 0
            for x in results[i][j]:
                deviation += (x - mean) * (x - mean)
            deviationResults[i].append(m.sqrt(deviation / 10))
            meanResults[i].append(mean)
    
    file.close()
    
    dataHBest = []
    dataHMean = []
    dataHDeviation = []
    dataHTimes = []
    for i in range(2):
        dataHBest.append([])
        dataHMean.append([])
        dataHDeviation.append([])
        dataHTimes.append([])
        for j in range(len(bestResults)):
            dataHBest[i].append(bestResults[j][i][-1])
            dataHMean[i].append(meanResults[j][i])
            dataHDeviation[i].append(deviationResults[j][i])
            dataHTimes[i].append(times[j][i])
    
    plt.figure()
    for i in range(2):
        plt.plot(dataHBest[i], 'o', label = algorithmsNames[0][i][:-1] + '-swap elem') 
    for i in range(2):
        plt.plot(dataBest[i + 3], 'o', label = algorithmsNames[0][i][:-1] + '-arc inversion')
    plt.title('Best results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/SwapBest.pdf', bbox_inches='tight')
    
    plt.figure()
    x = [i for i in range(instanceNo)]
    for i in range(2):
        plt.errorbar(x, dataHMean[i], yerr=dataHDeviation[i], fmt='o', label = algorithmsNames[0][i][:-1] + '-swap elem')
    for i in range(2):
        plt.errorbar(x, dataMean[i + 3], yerr=dataDeviation[i + 3], fmt='o', label = algorithmsNames[0][i][:-1] + '-arc inversion') 
    plt.xticks(range(instanceNo), instanceNames)
    plt.title('Mean results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.yscale('log')
    plt.savefig('charts/SwapMean.pdf', bbox_inches='tight')
    
    plt.figure()
    for i in range(2):
        plt.plot(dataHTimes[i], 'o', label = algorithmsNames[0][i][:-1] + '-swap elem') 
    for i in range(2):
        plt.plot(dataTimes[i + 3], 'o', label = algorithmsNames[0][i][:-1] + '-arc inversion') 
    plt.title('Mean time of 1 run')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(instanceNo), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/SwapTime.pdf', bbox_inches='tight')
    

optimal, dataBest, dataMean, dataDeviation, dataTimes = result2()
#result3(optimal)
#result4(optimal)
#result5(optimal)
#startH(optimal, dataBest, dataMean, dataDeviation, dataTimes)
#compareSwap(optimal, dataBest, dataMean, dataDeviation, dataTimes)