import matplotlib.pyplot as plt
import math as m

def result2():
    file = open("results2.txt", "r")
    instanceNo = int(file.readline())
    instanceNames = []
    optimal = []
    algorithmsNames = []
    results = []
    bestResults = []
    worstResults = []
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
        worstResults.append([])
        solutionNo.append([])
        steps.append([])
        times.append([])
        iterationTimes.append([])
        meanResults.append([])
        deviationResults.append([])
        for j in range(5):
            algorithmsNames[i].append(file.readline())
            results[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            bestResults[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
            worstResults[i].append([(float(x) - optimal[i]) / optimal[i] for x in file.readline().split()])
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
    dataWorst = []
    dataDeviation = []
    dataSolutionNo = []
    dataSteps = []
    dataTimes = []
    dataIterationTimes = []
    for i in range(5):
        dataBest.append([])
        dataMean.append([])
        dataDeviation.append([])
        dataWorst.append([])
        dataSolutionNo.append([])
        dataSteps.append([])
        dataTimes.append([])
        dataIterationTimes.append([])
        for j in range(len(bestResults)):
            dataBest[i].append(bestResults[j][i][-1])
            dataMean[i].append(meanResults[j][i])
            dataDeviation[i].append(deviationResults[j][i])
            dataWorst[i].append(worstResults[j][i][-1])
            dataSolutionNo[i].append(solutionNo[j][i])
            dataSteps[i].append(steps[j][i])
            dataTimes[i].append(times[j][i])
            dataIterationTimes[i].append(iterationTimes[j][i])
            
    
    plt.figure()
    for i in range(5):
        plt.plot(dataBest[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Best results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.xticks(range(9), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2best.pdf', bbox_inches='tight')
    #plt.show()
    

    plt.figure()
    x = [i for i in range(9)]
    for i in range(5):
        plt.errorbar(x, dataMean[i], yerr=dataDeviation[i], fmt='o', label = algorithmsNames[0][i]) 
    plt.xticks(range(9), instanceNames)
    plt.title('Mean results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.yscale('log')
    plt.savefig('charts/2mean.pdf', bbox_inches='tight')
    #plt.show()
    
    plt.figure()
    for i in range(5):
        plt.plot(dataWorst[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Worst results')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    plt.ylabel('solution quality')
    plt.xticks(range(9), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2worst.pdf', bbox_inches='tight')
    #plt.show()
    
    plt.figure()
    for i in range(5):
        plt.plot(dataSolutionNo[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Solution number')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(9), instanceNames)
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
    plt.xticks(range(9), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2steps.pdf', bbox_inches='tight')
    #plt.show()
    
    plt.figure()
    for i in range(5):
        plt.plot(dataTimes[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Mean time of 1 run')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(9), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2time.pdf', bbox_inches='tight')
    #plt.show()
    
    plt.figure()
    for i in range(5):
        plt.plot(dataIterationTimes[i], 'o', label = algorithmsNames[0][i]) 
    plt.title('Mean iteration time')
    plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('instance')
    #plt.ylabel('solution quality')
    plt.xticks(range(9), instanceNames)
    plt.yscale('log')
    plt.savefig('charts/2iterationTime.pdf', bbox_inches='tight')
    #plt.show()
    
    return optimal
    
    
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

optimal = result2()
result3(optimal)
result4(optimal)