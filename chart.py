import matplotlib.pyplot as plt
import math as m

def readResult2():
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
            solutionNo[i].append([int(x) for x in file.readline().split()])
            steps[i].append([int(x) for x in file.readline().split()])
            times[i].append(float(file.readline()))
            iterationTimes[i].append(float(file.readline()))
        
            #oblicza srednia i odchylenie z rozwiazan
            meanResults[i].append([])
            deviationResults[i].append([])
            sumResults = 0
            for l in range(1, 11):
                sumResults += results[i][j][l - 1]
                meanResults[i][j].append(sumResults / l)
                deviation = 0
                for k in range(l):
                    deviation += (results[i][j][k] - meanResults[i][j][l - 1]) * (results[i][j][k] - meanResults[i][j][l - 1])
                deviationResults[i][j].append(m.sqrt(deviation / l))
                
            
        #stara wersja - osobne wykresy dla kazdej instancji
        """plt.figure()
        xAxis = [i for i in range(1, 11)]
        for j in range(5):
            plt.plot(xAxis, bestResults[j], label = algorithmsNames[j][:-1] + ' best')
            plt.errorbar(xAxis, meanResults[j], yerr=deviationResults[j], fmt='-o', label = algorithmsNames[j][:-1] + ' mean')
            #plt.plot(xAxis, worstResults[j], label = algorithmsNames[j][:-1] + ' worst')
        plt.title(instance + ' Compare with optimum')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylim(0, 1)
        plt.xlim(1, 10)
        #plt.show()
        plt.savefig('charts/' + instance + '_result.pdf', bbox_inches='tight')
        
        #Porównanie algorytmów H, S, G
        alg = [0, 3, 4]
        plt.figure()
        xAxis = [i for i in range(1, 11)]
        for j in alg:
            plt.plot(xAxis, bestResults[j], 'o-', label = algorithmsNames[j][:-1] + ' best')
            plt.errorbar(xAxis, meanResults[j], yerr=deviationResults[j], fmt='-o', label = algorithmsNames[j][:-1] + ' mean')
            plt.plot(xAxis, worstResults[j], 'o-', label = algorithmsNames[j][:-1] + ' worst')
        plt.title(instance + ' Compare with optimum')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xlim(1, 10)
        #plt.show()
        plt.savefig('charts/' + instance + '_G_S_H_result.pdf', bbox_inches='tight')
        
        #Porównanie algorytmów R, RW
        alg = [1, 2]
        plt.figure()
        xAxis = [i for i in range(1, 11)]
        for j in alg:
            plt.plot(xAxis, bestResults[j], 'o-', label = algorithmsNames[j][:-1] + ' best')
            plt.errorbar(xAxis, meanResults[j], yerr=deviationResults[j], fmt='-o', label = algorithmsNames[j][:-1] + ' mean')
            plt.plot(xAxis, worstResults[j], 'o-', label = algorithmsNames[j][:-1] + ' worst')
        plt.title(instance + ' Compare with optimum')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xlim(1, 10)
        #plt.show()
        plt.savefig('charts/' + instance + '_R_RW_result.pdf', bbox_inches='tight')
        
        #Liczba krokow G i S
        plt.figure()
        xAxis = [i for i in range(1, 11)]
        for j in range(3, 5):
            plt.plot(xAxis, steps[j], label = algorithmsNames[j][:-1])
        plt.title(instance + ' Number of steps')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xlim(1, 10)
        #plt.show()
        plt.savefig('charts/' + instance + '_steps.pdf', bbox_inches='tight')
        
        #Liczba przejrzanyc rozwiazan
        plt.figure()
        xAxis = [i for i in range(1, 11)]
        for j in range(5):
            plt.plot(xAxis, solutionNo[j], 'o', label = algorithmsNames[j][:-1])
        plt.title(instance + ' Number of solutions')
        plt.legend(prop={'size': 10}, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xlim(1, 10)
        #plt.show()
        plt.savefig('charts/' + instance + '_solution_No.pdf', bbox_inches='tight')"""
    
    file.close()
    
def readResult3():
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
            startSolutionG[i].append(float(line[0]))
            finishSolutionG[i].append(float(line[1]))
        
        file.readline()
        startSolutionS.append([])
        finishSolutionS.append([])
        for j in range(iterationNo):
            line = file.readline().split()
            startSolutionS[i].append(float(line[0]))
            finishSolutionS[i].append(float(line[1]))
    file.close()
    
def readResult4():
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
            bestSolutionG[i].append(float(line[0]))
            meanSolutionG[i].append(float(line[1]))
        
        file.readline()
        bestSolutionS.append([])
        meanSolutionS.append([])
        for j in range(iterationNo):
            line = file.readline().split()
            bestSolutionS[i].append(float(line[0]))
            meanSolutionS[i].append(float(line[1]))
    file.close()

readResult2()
readResult3()
readResult4()