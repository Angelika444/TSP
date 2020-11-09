import matplotlib.pyplot as plt
import math as m

#wczytanie danych
file = open("results.txt", "r")
instanceNo = int(file.readline())
for i in range(instanceNo):
    algorithmsNames = []
    results = []
    bestResults = []
    worstResults = []
    meanResults = []
    deviationResults = []
    solutionNo = []
    steps = []
    instance = file.readline()[:-5]
    optimal = float(file.readline())
    for j in range(5):
        algorithmsNames.append(file.readline())
        results.append([(float(x) - optimal) / float(x) for x in file.readline().split()])
        bestResults.append([(float(x) - optimal) / float(x) for x in file.readline().split()])
        worstResults.append([(float(x) - optimal) / float(x) for x in file.readline().split()])
        solutionNo.append([int(x) for x in file.readline().split()])
        steps.append([int(x) for x in file.readline().split()])
    
        #oblicza srednia i odchylenie z rozwiazan
        meanResults.append([])
        deviationResults.append([])
        sumResults = 0
        for l in range(1, 11):
            sumResults += results[j][l - 1]
            meanResults[j].append(sumResults / l)
            deviation = 0
            for k in range(l):
                deviation += (results[j][k] - meanResults[j][l - 1]) * (results[j][k] - meanResults[j][l - 1])
            deviationResults[j].append(m.sqrt(deviation / l))
        
    #rysowanie wykresow - osobno dla kazdej instancji
    #Miara jakoci rozwiązania - ile razy gorsze od optimum
    #Porownanie wszystkich algorytmow
    plt.figure()
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
    plt.savefig('charts/' + instance + '_solution_No.pdf', bbox_inches='tight')
    
file.close()
