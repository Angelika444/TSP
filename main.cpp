#include<iostream>
#include <vector>
#include<cstdlib>
#include<ctime>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
using namespace std;


struct point{
    int index;
    float x;
    float y;
};

//liczba zapamitanych rozwiazan dla p. 3
const int iterations_p3 = 200;
//liczba zapamitanych rozwiazan dla p. 4
const int iterations_p4 = 300;
//liczba zapamitanych rozwiazan dla p. 5
const int iterations_p5 = 100;

struct algorithmResult {
    double best = 1000000000;
    double worst = -1;
    double sum = 0;
    //liczba krokow - dla greedy i steepest
    int steps = 0;
    //liczba przejrzanych rozwiazan
    long long int solutionNo = 0;
    //aktualny sumaryczny dystans
    double actual;
    int iterations = 0;
    double time;
    double iterationTime;
    string name;
    double firstBestResults[10];
    double firstWorstResults[10];
    double first10Results[10];
    int firstSteps[10] = {0};
    long long int firstSolutionNo[10];
    double startSolution[iterations_p3];
    double finishSolution[iterations_p3];
    double bestSolution[iterations_p4];
    double meanSolution[iterations_p4];
};

int extractIntegerWords(string str)
{
    stringstream ss;
    /* Storing the whole string into string stream */
    ss << str;
    /* Running loop till the end of the stream */
    string temp;
    int found;
    while (!ss.eof()) {
        /* extracting word by word from stream */
        ss >> temp;
        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found)
            //cout << found << " ";
            ;
        /* To save from space at the end of string */
        temp = "";
    }
    return found;
}

double getDistance(int i, int j, double** distanceMatrix)
{
    return distanceMatrix[i][j];
}

double calcDistance(point p1, point p2)
{
    //dystans euklidesowy (sa tez inne)
    double x1 = p1.x - p2.x, x2 = p1.y - p2.y;
    return sqrt(x1 * x1 + x2 * x2);
}

double createDistanceMatrix(int dim, double** matrix, point tab[])
{
    double mini = 1000000000000;
    for(int i = 0; i<dim-1; i++)
    {
        matrix[i][i] = 0;
        for(int j = i+1; j<dim; j++)
        {
            matrix[i][j] = calcDistance(tab[i], tab[j]);
            matrix[j][i] = matrix[i][j];
            if (mini > matrix[i][j]) mini = matrix[i][j];
        }
    }
    return mini;
}

double calcSolutionDistance(int dim, int solution[], double** distanceMatrix)
{
    double distance = getDistance(solution[0], solution[dim - 1], distanceMatrix);
    for (int i = 0; i < dim - 1; i++)
    {
        distance += getDistance(solution[i], solution[i + 1], distanceMatrix);
    }
    return distance;
}

void initResult(algorithmResult &result, int dim, int solution[], double** distanceMatrix, string name)
{
    double distance = calcSolutionDistance(dim, solution, distanceMatrix);
    result.best = distance;
    result.worst = 0;
    result.actual = distance;
    result.name = name;
}

void actualizeResult(algorithmResult &result, int dim, int solution[], int counter, int steps)
{
    result.sum += result.actual;
    if (result.actual < result.best) result.best = result.actual;
    if (result.actual > result.worst) result.worst = result.actual;
    result.solutionNo += counter;
    result.steps += steps;
}

void randomSolution (int dim, int solution[])
{
    for (int i = dim; i > 0; i--)
    {
        int index = rand() % i;
        swap(solution[index], solution[i - 1]);
    }
}

void HSolution(point tab[], int dim, int solution[], double** distanceMatrix, vector < int > &startCity, algorithmResult &result)
{
    //solution is a permutation with cities index in tab, not field index in struct, so from 0 to dim-1
    int i, j, solution_index, index_min_dinstance, startCityIndex;
    double distance = 0, current_distance, min_distance;
    //startCity to wektor z indeksami miast, usuwamy miasto, ktore zostalo wylosowane jako poczatkowe
    //w ten sposob przy kazdym odpaleniu algorytmu startujemy z innego miasta
    startCityIndex = rand() % startCity.size();
    solution_index = startCity[startCityIndex];
    startCity.erase(startCity.begin() + startCityIndex);
    solution[0] = solution_index;
    //wektor nie odwiedzonych jeszcze miast, odwiedzone miasta usuwamy, dzieki czemu liczymy dystans tylko do tych, ktore mozemy jeszcze odwiedzic
    vector <int> not_visited;
    for (i = 0; i < dim; i++)
    {
        not_visited.push_back(i);
    }
    not_visited.erase(not_visited.begin() + solution_index);
    for (i = 1; i < dim ; i++)
    {
        index_min_dinstance = 0;
        min_distance = getDistance(solution_index, not_visited[index_min_dinstance], distanceMatrix);
        for (j = 1; j < not_visited.size(); j++)
        {
            current_distance = getDistance(solution_index, not_visited[j], distanceMatrix);
            if (current_distance < min_distance) {
                min_distance = current_distance;
                index_min_dinstance = j;
            }
        }
        solution[i] = not_visited[index_min_dinstance];
        solution_index = not_visited[index_min_dinstance];
        not_visited.erase(not_visited.begin() + index_min_dinstance);
        distance += min_distance;
    }
    distance += getDistance(solution[0], solution[dim-1], distanceMatrix);
    result.actual = distance;
    actualizeResult(result, dim, solution, 1, 1);

}

double calcDistDif(int solution[], int i, int j, int dim, double** distanceMatrix)
{
    //obliczam roznice w jakosci rozwiazania zamiast liczyc wszystko od nowa
    double sol_diff = 0.0;
    //odwrocenie kolejnosci luku - odejmuemy 2 stare drogi i dodajemy 2 nowe drogi
    //sprawdzamy przy tym czy nie trafilismy na skrajne wartosci i=0, j=dim-1
    if (i != 0)
    {
        sol_diff -= getDistance(solution[i-1], solution[i], distanceMatrix);
        sol_diff += getDistance(solution[i-1], solution[j], distanceMatrix);
    }
    else
    {
        sol_diff -= getDistance(solution[dim-1], solution[i], distanceMatrix);
        sol_diff += getDistance(solution[dim-1], solution[j], distanceMatrix);
    }

    if (j != dim-1)
    {
        sol_diff -= getDistance(solution[j], solution[j+1], distanceMatrix);
        sol_diff += getDistance(solution[i], solution[j+1], distanceMatrix);
    }
    else if (i != 0)
    {
        sol_diff -= getDistance(solution[j], solution[0], distanceMatrix);
        sol_diff += getDistance(solution[i], solution[0], distanceMatrix);
    }
    //kiedy i=0 i j=dim-1
    else
    {
        sol_diff = 0;
        //sol_diff -= getDistance(solution[j], solution[j-1], distanceMatrix);
        //sol_diff += getDistance(solution[i], solution[j-1], distanceMatrix);
        //sol_diff -= getDistance(solution[i], solution[i+1], distanceMatrix);
        //sol_diff += getDistance(solution[j], solution[i+1], distanceMatrix);
    }

    return sol_diff;
}

double calcDistDifSwap(int solution[], int i, int j, int dim, double** distanceMatrix)
{
    //obliczam roznice w jakosci rozwiazania zamiast liczyc wszystko od nowa
    double sol_diff = 0.0;
    //kiedy i=0 i j=dim-1
    if (i == 0 && j == dim - 1)
    {
        sol_diff = 0;
        sol_diff -= getDistance(solution[j], solution[j-1], distanceMatrix);
        sol_diff += getDistance(solution[i], solution[j-1], distanceMatrix);
        sol_diff -= getDistance(solution[i], solution[i+1], distanceMatrix);
        sol_diff += getDistance(solution[j], solution[i+1], distanceMatrix);
    }
    else if (j - i == 1)
    {
        if (i != 0)
        {
            sol_diff -= getDistance(solution[i-1], solution[i], distanceMatrix);
            sol_diff += getDistance(solution[i-1], solution[j], distanceMatrix);
        }
        else
        {
            sol_diff -= getDistance(solution[dim-1], solution[i], distanceMatrix);
            sol_diff += getDistance(solution[dim-1], solution[j], distanceMatrix);
        }

        if (j != dim-1)
        {
            sol_diff -= getDistance(solution[j], solution[j+1], distanceMatrix);
            sol_diff += getDistance(solution[i], solution[j+1], distanceMatrix);
        }
        else
        {
            sol_diff -= getDistance(solution[j], solution[0], distanceMatrix);
            sol_diff += getDistance(solution[i], solution[0], distanceMatrix);
        }
    }
    else
    {
        if (i != 0)
        {
            sol_diff -= getDistance(solution[i-1], solution[i], distanceMatrix);
            sol_diff += getDistance(solution[i-1], solution[j], distanceMatrix);
        }
        else
        {
            sol_diff -= getDistance(solution[dim-1], solution[i], distanceMatrix);
            sol_diff += getDistance(solution[dim-1], solution[j], distanceMatrix);
        }
        sol_diff -= getDistance(solution[i+1], solution[i], distanceMatrix);
        sol_diff += getDistance(solution[i+1], solution[j], distanceMatrix);

        if (j != dim-1)
        {
            sol_diff -= getDistance(solution[j], solution[j+1], distanceMatrix);
            sol_diff += getDistance(solution[i], solution[j+1], distanceMatrix);
        }
        else
        {
            sol_diff -= getDistance(solution[j], solution[0], distanceMatrix);
            sol_diff += getDistance(solution[i], solution[0], distanceMatrix);
        }
        sol_diff -= getDistance(solution[j], solution[j-1], distanceMatrix);
        sol_diff += getDistance(solution[i], solution[j-1], distanceMatrix);
    }


    return sol_diff;
}

//zamiana wszystkich elementow w luku - ostatni z pierwszym az nie dojdziemy do srodka
void swapElements(int i, int j, int dim, int solution[])
{
    //jesli i=0, j=dim-1, to zamieniamy je miejscami, bo odwrocenie wszystkich elementow pomiedzy nimi daloby dokladnie to samo sasiedztwo
    //if (i == 0 && j == dim - 1) swap(solution[i], solution[j]);
    //else
    //{
       for (int k = 0; k < ceil(double(j + i) / 2) - i; k++)
        {
            swap(solution[i + k], solution[j - k]);
        }
    //}
}

void greedySwap(int dim, int solution[], double** distanceMatrix, algorithmResult &result, time_t time_start, int algorithmTime)
{
    //greedy
    bool swaped;
    time_t current_time;
    int counter = 1, steps = 0;
    if (result.iterations < iterations_p3) result.startSolution[result.iterations] = result.actual;
    do
    {
        swaped = false;
        //generuj po kolei elementy sasiedztwa
        // zewnetrzna petla musi byc do dim-1, a wewnetrzna do dim, inaczej nie bedziemy miec nigdy pary 0, dim-1
        for(int i = 0; i < dim - 1; i++)
        {
            for(int j = i+1; j < dim; j++)
            {
                counter++;
                //double sol_diff = calcDistDif(solution, i, j, dim, distanceMatrix);
                //przy zamianie elementow
                double sol_diff = calcDistDifSwap(solution, i, j, dim, distanceMatrix);

                //jesli uzyskalismy poprawe, to nalezy dokonac zamiany od razu
                if (sol_diff < 0)
                {
                    //odwrocenie luku - trzeba odwrocic wszyskie elementy pomiedzy i a j
                    //swapElements(i, j, dim, solution);
                    //przy zamianie par - w elementach nieobowiazkowych jest, zeby porownac te 2 sasiedztwa
                    swap(solution[i], solution[j]);

                    swaped = true;
                    steps++;
                    //aktualizujemy aktualna dlugosc sciezki, zeby potem nie trzeba jej bylo ponownie obliczac
                    result.actual += sol_diff;
                    i = dim;
                    break;
                }
            }
        }
        current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
    } while (/*current_time < algorithmTime && */swaped);

    //wynik koncowy jest najlepszym dotychczasowym, dlatego aktualizujemy go na koncu
    actualizeResult(result, dim, solution, counter, steps);
    if (result.iterations < 10)
    {
        result.firstBestResults[result.iterations] = result.best;
        //result.firstWorstResults[result.iterations] = result.worst;
        result.first10Results[result.iterations] = result.actual;
        result.firstSteps[result.iterations] = steps;
        result.firstSolutionNo[result.iterations] = counter;
    }
    if (result.iterations < iterations_p3) result.finishSolution[result.iterations] = result.actual;
    if (result.iterations < iterations_p4) result.bestSolution[result.iterations] = result.best;
    if (result.iterations < iterations_p4) result.meanSolution[result.iterations] = result.sum / (result.iterations + 1);

}

void steepestSwap(point tab[], int dim, int solution[], double** distanceMatrix, algorithmResult &result, time_t time_start, int algorithmTime)
{
    //steepest
    bool swaped = true;
    time_t current_time;
    int counter = 0, numIterations = dim * (dim - 1) / 2 + 1, steps = 0;
    double best_diff = 0.0;
    double sol_diff = 0.0;
    int best_i=0, best_j=0;
    if (result.iterations < iterations_p3) result.startSolution[result.iterations] = result.actual;

    do
    {
        counter += numIterations;
        best_diff = 0.0;
        //generuj po kolei elementy sasiedztwa
        for(int i = 0; i < dim - 1; i++)
        {
            for(int j = i+1; j < dim; j++)
            {
                //double sol_diff = calcDistDif(solution, i, j, dim, distanceMatrix);
                //przy zamianie elementow
                double sol_diff = calcDistDifSwap(solution, i, j, dim, distanceMatrix);

                //jesli uzyskalismy poprawe, to nalezy zapisac jej wartosc do porownania z innymi
                if(sol_diff < best_diff)
                {
                    best_diff = sol_diff;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        //dopiero teraz po przejrzeniu ca³ego sasiedztwa robimy swapa
        if(best_diff < 0)   //zabezpieczam sie przed wywaleniem algorytmu, jezeli bylibysmy ju¿ w optimum, i nie ma lepszego rozwiazania
        {
            //odwrocenie luku - trzeba odwrocic wszyskie elementy pomiedzy i a j,
            //swapElements(best_i, best_j, dim, solution);
            //przy zamianie par - w elementach nieobowiazkowych jest, zeby porownac te 2 sasiedztwa
            swap(solution[best_i], solution[best_j]);

            steps++;
            result.actual += best_diff;
        }
        else
        {
            swaped = false;
        }
        current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
    }
    while (/*current_time < algorithmTime && */swaped);

    //wynik koncowy jest najlepszym dotychczasowym, dlatego aktualizujemy go na koncu
    actualizeResult(result, dim, solution, counter, steps);
    if (result.iterations < 10)
    {
        result.firstBestResults[result.iterations] = result.best;
        //result.firstWorstResults[result.iterations] = result.worst;
        result.first10Results[result.iterations] = result.actual;
        result.firstSteps[result.iterations] = steps;
        result.firstSolutionNo[result.iterations] = counter;
    }
    if (result.iterations < iterations_p3) result.finishSolution[result.iterations] = result.actual;
    if (result.iterations < iterations_p4) result.bestSolution[result.iterations] = result.best;
    if (result.iterations < iterations_p4) result.meanSolution[result.iterations] = result.sum / (result.iterations + 1);

}

void tabuSearch(point tab[], int dim, int solution[], double** distanceMatrix, algorithmResult &result, time_t time_start, int algorithmTime)
{
    //tabu search
    //rozmiary zgodne z tym, co jest na stronie
    int tabuListSize = int (dim / 4), masterListSize = int (dim / 10);
    //pomocnicze struktury danych
    //do masterList, masterListIndex zapisujemy roznice w rozwiazaniach i odpowiadajacy indeks zamiany
    //potem wybieramy masterListSize najlepszych rozwiazan i je przerzucamy do masterlistChosen i masterListIndexChosen
    int tabuList[tabuListSize][2], tabuListIndex = 0, i, j, ii, s, best_index, getBest = 0;
    vector <double> masterList, masterListChosen;
    vector <int> masterListIndexI, masterListIndexJ, masterListIndexIChosen, masterListIndexJChosen;
    double bestRejected, worstElem, actualBest = result.actual;
    bool tabu;
    //limit okresla ile iteracji wykonujemy bez poprawy najlepszego rozwiazania, dla >1000 jest 100, bo to i tak dlugo dziala
    int limit = 300;
    if (dim > 1000) limit = 100;
    else if (dim <= 100) limit = 500;

    //liczba przejrzanych rozwiazan przy przeskaniu wszystkich par i,j
    int numIterations = dim * (dim - 1) / 2 + 1;
    double best_diff = 0.0;
    double sol_diff = 0.0;
    int best_i = 0, best_j = 0, last_i = -1, last_j = -1;

    do
    {
        masterList.clear();
        masterListChosen.clear();
        masterListIndexI.clear();
        masterListIndexJ.clear();
        masterListIndexIChosen.clear();
        masterListIndexJChosen.clear();
        //generuj po kolei elementy sasiedztwa
        //tak jak dla steepest
        for(i = 0; i < dim - 1; i++)
        {
            for(j = i+1; j < dim; j++)
            {
                tabu = false;
                //sol_diff = calcDistDif(solution, i, j, dim, distanceMatrix);
                //przy zamianie elementow
                sol_diff = calcDistDifSwap(solution, i, j, dim, distanceMatrix);

                //jesli rozwiazanie jest lepsze od najlepszego, to łamiemy liste tabu
                if (result.actual + sol_diff < actualBest - 0.000001)
                {
                    masterList.push_back(sol_diff);
                    masterListIndexI.push_back(i);
                    masterListIndexJ.push_back(j);
                }
                else
                {
                    //sprawdzamy czy ruch jest zgodny z lista tabu
                    for(ii = 0; ii < tabuListSize; ii++)
                    {
                        if((tabuList[ii][0] == solution[i] && tabuList[ii][1] == solution[j]) || (tabuList[ii][0] == solution[j] && tabuList[ii][1] == solution[i])) tabu = true;
                    }
                    //jesli zamiana jest zgodna z tabu
                    //i mamy mniej elementów na liscie master niz jest wymagane lub rozwiazanie jest lepsze od najgorszego z tabu (zeby nie wpisywac wszystkich rozwiazan na liste tabu)
                    if (!tabu && (masterList.size() < masterListSize || sol_diff < worstElem))
                    {
                        masterList.push_back(sol_diff);
                        masterListIndexI.push_back(i);
                        masterListIndexJ.push_back(j);
                        if (sol_diff > worstElem) worstElem = sol_diff;
                    }
                }
            }
        }

        //masterListSize najlepszych rozwiazan przepisujemy na liste masterListChosen
        for (j = 0; j < min (masterListSize, int(masterList.size())); j++)
        {
            best_diff = masterList[0];
            best_index = 0;
            for (i = 1; i < masterList.size(); i++)
            {
                if (masterList[i] < best_diff)
                {
                    best_diff = masterList[i];
                    best_index = i;
                }
            }
            masterListChosen.push_back(best_diff);
            masterListIndexIChosen.push_back(masterListIndexI[best_index]);
            masterListIndexJChosen.push_back(masterListIndexJ[best_index]);

            masterList.erase(masterList.begin() + best_index);
            masterListIndexI.erase(masterListIndexI.begin() + best_index);
            masterListIndexJ.erase(masterListIndexJ.begin() + best_index);
        }

        //bestRejected - najlepsze rozwiazanie, ktorego nie wybralismy
        if (masterList.size() > 0)
        {
            bestRejected = masterList[0];
            for (i = 1; i < masterList.size(); i++)
            {
                if (masterList[i] < bestRejected) bestRejected = masterList[i];
            }
        }
        else bestRejected = 100000000;

        for (j = 0; j < int(masterListChosen.size()); j++)
        {
            best_i = masterListIndexIChosen[j];
            best_j = masterListIndexJChosen[j];

            //sprawdzamy, czy nie bylo zamiany w stylu (0, 10) a teraz jest (11, dim-1), bo wtedy dostalibysmy oryginalny uklad miast, przez co mozna by sie zapetlic
            if ((best_i == 0 && last_j == dim - 1 && last_i - best_j == 1) || (last_i == 0 && best_j == dim - 1 && best_i - last_j == 1));
            else
            {
                //liczymy jeszcze raz roznice, bo po pierwszej zamianie mogla zmienic sie wartosc rozwiazania dla innych par wierzcholkow z listy master
                //best_diff = calcDistDif(solution, best_i, best_j, dim, distanceMatrix);
                //przy zamianie elementow
                best_diff = calcDistDifSwap(solution, i, j, dim, distanceMatrix);

                //jezeli rozwiazanie jest gorsze od najlepszego odrzuconego, to wychodzimy z petli - tworzymy ponownie masterList
                if (best_diff > bestRejected) j = int(masterListChosen.size());
                //w przeciwnym wypadku zamieniamy wierzcholki
                else
                {
                    //odwrocenie luku - trzeba odwrocic wszyskie elementy pomiedzy i a j,
                    //swapElements(best_i, best_j, dim, solution);
                    //przy zamianie par - w elementach nieobowiazkowych jest, zeby porownac te 2 sasiedztwa
                    swap(solution[best_i], solution[best_j]);

                    //na liste tabu wpisujemy pare wierzcholkow, ktorych nie mozna zamieniac
                    //z listy tabu znikaja najstarsze elementy
                    tabuList[tabuListIndex][0] = solution[best_i];
                    tabuList[tabuListIndex][1] = solution[best_j];
                    tabuListIndex++;
                    if (tabuListIndex == tabuListSize) tabuListIndex = 0;

                    //aktualizujemy obecny wynik funcji celu
                    result.actual += best_diff;
                    //liczymy ile iteracji mija od ostatniego polepszenia wyniku
                    //actualBest jest najlepszym wynikiem osiagnietym podczas tego uruchomienia algorytmu
                    getBest++;
                    if (result.actual < actualBest - 0.000001)
                    {
                        //cout<<getBest<<" "<<actualBest<<endl;
                        actualBest = result.actual;
                        getBest = 0;
                    }

                    //if j==0, to znaczy ze przejrzelismy numIterations (tyle co w steepest)
                    // else przejrzelismy 1, bo ponownie musielismy obliczyc rozwiazanie
                    if (j == 0) actualizeResult(result, dim, solution, numIterations, 1);
                    else actualizeResult(result, dim, solution, 1, 1);

                    last_i = best_i;
                    last_j = best_j;
                }
            }
        }
    }
    //jak rozwiazanie nie poprawi sie przez limit iteracji, to przerywamy algorytm
    while (getBest < limit);

    if (result.iterations < 10)
    {
        result.firstBestResults[result.iterations] = result.best;
        //nie zapisujemy najgorsze rozwiązania, bo to trzeba zrobic inaczej dla tabu, a to nie jest obowiązkowe
        //result.firstWorstResults[result.iterations] = result.worst;
        result.first10Results[result.iterations] = actualBest;
        result.firstSolutionNo[result.iterations] = result.solutionNo;
    }
}

void randomWalkSwap(point tab[], int dim, int solution[], double** distanceMatrix, algorithmResult &result, double algorithmTime)
{
    //random walk
    double sol_diff;
    int i=0, j=0;
    bool swaped = false;
    int counter = 0;
    time_t time_start = clock();
    //dostajemy pierwsze losowe rozwiazanie, ktore uwzgledniamy, bo moze sie okazac najlepsze
    actualizeResult(result, dim, solution, 1, 1);
    //idziemy dana sciezka 1000 razy, potem losujemy nowe rozwiazanie - nowa sciezka
    //while(counter < 1000)
    //idziemy tyle czasu ile trwalo uruchomienie greedy
    do
    {
        counter++;
        //losuj pary i, j, nie chcemy takich samych wartosci bo taka zamiana nic
        do{
            i = rand() % dim;
            j = rand() % dim;
        } while(i==j);
        //robie swapa jezeli i > j zeby nie sprawdzac dodatkowych warunkow brzegowych w liczeniu funkcji kosztu
        if (i > j){
            swap(i, j);
        }

        sol_diff = calcDistDif(solution, i, j, dim, distanceMatrix);
        swapElements(i, j, dim, solution);
        //swap(solution[i], solution[j]);
        result.actual += sol_diff;
        //w RW kazdy posredni wynik moze byc najlepszy, dlatego aktualizujemy wynik przy kazdym rozwiazaniu
        actualizeResult(result, dim, solution, 1, 1);
        swaped = true;
    } while((clock() - time_start) / double(CLOCKS_PER_SEC) < algorithmTime);

    if (result.iterations < 10)
    {
        result.firstBestResults[result.iterations] = result.best;
        //result.firstWorstResults[result.iterations] = result.worst;
        //jako rozwiazanie z iteracji zapisywana jest srednia rozwiazan z tej iteracji
        result.first10Results[result.iterations] = result.sum / result.steps;
        result.firstSolutionNo[result.iterations] = counter + 1;
        result.sum = 0;
        result.steps = 0;
    }
}

double temp(double prev_temp, double alfa)
{
    return alfa * prev_temp;
}

void SimulatedAnnealing(point tab[], int dim, int solution[], double** distanceMatrix, algorithmResult &result, time_t time_start, int algorithmTime, double min_diff)
{
    //deklarujemy hiperparametry

    //temperatura startowa taka, żeby akceptowane było 95% ruchów.
    //liczę to tak szacując, żę 1/2 rozwiązań jest lepsza a 1/2 rozwiązań jest gorsza
    //lepsze 0.5 będzie zawsze akceptowane
    //0.9 z pozostałych 0.5 też będzie zawsze akceptowane, czyli na początku wartość funkcji temperatury musi wynosić 0.9 (dobrze???)
    double c_start = 0.9;
    //długość łańcucha Markowa L (liczba iteracji na jednym poziomie temp.) proporcjonalna do średniego rozmiaru sąsiedztwa (chodzi o to, by dać szanse odwiedzenia każdego sąsiada), a zatem zwykle zależna od rozmiaru instancji
    //int instance_size = sizeof(solution)/sizeof(solution[0]);
    int L = dim * dim * 2;
    //funkcja dekrementacji temperatury – c(k+1) = alfa * c(k), 0.8 < alfa <0.99, niech alfa = 0.9
    //zdefiniowana wyżej - tutaj ustalam współczynnik alfa
    double alfa = 0.95;
    //warunek stopu (2 warunki):
    //brak poprawy po P*L iteracjach (np. P=10)
    int P = 10;
    //temperatura spada do poziomu na którym prawdopodobieństwo akceptacji ruchów pogarszających jest bliskie zeru, np. 0.01.
    double min_probability = 0.001;

    //inne potrzebne zmienne
    double c = c_start;
    double prob, random_number;
    int currentL = 0;
    int noIterationsNoImprovement = 0;
    bool noImprovement = false;
    double actualBest = result.actual;



    //zaczynamy jak greedy
    time_t current_time;
    int counter = 1;
    do
    {
        //generuj po kolei elementy sasiedztwa
        for(int i = 0; i < dim - 1; i++)
        {
            for(int j = i+1; j < dim; j++)
            {
                counter++;
                //double sol_diff = calcDistDif(solution, i, j, dim, distanceMatrix);
                //przy zamianie elementow
                double sol_diff = calcDistDifSwap(solution, i, j, dim, distanceMatrix);

                //jesli uzyskalismy poprawe, to nalezy dokonac zamiany od razu tak samo jak w greedy
                if (sol_diff < 0)
                {
                    //odwrocenie luku - trzeba odwrocic wszyskie elementy pomiedzy i a j
                    //swapElements(i, j, dim, solution);
                    //przy zamianie par - w elementach nieobowiazkowych jest, zeby porownac te 2 sasiedztwa
                    swap(solution[i], solution[j]);

                    //aktualizujemy aktualna dlugosc sciezki, zeby potem nie trzeba jej bylo ponownie obliczac
                    result.actual += sol_diff;
                    actualizeResult(result, dim, solution, counter, 1);
                    counter = 0;
                    //zapamietujemy najlepszy wynik  witeracji, bo niekoniecznie konczymy w optimum
                    if (result.actual < actualBest) actualBest = result.actual;
                    noIterationsNoImprovement = 0;

                }
                else    //jesli nie uzyskalismy poprawy to z pewnym prawdopodobienstwem i tak dokonujemy zamiany
                {
                    //licze prawdopodobienstwo zamiany na podstawie roznicy w jakosci oraz temperatury c
                    prob = exp(-(sol_diff/min_diff)/c);
                    //losuje liczbe od 0 do 1
                    random_number = ((double) rand() / (RAND_MAX));
                    //jezeli losowanie zadziałało robie zamiane
                    if(prob > random_number)
                    {
                        //robie wszystko to samo co w zamianie greedy:
                        //odwrocenie luku - trzeba odwrocic wszyskie elementy pomiedzy i a j
                        //swapElements(i, j, dim, solution);
                        //przy zamianie par - w elementach nieobowiazkowych jest, zeby porownac te 2 sasiedztwa
                        swap(solution[i], solution[j]);

                        //aktualizujemy aktualna dlugosc sciezki, zeby potem nie trzeba jej bylo ponownie obliczac
                        result.actual += sol_diff;
                        actualizeResult(result, dim, solution, counter, 1);
                        counter = 0;
                    }
                    //zwiekszam liczbe iteracji bez porpawy (nawet jesli zamieniamy, to dostajemy rozwiazanie gorsze)
                    noIterationsNoImprovement++;
                }

                //niezaleznie od tego co sie wydarzylo trzeba zmienic stan parametrow
                //zmniejszenie temperatury co L iteracji i reset wartośc currentL
                currentL++;
                if (currentL == L)
                {
                    c = temp(c, alfa);
                    currentL = 0;
                }
                //przerwanie algorytmu po P * L iteracjach bez poprawy
                if (noIterationsNoImprovement >= L * P)
                {
                    noImprovement = true;
                    i = dim;
                    break;
                }
            }
        }
        current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
    } while (/*current_time < algorithmTime && */!noImprovement && c >= min_probability);

    //if(noImprovement) cout<<"n"<<endl;
    //if (c < min_probability) cout<<"c"<<endl;

    if (result.iterations < 10)
    {
        result.firstBestResults[result.iterations] = result.best;
        //result.firstWorstResults[result.iterations] = result.worst;
        result.first10Results[result.iterations] = actualBest;
        result.firstSolutionNo[result.iterations] = result.solutionNo;
    }
    //cout<<result.best<<endl;
}



int main()
{
    srand(time(0));
    double current_time, inner_current_time;
    int n = 10000;
    int permutation [n];
    int i, j;
    for (i = 0; i < n; i++)
    {
        permutation[i]=i;
    }


    string instances[9] = {"berlin52.tsp", "st70.tsp", "eil76.tsp", "rd100.tsp", "kroA100.tsp", "lin105.tsp", "ch150.tsp", "pcb442.tsp", "pr1002.tsp"};
    double optimal[9] = {7542, 675, 538, 7910, 21282, 14379, 6528, 50778, 259045};
    ifstream datafile;
    ofstream output, output3, output4, output5, outputH, outputSwap;
    output.open ("results2.txt");
    output3.open("results3.txt");
    output4.open("results4.txt");
    output5.open("results5.txt");
    outputH.open("resultsH.txt");
    outputSwap.open("resultsSwap.txt");
    string line;
    string name, type, comment, dimension, edge_weight_type, node_coord_section;
    int instanceNo = 9;
    output << instanceNo << endl;
    //ostatnie 2 instancje za dlugo dzialaja, zeby dla nich analizowac punkt 3 i 4
    output3 << min(instanceNo, 7) << endl;
    output4 << min(instanceNo, 6) << endl;
    output5 << min(instanceNo, 7) << endl;
    outputH << instanceNo << endl;
    outputSwap << instanceNo << endl;
    output3 << iterations_p3 << endl;
    output4 << iterations_p4 << endl;
    output5 << iterations_p5 << endl;
    //TODO
    //pozniej zmienic k, bo na razie tylko berlin jest wczytywany
    for (int k = 0; k < instanceNo; k++)
    {
        output << instances[k] << endl;
        output << optimal[k] << endl;
        output3 << instances[k] << endl;
        output4 << instances[k] << endl;
        output5 << instances[k] << endl;
        outputH << instances[k] << endl;
        outputSwap << instances[k] << endl;
        cout << instances[k] << endl;
        datafile.open("instances/" + instances[k]);
        getline(datafile, name);
        getline(datafile, type);
        getline(datafile, comment);
        getline(datafile, dimension);
        getline(datafile, edge_weight_type);
        getline(datafile, node_coord_section);

        int dim = extractIntegerWords(dimension);
        cout<<dim<<endl;

        point tab[dim];
        int point_it = 0;
        while (!getline(datafile, line).eof())
        {
            if(line.compare("EOF")!=0)
            {
                //cout<<line<<endl;
                stringstream ss;
                ss << line;
                ss>>tab[point_it].index;
                ss>>tab[point_it].x;
                ss>>tab[point_it++].y;
            }
        }
        datafile.close();

        int solution[dim];

        double **distanceMatrix = new double*[dim];
        for(i = 0; i < dim; i++) {
            distanceMatrix[i] = new double[dim];
            solution[i] = i;
        }

        double min_diff;
        min_diff = createDistanceMatrix(dim, distanceMatrix, tab);
        int solutionGS[2][100][dim];

        //przechowuje info o najlepszych, najgorszych rozw. itp.
        algorithmResult results[7];

        //kazdy algorytm ma swoja petle, bo kazdy ma dzialac tyle samo czasu
        //RW, G, S - algorytm jest puszczany z losowa poczatkowa permutacja, wykonuje sie okreslona liczbe razy (RW) albo poki nie dojdzie do optimum (G, S)
        //potem jest puszczany ponownie z inna permutacja az nie skonczy sie czas
        //H losuje poczatkowe miasto - zrobione, zeby za kazdym razem bylo ono inne (dzieki startCity)
        //jesli skonstruuje rozwiazanie dla wszystkich mozliwych miast poczatkowych, to algorytm jest powtarzany, aby zmierzyc czas
        int algorithmTime = 10, iterationNo = 10;


        time_t time_start;
        //greedy - 3
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[3], dim, solution, distanceMatrix, "G");
        do
        {
            greedySwap(dim, solution, distanceMatrix, results[3], time_start, algorithmTime);
            randomSolution(dim, solution);
            results[3].actual = calcSolutionDistance(dim, solution, distanceMatrix);
            if (results[3].iterations < iterations_p5)
            {
                for(i=0; i<dim; i++) solutionGS[0][results[3].iterations][i] = solution[i];
            }

            results[3].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[3].iterations < iterationNo);
        results[3].time = current_time / double(results[3].iterations);
        results[3].iterationTime = current_time / double(results[3].iterations);

        //steepest - 4
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[4], dim, solution, distanceMatrix, "S");
        do
        {
            steepestSwap(tab, dim, solution, distanceMatrix, results[4], time_start, algorithmTime);
            randomSolution(dim, solution);
            results[4].actual = calcSolutionDistance(dim, solution, distanceMatrix);
            if (results[4].iterations < iterations_p5)
            {
                for(i=0; i<dim; i++) solutionGS[1][results[4].iterations][i] = solution[i];
            }

            results[4].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[4].iterations < iterationNo);
        results[4].time = current_time / double(results[4].iterations);
        results[4].iterationTime = current_time / double(results[4].iterations);

        //tabusearch - 5
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[5], dim, solution, distanceMatrix, "TS");
        do
        {
            tabuSearch(tab, dim, solution, distanceMatrix, results[5], time_start, algorithmTime);
            randomSolution(dim, solution);
            results[5].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[5].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[5].iterations < iterationNo);
        results[5].time = current_time / double(results[5].iterations);
        results[5].iterationTime = current_time / double(results[5].iterations);

        double iteration_time = results[3].time * 0.8;
        time_t innerLoopTime;

        //Simulated Annealing - 6
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[6], dim, solution, distanceMatrix, "SA");
        do
        {
            SimulatedAnnealing(tab, dim, solution, distanceMatrix, results[6], time_start, algorithmTime, min_diff);
            randomSolution(dim, solution);
            results[6].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[6].iterations ++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[6].iterations < iterationNo);
        results[6].time = current_time / double(results[6].iterations);
        results[6].iterationTime = current_time / double(results[6].iterations);


        //H heuristics - 0
        time_start = clock();
        results[0].name = "H";
        vector <int> startCity;
        do
        {
            innerLoopTime = clock();
            startCity.clear();
            //do
            //{
                for (i = 0; i < dim; i++) startCity.push_back(i);
                do
                {
                    HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
                    current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
                    inner_current_time = (clock() - innerLoopTime) / double(CLOCKS_PER_SEC);
                } while (current_time < algorithmTime && (startCity.size() > 0) && inner_current_time < iteration_time);
            //} while (inner_current_time < iteration_time);

            if (results[0].iterations < 10)
            {
                results[0].firstBestResults[results[0].iterations] = results[0].best;
                //results[0].firstWorstResults[results[0].iterations] = results[0].worst;
                //jako rozwiazanie z iteracji zapisywana jest srednia rozwiazan z tej iteracji
                results[0].first10Results[results[0].iterations] = results[0].sum / results[0].steps;
                results[0].firstSolutionNo[results[0].iterations] = results[0].steps;
                results[0].sum = 0;
                results[0].steps = 0;
            }
            results[0].iterations ++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[0].iterations < iterationNo);
        results[0].time = current_time / double(results[0].solutionNo);
        results[0].iterationTime = current_time / double(results[0].iterations);

        //random - 1
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[1], dim, solution, distanceMatrix, "R");
        double distance = calcSolutionDistance(dim, solution, distanceMatrix);
        results[1].sum = distance;
        results[1].actual = distance;
        results[1].solutionNo++;
        results[1].steps++;
        //results[1].worst = distance;
        do
        {
            innerLoopTime = clock();
            do
            {
                randomSolution(dim, solution);
                results[1].actual = calcSolutionDistance(dim, solution, distanceMatrix);
                actualizeResult(results[1], dim, solution, 1, 1);
                current_time = (clock() - innerLoopTime) / double(CLOCKS_PER_SEC);
            } while (current_time < iteration_time);

            if (results[1].iterations < 10)
            {
                results[1].firstBestResults[results[1].iterations] = results[1].best;
                //results[1].firstWorstResults[results[1].iterations] = results[1].worst;
                //jako rozwiazanie z iteracji zapisywana jest srednia rozwiazan z tej iteracji
                results[1].first10Results[results[1].iterations] = results[1].sum / results[1].steps;
                results[1].firstSolutionNo[results[1].iterations] = results[1].steps;
                results[1].sum = 0;
                results[1].steps = 0;
            }
            results[1].iterations ++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[1].iterations < 10);
        results[1].time = current_time / double(results[1].solutionNo);
        results[1].iterationTime = current_time / double(results[1].iterations);

        //random walk - 2
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[2], dim, solution, distanceMatrix, "RW");
        do
        {
            randomWalkSwap(tab, dim, solution, distanceMatrix, results[2], iteration_time);
            randomSolution(dim, solution);
            results[2].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[2].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[2].iterations < iterationNo);
        results[2].time = current_time / double(results[2].solutionNo);
        results[2].iterationTime = current_time / double(results[2].iterations);

        //wyswietlenie statystyk
        for (i = 0; i < 7; i++)
        {
            cout << results[i].name << ": ";
            cout << "best: " << results[i].best << ", ";
            cout << "best opt: " << optimal[k] / results[i].best<< ", ";
            //cout << "worst: " << results[i].worst << ", ";
            //W G i S liczony jest tylko koncowy wynik dla kazdej iteracji, w pozostalych algorytmach - kazady posredni wynik
            if (i >= 3) cout << "mean: " << results[i].sum / double(results[i].iterations) << ", ";
            else cout << "mean: " << results[i].sum / double(results[i].solutionNo) << ", ";
            cout << "iterations: " << results[i].iterations << ", ";
            //cout << "steps: " << results[i].steps << ", ";
            cout << "number reviewed solutions: " << results[i].solutionNo << ", ";
            cout << "mean time: " << results[i].time << ", ";
            cout<<endl;
        }

        //zapis do pliku punkt 2:
        for (i = 0; i < 7; i++)
        {
            //zapis nazwy algorytmu
            output << results[i].name << endl;
            //zapis rozwiazania w kolejnych iteracjach
            for (j = 0; j < 10; j++)
            {
                output << results[i].first10Results[j] << " ";
            }
            output << endl;
            //zapis najlepszych rozwiazan
            for (j = 0; j < 10; j++)
            {
                output << results[i].firstBestResults[j] << " ";
            }
            output << endl;
            //zapis liczby przejrzanych rozwiązań
            for (j = 0; j < 10; j++)
            {
                output << results[i].firstSolutionNo[j] << " ";
            }
            output << endl;
            //zapis liczby krokow - wazne tylko dla G i S
            for (j = 0; j < 10; j++)
            {
                output << results[i].firstSteps[j] << " ";
            }
            output << endl;
            //zapis sredniego czasu 1 puszczenia algorytmu
            output << results[i].time << endl;
            //zapis sredniego czasu 1 iteracji algorytmu (dla G i S to to samo co wyzej, dla pozostalych w 1 iteracji wykonywane bylo wiele powtorzen algorytmu
            output << results[i].iterationTime << endl;
        }

        if (k < 7)
        {
            //zapis do pliku punkt 3:
            for (i = 3; i < 5; i++)
            {
                //zapis nazwy algorytmu
                output3 << results[i].name << endl;
                //zapis wartosci: poczatkowe koncowe rozwiazanie
                for (j = 0; j < iterations_p3; j++)
                {
                    output3 << results[i].startSolution[j] << " " << results[i].finishSolution[j] << endl;
                }
            }

            //zapis do pliku punkt 4:
            for (i = 3; i < 5; i++)
            {
                //zapis nazwy algorytmu
                output4 << results[i].name << endl;
                //zapis wartosci: poczatkowe koncowe rozwiazanie
                for (j = 0; j < iterations_p4; j++)
                {
                    output4 << results[i].bestSolution[j] << " " << results[i].meanSolution[j] << endl;
                }
            }

            //zapis do pliku punkt 5:
            for (i = 3; i < 5; i++)
            {
                //zapis nazwy algorytmu
                output5 << results[i].name << endl;
                //zapis wartosci: poczatkowe koncowe rozwiazanie
                //zapis permutacji (do kazdego elementu +1)
                for (j = 0; j < iterations_p5; j++)
                {
                    output5 << results[i].finishSolution[j] << endl;
                    for (int l=0; l<dim; l++)
                    {
                        output5 << solutionGS[i-3][j][l] + 1 << " ";
                    }
                    output5 << endl;
                }
            }
        }

/*
        //start from H solution
        vector <int> startCity;
        for (i = 0; i < dim; i++) startCity.push_back(i);

        //greedy - 3
        time_t time_start = clock();
        HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
        initResult(results[3], dim, solution, distanceMatrix, "G");
        do
        {
            greedySwap(dim, solution, distanceMatrix, results[3], time_start, algorithmTime);
            startCity.clear();
            for (i = 0; i < dim; i++) startCity.push_back(i);
            HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
            results[3].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[3].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[3].iterations < 10);
        results[3].time = current_time / double(results[3].iterations);
        results[3].iterationTime = current_time / double(results[3].iterations);

        //steepest - 4
        time_start = clock();
        HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
        initResult(results[4], dim, solution, distanceMatrix, "S");
        do
        {
            steepestSwap(tab, dim, solution, distanceMatrix, results[4], time_start, algorithmTime);
            startCity.clear();
            for (i = 0; i < dim; i++) startCity.push_back(i);
            HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
            results[4].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[4].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[4].iterations < 10);
        results[4].time = current_time / double(results[4].iterations);
        results[4].iterationTime = current_time / double(results[4].iterations);

        //tabusearch - 5
        time_start = clock();
        HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
        initResult(results[5], dim, solution, distanceMatrix, "TS");
        do
        {
            tabuSearch(tab, dim, solution, distanceMatrix, results[5], time_start, algorithmTime);
            startCity.clear();
            for (i = 0; i < dim; i++) startCity.push_back(i);
            HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
            results[5].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[5].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[5].iterations < iterationNo);
        results[5].time = current_time / double(results[5].iterations);
        results[5].iterationTime = current_time / double(results[5].iterations);

        double iteration_time = results[3].time * 0.8;
        time_t innerLoopTime;

        //Simulated Annealing - 6
        time_start = clock();
        HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
        initResult(results[6], dim, solution, distanceMatrix, "SA");
        do
        {
            SimulatedAnnealing(tab, dim, solution, distanceMatrix, results[6], time_start, algorithmTime, min_diff);
            startCity.clear();
            for (i = 0; i < dim; i++) startCity.push_back(i);
            HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
            results[6].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[6].iterations ++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[6].iterations < iterationNo);
        results[6].time = current_time / double(results[6].iterations);
        results[6].iterationTime = current_time / double(results[6].iterations);

        for (i = 3; i < 7; i++)
        {
            //zapis nazwy algorytmu
            outputH << results[i].name << endl;
            //zapis rozwiazania w kolejnych iteracjach
            for (j = 0; j < 10; j++)
            {
                outputH << results[i].first10Results[j] << " ";
            }
            outputH << endl;
            //zapis najlepszych rozwiazan
            for (j = 0; j < 10; j++)
            {
                outputH << results[i].firstBestResults[j] << " ";
            }
            outputH << endl;
            //zapis sredniego czasu 1 puszczenia algorytmu
            outputH << results[i].time << endl;
        }
*/
/*
        //kod dla swap vs arc inversion
        //greedy - 3
        time_t time_start = clock();
        randomSolution(dim, solution);
        initResult(results[3], dim, solution, distanceMatrix, "G");
        do
        {
            greedySwap(dim, solution, distanceMatrix, results[3], time_start, algorithmTime);
            randomSolution(dim, solution);
            results[3].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[3].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[3].iterations < 10);
        results[3].time = current_time / double(results[3].iterations);
        results[3].iterationTime = current_time / double(results[3].iterations);

        //steepest - 4
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[4], dim, solution, distanceMatrix, "S");
        do
        {
            steepestSwap(tab, dim, solution, distanceMatrix, results[4], time_start, algorithmTime);
            randomSolution(dim, solution);
            results[4].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[4].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime || results[4].iterations < 10);
        results[4].time = current_time / double(results[4].iterations);
        results[4].iterationTime = current_time / double(results[4].iterations);

        //zamiana elementow zamiast odwrocenia luku
        for (i = 3; i < 5; i++)
        {
            //zapis nazwy algorytmu
            outputSwap << results[i].name << endl;
            //zapis rozwiazania w kolejnych iteracjach
            for (j = 0; j < 10; j++)
            {
                outputSwap << results[i].first10Results[j] << " ";
            }
            outputSwap << endl;
            //zapis najlepszych rozwiazan
            for (j = 0; j < 10; j++)
            {
                outputSwap << results[i].firstBestResults[j] << " ";
            }
            outputSwap << endl;
            //zapis sredniego czasu 1 puszczenia algorytmu
            outputSwap << results[i].time << endl;
        }*/


    }
    output.close();
    output3.close();
    output4.close();
    output5.close();
    outputH.close();
    outputSwap.close();
    return 0;
}
