#include<iostream>
#include <vector>
#include<cstdlib>
#include<ctime>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
using namespace std;


struct point{
    int index;
    float x;
    float y;
};

struct algorithmResult {
    double best = -1;
    double worst;
    double sum = 0;
    //liczba krokow - dla greedy i steepest
    int steps = 0;
    //liczba przejrzanych rozwiazan
    int numSolutions = 0;
    //aktualny sumaryczny dystans
    double actual;
    int iterations = 0;
    double time;
    string name;
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

int** createDistanceMatrix(int dim, double** matrix, point tab[])
{
    for(int i = 0; i<dim-1; i++)
    {
        matrix[i][i] = 0;
        for(int j = i+1; j<dim; j++)
        {
            matrix[i][j] = calcDistance(tab[i], tab[j]);
            matrix[j][i] = matrix[i][j];
        }
    }
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
    result.worst = distance;
    result.actual = distance;
    result.name = name;
}

void actualizeResult(algorithmResult &result, int dim, int solution[], int counter, int steps)
{
    result.sum += result.actual;
    if (result.actual < result.best) result.best = result.actual;
    else if (result.actual > result.worst) result.worst = result.actual;
    result.numSolutions += counter;
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
    int i, j, solution_index, index_min_dinstance, min_distance, current_distance, distance = 0, startCityIndex;
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

    //jezeli jest to pierwsze nasze rozwiazanie
    if (result.best == -1)
    {
        result.best = distance;
        result.worst = distance;
    }
    else if (distance < result.best) result.best = distance;
    else if (distance > result.worst) result.worst = distance;
    result.sum += distance;
    result.numSolutions++;

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
        sol_diff -= getDistance(solution[j], solution[j-1], distanceMatrix);
        sol_diff += getDistance(solution[i], solution[j-1], distanceMatrix);
        sol_diff -= getDistance(solution[i], solution[i+1], distanceMatrix);
        sol_diff += getDistance(solution[j], solution[i+1], distanceMatrix);
    }

    return sol_diff;
}

//zamiana wszystkich elementow w luku - ostatni z pierwszym az nie dojdziemy do srodka
void swapElements(int i, int j, int dim, int solution[])
{
    //jesli i=0, j=dim-1, to zamieniamy je miejscami, bo odwrocenie wszystkich elementow pomiedzy nimi daloby dokladnie to samo sasiedztwo
    if (i == 0 && j == dim - 1) swap(solution[i], solution[j]);
    else
    {
       for (int k = 0; k < ceil(double(j + i) / 2) - i; k++)
        {
            swap(solution[i + k], solution[j - k]);
        }
    }
}

void greedySwap(int dim, int solution[], double** distanceMatrix, algorithmResult &result, time_t time_start, int algorithmTime)
{
    //greedy
    bool swaped;
    time_t current_time;
    int counter = 1, steps = 0;
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
                double sol_diff = calcDistDif(solution, i, j, dim, distanceMatrix);

                //jesli uzyskalismy poprawe, to nalezy dokonac zamiany od razu
                if (sol_diff < 0)
                {
                    //odwrocenie luku - trzeba odwrocic wszyskie elementy pomiedzy i a j
                    swapElements(i, j, dim, solution);
                    //przy zamianie par - w elementach nieobowiazkowych jest, zeby porownac te 2 sasiedztwa
                    //swap(solution[i], solution[j]);

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
    } while (current_time < algorithmTime && swaped);

    //wynik koncowy jest najlepszym dotychczasowym, dlatego aktualizujemy go na koncu
    actualizeResult(result, dim, solution, counter, steps);

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

    do
    {
        counter += numIterations;
        best_diff = 0.0;
        //generuj po kolei elementy sasiedztwa
        // zewnetrzna petla musi byc do dim-1, a wewnetrzna do dim, inaczej nie bedziemy miec nigdy pary 0, dim-1, jak w greedy
        for(int i = 0; i < dim - 1; i++)
        {
            for(int j = i+1; j < dim; j++)
            {
                double sol_diff = calcDistDif(solution, i, j, dim, distanceMatrix);

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
            swapElements(best_i, best_j, dim, solution);
            //przy zamianie par - w elementach nieobowiazkowych jest, zeby porownac te 2 sasiedztwa
            //swap(solution[i], solution[j]);

            steps++;
            result.actual += best_diff;
        }
        else
        {
            swaped = false;
        }
        current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
    }
    while (current_time < algorithmTime && swaped);

    //wynik koncowy jest najlepszym dotychczasowym, dlatego aktualizujemy go na koncu
    actualizeResult(result, dim, solution, counter, steps);

}

void randomWalkSwap(point tab[], int dim, int solution[], double** distanceMatrix, algorithmResult &result, time_t time_start, int algorithmTime)
{
    //random walk
    double sol_diff;
    int i=0, j=0;
    bool swaped = false;
    int counter = 0;
    //idziemy dana sciezka tyle razy, potem losujemy nowe rozwiazanie - nowa sciezka
    while(counter < 1000)
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
        actualizeResult(result, dim, solution, 1, 0);
        swaped = true;
    }
}


int main()
{
    srand(time(0));
    double current_time;
    int n = 10000;
    int permutation [n];
    int i, j;
    for (i = 0; i < n; i++)
    {
        permutation[i]=i;
    }


    string instances[9] = {"berlin52.tsp", "ch150.tsp", "eil76.tsp", "kroA100.tsp", "lin105.tsp", "pcb442.tsp", "pr1002.tsp", "rd100.tsp", "st70.tsp"};
    double optimal[9] = {7542, 6528, 538, 21282, 14379, 50778, 259045, 7910, 675};
    //TODO
    //pozniej zmienic k, bo na razie tylko berlin jest wczytywany
    for (int k = 1; k < 9; k++)
    {
        ifstream datafile;
        datafile.open("instances/" + instances[k]);
        string line;
        string name, type, comment, dimension, edge_weight_type, node_coord_section;
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
        }

        createDistanceMatrix(dim, distanceMatrix, tab);

        //przechowuje info o najlepszych, najgorszych rozw. itp.
        algorithmResult results[5];

        //kazdy algorytm ma swoja petle, bo kazdy ma dzialac tyle samo czasu
        //RW, G, S - algorytm jest puszczany z losowa poczatkowa permutacja, wykonuje sie okreslona liczbe razy (RW) albo poki nie dojdzie do optimum (G, S)
        //potem jest puszczany ponownie z inna permutacja az nie skonczy sie czas
        //H losuje poczatkowe miasto - zrobione, zeby za kazdym razem bylo ono inne (dzieki startCity)
        //jesli skonstruuje rozwiazanie dla wszystkich mozliwych miast poczatkowych, to algorytm jest powtarzany, aby zmierzyc czas

        //H heuristics - 0
        time_t time_start = clock();
        int algorithmTime = 1;
        results[0].name = "H";
        vector <int> startCity;
        do
        {
            for (i = 0; i < dim; i++) startCity.push_back(i);
            do
            {
                HSolution(tab, dim, solution, distanceMatrix, startCity, results[0]);
                current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
            } while (current_time < algorithmTime && (startCity.size() > 0));
        } while (current_time < algorithmTime);
        results[0].iterations = results[0].numSolutions;
        results[0].time = current_time / double(results[0].iterations);

        //random - 1
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[1], dim, solution, distanceMatrix, "R");
        double distance = calcSolutionDistance(dim, solution, distanceMatrix);
        results[1].sum = distance;
        results[1].actual = distance;
        results[1].numSolutions++;
        do
        {
            randomSolution(dim, solution);
            results[1].actual = calcSolutionDistance(dim, solution, distanceMatrix);
            actualizeResult(results[1], dim, solution, 1, 0);

            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime);
        results[1].iterations = results[1].numSolutions;
        results[1].time = current_time / double(results[1].iterations);

        //random walk - 2
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[2], dim, solution, distanceMatrix, "RW");
        do
        {
            randomWalkSwap(tab, dim, solution, distanceMatrix, results[2], time_start, algorithmTime);
            randomSolution(dim, solution);
            results[2].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[2].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime);
        results[2].time = current_time / double(results[2].iterations);

        //greedy - 3
        time_start = clock();
        randomSolution(dim, solution);
        initResult(results[3], dim, solution, distanceMatrix, "G");
        do
        {
            greedySwap(dim, solution, distanceMatrix, results[3], time_start, algorithmTime);
            randomSolution(dim, solution);
            results[3].actual = calcSolutionDistance(dim, solution, distanceMatrix);

            results[3].iterations++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < algorithmTime);
        results[3].time = current_time / double(results[3].iterations);

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
        } while (current_time < algorithmTime);
        results[4].time = current_time / double(results[4].iterations);

        //wyswietlenie statystyk
        for (i = 0; i < 5; i++)
        {
            cout << results[i].name << ": ";
            cout << "best: " << results[i].best << ", ";
            cout << "worst: " << results[i].worst << ", ";
            //W RW kazdy posredni wynik moze byc najlepszy, wiec srednia liczona jest ze wszystkich przejrzanych rozwiaqzan
            if (i == 2) cout << "mean: " << results[i].sum / double(results[i].numSolutions) << ", ";
            else cout << "mean: " << results[i].sum / double(results[i].iterations) << ", ";
            cout << "iterations: " << results[i].iterations << ", ";
            cout << "steps: " << results[i].steps << ", ";
            cout << "number reviewed solutions: " << results[i].numSolutions << ", ";
            cout << "mean time: " << results[i].time << ", ";
            cout<<endl;
        }
    }
    return 0;
}
