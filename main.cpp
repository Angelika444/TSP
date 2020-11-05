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

void generatePermutation (int n, point permutation[])
{
    for (int i = n; i > 0; i--)
    {
        int index = rand() % i;
        swap(permutation[index], permutation[i - 1]);
    }
}

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

float calcDistance(point p1, point p2)
{
    //dystans euklidesowy (sa tez inne)
    double x1 = p1.x - p2.x, x2 = p1.y - p2.y;
    return sqrt(x1 * x1 + x2 * x2);
}

void generateFirstSolution(point tab[], int dim, int solution[])
{
    //solution is a permutation with cities index in tab, not field index in struct, so from 0 to dim-1
    int i, j, solution_index = rand() % dim, index_min_dinstance, min_distance, current_distance, sum_distance = 0;
    solution[0] = solution_index;
    vector <int> not_visited;
    for (i = 0; i < dim; i++)
    {
        not_visited.push_back(i);
    }
    not_visited.erase(not_visited.begin() + solution_index);
    for (i = 1; i < dim ; i++)
    {
        index_min_dinstance = 0;
        min_distance = calcDistance(tab[solution_index], tab[not_visited[index_min_dinstance]]);
        for (j = 1; j < not_visited.size(); j++)
        {
            current_distance = calcDistance(tab[solution_index], tab[not_visited[j]]);
            if (current_distance < min_distance) {
                min_distance = current_distance;
                index_min_dinstance = j;
            }
        }
        solution[i] = not_visited[index_min_dinstance];
        not_visited.erase(not_visited.begin() + index_min_dinstance);
        sum_distance += min_distance;
    }
    sum_distance += calcDistance(tab[solution[0]], tab[solution[dim - 1]]);
    cout << "distance: " << sum_distance <<endl;
}


void greedySwap(point tab[], int dim, int solution[])
{
    //greedy
    bool swaped = false;
    //generuj po kolei elementy sasiedztwa
    for(int i = 0; i < dim; i++)
    {
        for(int j = i+1; j < dim-1; j++)
        {
            //cout do usuniecia
            cout<<"i, j: "<<i<<" "<<j<<" "<<endl;
            //obliczam roznice w jakosci rozwiazania zamiast liczyc wszystko od nowa
            double sol_diff = 0.0;
            //odwrocenie kolejnosci luku - odejmuemy 2 stare drogi i dodajemy 2 nowe drogi
            //sprawdzamy przy tym czy nie trafilismy na skrajne wartosci i=0, j=dim-1

            if (i != 0){
                    sol_diff -= calcDistance(tab[solution[i-1]], tab[solution[i]]);
                    sol_diff += calcDistance(tab[solution[i-1]], tab[solution[j]]);
            }
            else {
                    sol_diff -= calcDistance(tab[solution[dim-1]], tab[solution[i]]);
                    sol_diff += calcDistance(tab[solution[dim-1]], tab[solution[j]]);
            }

            if (j != dim-1){
                    sol_diff -= calcDistance(tab[solution[j]], tab[solution[j+1]]);
                    sol_diff += calcDistance(tab[solution[i]], tab[solution[j+1]]);
            }
            else {
                    sol_diff -= calcDistance(tab[solution[j]], tab[solution[0]]);
                    sol_diff += calcDistance(tab[solution[i]], tab[solution[0]]);
            }


            //jesli uzyskalismy poprawe, to nalezy dokonac zamiany od razu
            //cout do usuniecia
            cout<<"sol_diff: "<<sol_diff<<endl;
            if (sol_diff < 0)
            {
                swap(solution[i], solution[j]);
                swaped = true;
                break;
            }
        }
        //cout do usuniecia
        for(int h = 0; h<dim; h++) cout<<solution[h]<<" ";
        cout<<endl;
        if (swaped) break;
    }
    //cout do usuniecia
    cout<<endl;
}


void steepestSwap(point tab[], int dim, int solution[])
{
    //steepest
    double best_diff = 0.0;
    double sol_diff = 0.0;
    int best_i=0, best_j=0;
    //generuj po kolei elementy sasiedztwa
    for(int i = 0; i < dim; i++)
    {
        for(int j = i+1; j < dim-1; j++)
        {
            //cout do usuniecia
            cout<<"i, j: "<<i<<" "<<j<<" "<<endl;
            //obliczam roznice w jakosci rozwiazania zamiast liczyc wszystko od nowa
            sol_diff = 0.0;
            //odwrocenie kolejnosci luku - odejmuemy 2 stare drogi i dodajemy 2 nowe drogi
            //sprawdzamy przy tym czy nie trafilismy na skrajne wartosci i=0, j=dim-1

            if (i != 0){
                    sol_diff -= calcDistance(tab[solution[i-1]], tab[solution[i]]);
                    sol_diff += calcDistance(tab[solution[i-1]], tab[solution[j]]);
            }
            else {
                    sol_diff -= calcDistance(tab[solution[dim-1]], tab[solution[i]]);
                    sol_diff += calcDistance(tab[solution[dim-1]], tab[solution[j]]);
            }

            if (j != dim-1){
                    sol_diff -= calcDistance(tab[solution[j]], tab[solution[j+1]]);
                    sol_diff += calcDistance(tab[solution[i]], tab[solution[j+1]]);
            }
            else {
                    sol_diff -= calcDistance(tab[solution[j]], tab[solution[0]]);
                    sol_diff += calcDistance(tab[solution[i]], tab[solution[0]]);
            }
            //cout<<sol_diff<<endl;


            //jesli uzyskalismy poprawe, to nalezy zapisac jej wartosc do porownania z innymi
            if(sol_diff < best_diff)
            {
                best_diff = sol_diff;
                best_i = i;
                best_j = j;
            }
            //cout do usuniecia
            cout<<"sol_diff: "<<sol_diff<<endl;
        }
    }
    //dopiero teraz po przejrzeniu ca³ego sasiedztwa robimy swapa
    if(best_diff < 0)   //zabezpieczam sie przed wywaleniem algorytmu, jezeli bylibysmy ju¿ w optimum, i nie ma lepszego rozwiazania
    {
        swap(solution[best_i], solution[best_j]);
        //cout do usuniecia
        cout<<"best_i, best_j, best_diff: "<<" "<<best_i<<", "<<best_j<<" "<<best_diff<<endl;
        for(int h = 0; h<dim; h++) cout<<solution[h]<<" ";
        cout<<endl;
    }
    else
    {
        //st¹d mo¿na w sumie przerywac dzi³anie zewnêtrznej pêtli while ale nie wiem jak to w eksperymentach potem dobrze opisaæ
        //w greedym w sumie tez by mozna gdzies taka logike zaimplementowac
        ;
    }
}

void randomWalkSwap(point tab[], int dim, int solution[])
{
    //random walk
    double sol_diff = 0.0;
    int i=0, j=0;
    bool swaped = false;
    int counter = 0;
    //dopoki nie zamienisz (lub nie minie okreslona liczba losowan zeby algorytm nie zapetlil sie w optimum)
    while(!swaped | counter > 1000)
    {
        counter++;
        //losuj pary i, j

        /*
        //taki sposob moze byc szybsz, ale nie wiem czy idealnie losowy
        i = rand() % (dim-1);
        cout<<"division_by i "<<i<<endl;
        j = rand() % (dim-i-1) + (i+1);
        */

        //alternatywny sposob losowania
        //nie chcemy takich samych wartosci bo taka zamiana nic
        do{
            i = rand() % dim;
            j = rand() % dim;
        }while(i==j);
        //robie swapa jezeli i > j zeby nie sprawdzac dodatkowych warunkow brzegowych w liczeniu funkcji kosztu
        if (i>j){
            i=i+j;
            j=i-j;
            i=i-j;
        }
        //cout do usuniecia
        cout<<"i, j: "<<i<<" "<<j<<" "<<endl;
        //obliczam roznice w jakosci rozwiazania zamiast liczyc wszystko od nowa
        sol_diff = 0.0;
        //odwrocenie kolejnosci luku - odejmuemy 2 stare drogi i dodajemy 2 nowe drogi
        //sprawdzamy przy tym czy nie trafilismy na skrajne wartosci i=0, j=dim-1

        if (i != 0){
                sol_diff -= calcDistance(tab[solution[i-1]], tab[solution[i]]);
                sol_diff += calcDistance(tab[solution[i-1]], tab[solution[j]]);
        }
        else {
                sol_diff -= calcDistance(tab[solution[dim-1]], tab[solution[i]]);
                sol_diff += calcDistance(tab[solution[dim-1]], tab[solution[j]]);
        }

        if (j != dim-1){
                sol_diff -= calcDistance(tab[solution[j]], tab[solution[j+1]]);
                sol_diff += calcDistance(tab[solution[i]], tab[solution[j+1]]);
        }
        else {
                sol_diff -= calcDistance(tab[solution[j]], tab[solution[0]]);
                sol_diff += calcDistance(tab[solution[i]], tab[solution[0]]);
        }
        //cout<<sol_diff<<endl;


        //jesli uzyskalismy poprawe, to nalezy natychmiast zamienic
        cout<<"sol_diff: "<<sol_diff<<endl;
        if (sol_diff < 0)
        {
            cout<<"if\n";
            swap(solution[i], solution[j]);
            swaped = true;
        }
    }
    //cout do usuniecia
    for(int h = 0; h<dim; h++) cout<<solution[h]<<" ";
    cout<<endl;
}



int main()
{
    srand(time(0));
    int counter = 0, current_time, n = 10000;
    int permutation [n];
    for (int i = 0; i < n; i++)
    {
        permutation[i]=i;
    }

    string instances[9] = {"berlin52.tsp", "ch150.tsp", "eil76.tsp", "kroA100.tsp", "lin105.tsp", "pcb442.tsp", "pr1002.tsp", "rd100.tsp", "st70.tsp"};
    for (int k = 0; k < 1; k++)
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
        generateFirstSolution(tab, dim, solution);
        for(int i =0; i<dim; i++) cout<<solution[i]<<" ";

        time_t time_start = clock();
        //petla zewnetrzna zawierajaca caly algorytm
        do
        {
            /**
            //random
            generatePermutation(dim, tab);
            */


            //random walk
            randomWalkSwap(tab, dim, solution);


            /**
            //greedy
            greedySwap(tab, dim, solution);
            */

            /**
            //steepest
            steepestSwap(tab, dim, solution);
            */




            counter++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < 100);

        cout << "mean algorithm time: " << current_time / double(counter) << endl;
    }
    return 0;
}
