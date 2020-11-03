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

int main()
{
    srand(time(0));
    int counter = 0, current_time, n = 10000, i;
    int permutation [n];
    for (i = 0; i < n; i++)
    {
        permutation[i]=i;
    }

    string instances[9] = {"berlin52.tsp", "ch150.tsp", "eil76.tsp", "kroA100.tsp", "lin105.tsp", "pcb442.tsp", "pr1002.tsp", "rd100.tsp", "st70.tsp"};
    for (int k = 0; k < 9; k++)
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

        time_t time_start = clock();
        do
        {
            generatePermutation(dim, tab);
            counter++;
            current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
        } while (current_time < 1);

        cout << "mean algorithm time: " << current_time / double(counter) << endl;
    }
    return 0;
}
