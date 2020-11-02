#include<iostream>
#include <vector>
#include<cstdlib>
#include<ctime>
#include <fstream>
#include <sstream>
#include <math.h>
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

float distance(point p1, point p2)
{
    //dystans euklidesowy (sa tez inne)
    return sqrt((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
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

    ifstream datafile;
    datafile.open("instances/berlin52.tsp");
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

    time_t time_start = clock();
    do
    {
        generatePermutation(dim, tab);
        counter++;
        current_time = (clock() - time_start) / double(CLOCKS_PER_SEC);
    }while (current_time < 1);

    cout << "mean algorithm time: " << current_time / double(counter) << endl;
    //for (int i = 0; i < n; i++) cout << permutation[i] << " ";

    return 0;
}
