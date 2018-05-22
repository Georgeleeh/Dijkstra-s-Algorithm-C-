#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dijkstra.h"

int dijkstraMain()
{
    // initialize filenames
    // edges
    char *distancesFile = "ukcities.txt";
    // test-pairs
    char *citypairsFile = "citypairs.txt";
    // output file
    char *outputFile = "testout.txt";

    //create array of city edge structures, with length equaling the number of lines in the file
    intertown_distance edges[countLines(distancesFile)];
    intertown_distance testpairs[countLines(citypairsFile)];

    // read files into structures
    readFileTointertown_distance(distancesFile, edges);
    //printintertown_distance(edges);
    readFileTointertown_distance(citypairsFile, testpairs);

    // create array of city structures, with length equaling the number of unique cities in edges[]
    City cities[countCities(edges)];

    // create array of city indexes in the shortest path
    int pathArr[countCities(edges)];

    // turn edges file into unordered list of unique cities and index numbers
    populateCityNames(edges, cities, distancesFile);

    // use city index numbers to convert all cities in edges() to their index numbers
    cityNamesToNumbers(edges, cities);

    // create graph with one node for each unique city
    Graph* graph = graphConstructor(countCities(edges));

    // add edges form struct to graph
    for (int i = 0; i < countLines(distancesFile); i++) {
        addEdge(graph, atoi(edges[i].start), atoi(edges[i].end), edges[i].distance);
    }

    // obtain index number for each city to test, enter it into dijkstra function
    for (int i = 0; i < countLines(citypairsFile); i++) {
        testpairs[i].distance = dijkstra(graph, cityNumberFromName(testpairs[i].start, edges, cities), cityNumberFromName(testpairs[i].end, edges, cities), pathArr);
        // print the paths calculated by Dijkstra's algorithm
        printArr(pathArr,countCities(edges),cityNumberFromName(testpairs[i].end, edges, cities),cityNumberFromName(testpairs[i].start, edges, cities),cities);
        printf("\n");
    }

    // print distances calculated by the dijkstra function
    printintertown_distance(testpairs);

    // write calculated distances to output file
    writeFileintertown_distance(outputFile,testpairs);

    return 0;
}

int timingMain()
{
    double totalTime = 0;
    double averageTime = 0;
    int samples = 1000;

    for(int i = 0; i < samples; i++)
    {
        double timeStart = (double)clock();
        dijkstraMain();
        double timeEnd = (double)clock();

        // Calculate how long we spent doing the sort
        double timeInSeconds = (timeEnd - timeStart)/CLOCKS_PER_SEC;

        totalTime += timeInSeconds;
    }

    averageTime = totalTime/samples;

    printf("\nsamples = %d, average time = %f, total time = %f", samples, averageTime, totalTime);

    return 0;
}

int main()
{

    dijkstraMain();

    return 0;

}
