#ifndef DIJKSTRA_H_INCLUDED
#define DIJKSTRA_H_INCLUDED

#define STRLEN 85
#define INFINITY 20000
#define ARRAYLEN 20

typedef struct City{
    char city[STRLEN];
    int index;
} City;

typedef struct Graph{
    int graphSize;
    struct List* array;
} Graph;

typedef struct HeapNode{
    int  cityInt;
    int distance;
} HeapNode;

typedef struct Heap{
    int capacity, size;
    int *position;
    struct HeapNode **array; // ** creates a pointer to a pointer
} Heap;

// Structure by Prof. John Robinson
typedef struct Intertown_distance{
    char start[STRLEN];
    char end[STRLEN];
    int distance;
    int arraylength;
} intertown_distance;

typedef struct ListNode{
    int destination, distance;
    struct ListNode* next;
} ListNode;

typedef struct List{
    struct ListNode *head;
} List;

// adds an edge to an undirected graph
void addEdge(Graph* graph, int origin, int destination, int distance);

// take a structure of edges and swap city names for their index numbers
void cityNamesToNumbers(intertown_distance edges[], City cities[]);

// take a city name an output its associated index
int cityNumberFromName(char cityName[], intertown_distance edges[], City cities[]);

// given an intertown_distance array, counts the number of unique cities in it
int countCities(intertown_distance structure[]);

// returns the number of lines in a file
// used to calculate structure array lengths
int countLines(char *filename);

// calculates distances of shortest path from origin to destination vertices. function is O(HLogN)
int dijkstra(Graph* graph, int origin, int destination, int pathArr[]);

// generate a valid heap at given index
void generateHeap(Heap* Heap, int index);

// creates a graph of graphSize vertices
Graph* graphConstructor(int graphSize);

// create a Heap
Heap* heapConstructor(int capacity);

// create a new Heap Node
HeapNode* HeapNodeConstructor(int cityInt, int distance);

// check if the given Heap is empty or not
int isEmpty(Heap* Heap);

//  check if given vertex index cityInt is in heap or not
int isInHeap(Heap *Heap, int cityInt);

// create a new list node and initialize it
ListNode* ListNodeConstructor(int destination, int distance);

// pop the smallest node from heap
HeapNode* popSmallest(Heap* Heap);

// take a list of edges and use it to populate an array of unique cities with index numbers
void populateCityNames(intertown_distance structure[], City cities[], char *filename);

// prints the shortest path from start to end
void printArr(int arr[], int n, int end, int start, City cityList[]);

// prints a city list with their index numbers
void printCities(City cities[], int size);

// prints an intertown_distance array
void printintertown_distance(intertown_distance structure[]);

// reads a compatible tab-delimited file into an intertown_distance structure
void readFileTointertown_distance(char * filename, intertown_distance structure[]);

// Swap two nodes of heap, used in generateHeap
void swapHeapNode(HeapNode** a, HeapNode** b);

// decrease distance value of given vertex index cityInt
void updateDist(Heap* Heap, int cityInt, int distance);

// writes an intertown_distance structure to a file
void writeFileintertown_distance(char * filename, intertown_distance structure[]);

#endif // DIJKSTRA_H_INCLUDED
