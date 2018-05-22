#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "dijkstra.h"

// adds an edge to an undirected graph
void addEdge(Graph* graph, int origin, int destination, int distance){
    // add edge from origin to destination
    ListNode* newNode = ListNodeConstructor(destination, distance);
    newNode->next = graph->array[origin].head;
    // add node at head
    graph->array[origin].head = newNode;

    // graph is undirected, add an edge from destination to origin also
    newNode = ListNodeConstructor(origin, distance);
    newNode->next = graph->array[destination].head;
    graph->array[destination].head = newNode;
}

// take a structure of edges and swap city names for their index numbers
void cityNamesToNumbers(intertown_distance edges[], City cities[]){
    char temp[STRLEN];

    for (int i = 0; i < edges[0].arraylength; i++){
        int n = cityNumberFromName(edges[i].start,edges,cities);
        sprintf(temp, "%d", n);
        strcpy(edges[i].start, temp);

        n = cityNumberFromName(edges[i].end,edges,cities);
        sprintf(temp, "%d", n);
        strcpy(edges[i].end, temp);
    }
}

// take a city name an output its associated index
int cityNumberFromName(char cityName[], intertown_distance edges[], City cities[]){
    int n;

    for (int i = 0; i < countCities(edges); i++){
        if (!strcmp(cities[i].city, cityName)){
            n = cities[i].index;
            break;
        }
    }
    return n;
}

// given an intertown_distance array, counts the number of unique cities in it
int countCities(intertown_distance structure[]){
    // holds every start and end from the intertown_distance structure
    char citylist[structure[0].arraylength*2][STRLEN];
    // holds number of cities currently thought to be unique
    int citiesnumber = 2 * structure[0].arraylength;

    // copy the start and end cities from the structure to adjacent locations
    // in the citylist array
    for (int i = 0; i < structure[0].arraylength; i++){
        strcpy(citylist[2*i], structure[i].start);
        strcpy(citylist[2*i+1], structure[i].end);
    }
    // compare every city in the list to every other city
    for (int i = 0; i < 2 * structure[0].arraylength; i++){
        for (int j = i+1; j < 2 * structure[0].arraylength; j++){
            // every time a duplicate city is found
            if (strcmp(citylist[i], citylist[j]) == 0){
                // empty the duplicate string
                memset(citylist[j],0,sizeof(citylist[j]));
                // decrement number of cities currently believed to be unique
                if (strcmp(citylist[i], "") != 0){
                    citiesnumber--;
                }
            }
        }
    }
    return citiesnumber;
}

// returns the number of lines in a file
// used to calculate structure array lengths
int countLines(char *filename){
    // count the number of lines in the file called filename
    FILE *infile = fopen(filename,"r");
    int ch=0;
    int lines=0;

    // break-soft if filename is null
    if (infile == NULL){
        printf("Error ( countlines() ): filename == NULL");
        return -1;
    }
    // loop until the end of the file
    while(!feof(infile))
    {
        // get each character one-by-one
        ch = fgetc(infile);
        // count each new line to get the total number of lines
        if(ch == '\n')
        {
            lines++;
        }
    }
    fclose(infile);
    return lines;
}

// calculates distances of shortest path from origin to destination vertices
//O(|H| log |N|)
int dijkstra(Graph* graph, int origin, int destination, int pathArr[]){

    // number of vertices in graph
    int graphSize = graph->graphSize;
    // distances used to pick smallest weight edge
    int distance[graphSize];
    // Heap represents H
    Heap* Heap = heapConstructor(graphSize);

    // Initialize heap with all vertices and infinite distance values
    for (int i = 0; i < graphSize; i++)
    {
        distance[i] = INFINITY;
        Heap->position[i] = i;
        Heap->array[i] = HeapNodeConstructor(i, distance[i]);
    }

    // set distance value of origin to 0 so that it is popped first
    Heap->array[origin] = HeapNodeConstructor(origin, distance[origin]);
    Heap->position[origin] = origin;
    distance[origin] = 0;
    updateDist(Heap, origin, distance[origin]);
    // at first, size of heap equals graphSize
    Heap->size = graphSize;

    // heap contains all nodes with infinite distance
    while (!isEmpty(Heap)){
        // pop the next node from the heap and store its index
        HeapNode* HeapNode = popSmallest(Heap);
        int index = HeapNode->cityInt;
        ListNode* dijkList = graph->array[index].head;

        // traverse through all adjacent vertices of the popped vertex and update their distance values
        while (dijkList != NULL){
            int cityInt = dijkList->destination;
            // if shortest distance to cityInt is still temporary, and distance to cityInt
            // through the popped vertex is less than its previously calculated distance
            if (isInHeap(Heap, cityInt) && distance[index] != INFINITY && dijkList->distance + distance[index] < distance[cityInt]){
                distance[cityInt] = distance[index] + dijkList->distance;
                // update heap distance value too
                updateDist(Heap, cityInt, distance[cityInt]);
                // add city index of current finalized node to path array
                pathArr[cityInt] = index;
                // check current HeapNode index against destination
                if (cityInt == destination){
                    // return final distance of destination from origin
                    return distance[destination];
                }
            }
            dijkList = dijkList->next;
        }
    }
    // return -1 on failure to calculate distance to destination
    return -1;
}

// generate a valid heap at given index. Recursive function. O(NlogN)?
void generateHeap(Heap* Heap, int index){
    int smallest, left, right;
    smallest = index;
    left = 2 * index + 1;
    right = 2 * index + 2;

    if (left < Heap->size && Heap->array[left]->distance < Heap->array[smallest]->distance){
        smallest = left;
    }
    if (right < Heap->size && Heap->array[right]->distance < Heap->array[smallest]->distance){
        smallest = right;
    }

    if (smallest != index){
        // The nodes to be swapped in heap
        HeapNode *smallestNode = Heap->array[smallest];
        HeapNode *indexNode = Heap->array[index];

        // Swap positions
        Heap->position[smallestNode->cityInt] = index;
        Heap->position[indexNode->cityInt] = smallest;
        // Swap nodes
        swapHeapNode(&Heap->array[smallest], &Heap->array[index]);

        generateHeap(Heap, smallest);
    }
}

// creates a graph of graphSize vertices
Graph* graphConstructor(int graphSize){
    Graph* graph = (Graph*) malloc(sizeof(Graph));
    graph->graphSize = graphSize;
    // create an array of adjacency lists, of size graphSize
    graph->array = (List*) malloc(graphSize * sizeof(List));

    for (int i = 0; i < graphSize; i++){
        graph->array[i].head = NULL;
    }
    return graph;
}

// create a Heap
Heap* heapConstructor(int capacity){
    Heap* Heap = (struct Heap*)malloc(sizeof( Heap));
    Heap->position = (int *)malloc(capacity * sizeof(int));
    Heap->array = (HeapNode**) malloc(capacity * sizeof(HeapNode*));
    Heap->capacity = capacity;
    Heap->size = 0;
    return Heap;
}

// create a new Heap Node
HeapNode* HeapNodeConstructor(int cityInt, int distance){
    HeapNode* HeapNode = (struct HeapNode*) malloc(sizeof( HeapNode));
    HeapNode->cityInt = cityInt;
    HeapNode->distance = distance;
    return HeapNode;
}

// check if the given Heap is empty or not
int isEmpty(Heap* Heap){
    return Heap->size == 0;
}

//  check if given vertex index cityInt is in heap or not
int isInHeap(Heap *Heap, int cityInt){
   if (Heap->position[cityInt] < Heap->size){
     return 1;
   }
   return 0;
}

// create a new list node and initialize it
ListNode* ListNodeConstructor(int destination, int distance){
    ListNode* newNode = (ListNode*) malloc(sizeof(ListNode));
    newNode->destination = destination;
    newNode->distance = distance;
    newNode->next = NULL;
    return newNode;
}

// pop the smallest node from heap
// O(log N)
HeapNode* popSmallest(Heap* Heap){
    if (isEmpty(Heap)){
        return NULL;
    }

    // store the root node, which by design will be the smallest
    HeapNode* root = Heap->array[0];
    // move the last node to the root position
    HeapNode* lastNode = Heap->array[Heap->size - 1];
    Heap->array[0] = lastNode;
    // Update position of last node
    Heap->position[root->cityInt] = Heap->size-1;
    Heap->position[lastNode->cityInt] = 0;

    // Reduce heap size and make valid heap from root
    Heap->size--;
    generateHeap(Heap, 0);

    return root;
}

// take a list of edges and use it to populate an array of unique cities with index numbers
void populateCityNames(intertown_distance structure[], City cities[], char *filename){

    // holds every start and end from the intertown_distance structure
    char citylist[structure[0].arraylength*2][STRLEN];
    // holds number of number of confirmed unique cities
    int arrayposition = 0;

    // copy the start and end cities from the structure to adjacent locations
    // in the citylist array
    for (int i = 0; i < structure[0].arraylength; i++){
        strcpy(citylist[2*i], structure[i].start);
        strcpy(citylist[2*i+1], structure[i].end);
    }

    // compare every city in the list to every other city
    // every time a duplicate city is found empty the duplicate string
    for (int i = 0; i < 2 * structure[0].arraylength; i++){
        for (int j = i+1; j < 2 * structure[0].arraylength; j++){
            if (strcmp(citylist[i], citylist[j]) == 0){
                    strcpy(citylist[j], "");
            }
        }
    }

    // search list of cities with only unique cities and blank locations
    for (int i = 0; i < 2 * structure[0].arraylength; i++){
        //if location is not empty
        if (strcmp(citylist[i], "") != 0){
            // copy name of non-empty location to array of unique cities
            cities[arrayposition].index = arrayposition;
            strcpy(cities[arrayposition].city, citylist[i]);
            arrayposition++;
        }
    }
}

// prints the shortest path from start to end
void printArr(int arr[], int n, int end, int start, City cityList[]){
    int x = end;
    int arrayLen = 0;
    // array of city names in the path
    char temp[ARRAYLEN][STRLEN];
    // follow the path back from destination to origin
    while (x != start){
        // copy name of current city to the path names array
        strcpy(temp[arrayLen], cityList[x].city);
        x = arr[x];
        arrayLen++;
    }
    // print the array in reverse order, which will print
    // the cities origin to destination
    printf("%s -> ", cityList[start].city);
    for (int i = arrayLen-1; i > 0; i--){
        printf("%s -> ", temp[i]);
    }
    printf("%s\n", cityList[end].city);
}

// prints a city list with their index numbers
void printCities(City cities[], int size){
    //print unordered unique cities with index numbers
    for (int i = 0; i < size; i++) {
        printf("%d\t%s\n", i, cities[i].city);
    }
}

// prints an intertown_distance array
void printintertown_distance(intertown_distance structure[]){
    int i = 0;

    while(i < structure[0].arraylength){
        printf("%s\t", structure[i].start);
        printf("%s\t", structure[i].end);
        printf("%d\n", structure[i].distance);
        i++;
    }
};

// reads a compatible tab-delimited file into an intertown_distance structure
void readFileTointertown_distance(char * filename, intertown_distance structure[]) {
    // variable to hold entire line of file to be  read
    char line[STRLEN];
    // temp variables to hold important subjects from each line
    char  tempdistance[STRLEN], tempstart[STRLEN], tempend[STRLEN];
    FILE *infile;

    int i = 0;

    infile = fopen(filename, "r");

    if (!infile) {
        printf("Couldn't open %s for reading\n",filename);
    }
    while(i < countLines(filename) && fgets(line, sizeof(line), infile) != NULL){
        // reset variables
        strcpy(tempdistance, "");
        strcpy(tempdistance, "");
        strcpy(tempdistance, "");
        // scan one line of the file into 3 temporary variables
        sscanf(line, "%s\t%s\t%s", tempstart, tempend, tempdistance);

        if (atoi(tempdistance) == 0){
            printf("Warning: 0 distance read, %s to %s\n", tempstart, tempend);
        }

        structure[i].arraylength = countLines(filename);
        // move temp variables to structure
        structure[i].distance = atoi(tempdistance);
        strcpy(structure[i].start, tempstart);
        strcpy(structure[i].end, tempend);

        i++;
    }

    fclose(infile);
};

// Swap two nodes of heap, used in generateHeap
void swapHeapNode(HeapNode** node1,  HeapNode** node2){
     HeapNode* temp = *node1;
    *node1 = *node2;
    *node2 = temp;
}

// decrease distance value of given vertex index cityInt
// O(1) to search, O(log N) total
void updateDist(Heap* Heap, int cityInt, int distance){
    // Get the index of cityInt in  heap array
    int i = Heap->position[cityInt];

    // Get the node and update its distance value
    Heap->array[i]->distance = distance;

    // Travel up while the complete tree is not valid heap.
    // This is a O(Logn) loop
    while (i && Heap->array[i]->distance < Heap->array[(i - 1) / 2]->distance)
    {
        // Swap this node with its parent
        Heap->position[Heap->array[i]->cityInt] = (i-1)/2;
        Heap->position[Heap->array[(i-1)/2]->cityInt] = i;
        swapHeapNode(&Heap->array[i],  &Heap->array[(i - 1) / 2]);

        // move to parent index
        i = (i - 1) / 2;
    }
}

// writes an intertown_distance structure to a file
void writeFileintertown_distance(char * filename, intertown_distance structure[]) {
    FILE *infile;;
    int i = 0;

    infile = fopen(filename, "w");
    if (!infile) {
        printf("Couldn't open %s for writing\n",filename);
    }
    while(i < structure[0].arraylength){
        // scan one line of the file into 3 temporary variables
        fprintf(infile,"%s\t%s\t%d\n",structure[i].start,structure[i].end,structure[i].distance);
        i++;
    }
    fclose(infile);
    return;
};


