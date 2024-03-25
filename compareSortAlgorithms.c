// Marco Yglesia
// Assignment: Bonus Programming Assignment
// Language: C
// Date: Marco 20th, 2024
// Objective: You need to implement parseData, selectionSort, insertionSort, bubbleSort,
//            mergeSort and heapSort functions. Each sorting function will count the number of
//            extra memory allocated using the global variable extraMemoryAllocated.

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
    extraMemoryAllocated += sz;
    size_t* ret = malloc(sizeof(size_t) + sz);
    *ret = sz;
    printf("Extra memory allocated, size: %ld\n", sz);
    return &ret[1];
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

void DeAlloc(void* ptr)
{
    size_t* pSz = (size_t*)ptr - 1;
    extraMemoryAllocated -= *pSz;
    printf("Extra memory deallocated, size: %ld\n", *pSz);
    free(pSz);
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

size_t Size(void* ptr)
{
    return ((size_t*)ptr)[-1];
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated

void swap(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void heapify(int arr[], int n, int i) {
    // Find largest among root, left child and right child
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && arr[left] > arr[largest])
        largest = left;

    if (right < n && arr[right] > arr[largest])
        largest = right;

    // Swap and continue heapifying if root is not largest
    if (largest != i) {
        swap(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

void heapSort(int arr[], int n) {
// Build max heap
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // Heap sort
    for (int i = n - 1; i >= 0; i--) {
        swap(&arr[0], &arr[i]);

        // Heapify root element to get highest element at root again
        heapify(arr, i, 0);
    }
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

// Merges two subarrays of arr[].
// First subarray is arr[l...m]
// Second subarray is arr[m+1...r]
void merge(int arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
/* create temp arrays */
    int *L = (int*) Alloc(n1*sizeof(int));
    int *R = (int*) Alloc(n2*sizeof(int));
/* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];
/* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
/* Copy the remaining elements of L[], if there
are any */
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }
/* Copy the remaining elements of R[], if there
are any */
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
    DeAlloc(L);
    DeAlloc(R);
}

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
// Takes in an array of integers and 2 variables representing the left and right variables
void mergeSort(int pData[], int l, int r) {
    if (l < r)
    {
// get the mid point
        int m = (l+r)/2;
// Sort first and second halves
        mergeSort(pData, l, m);
        mergeSort(pData, m+1, r);
// printf("Testing l=%d r=%d m=%d\n", l, r, m);
        merge(pData, l, m, r);
    }
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------


// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n) {
    int i, item, j;
//i starts from 1 so that we can compare it with item at 0. See, j starts from i-1 in the next loop.
    for (i = 1; i < n; i++) {
        item = pData[i]; //take this item in hand
/* Move elements of arr[0..i-1], that are
greater than key, to one position ahead
of their current position */
        for(j=i-1; j>=0; j--)
        {
            if(pData[j]>item)
                pData[j+1] = pData[j]; //shift the item
            else
                break;
        }
        pData[j+1] = item;
    }
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n)
{
    int i, j;
    for (i = 0; i < n-1; i++)
// Last i elements are already in place. So, stop j before n-i
// As we will compare j with j+1, we will stop j before n-i-1
        for (j = 0; j < n-i-1; j++)
            if (pData[j] > pData[j+1])
                swap(&pData[j], &pData[j+1]);
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n) {
    int i;
    int j;
    int min_idx;
    int temp;
// One by one move boundary of unsorted subarray
    for (i = 0; i < n-1; i++) {
// Find the minimum element in unsorted array
        min_idx = i;
        for (j = i+1; j < n; j++){
            if (pData[j] < pData[min_idx]){
                min_idx = j;
            }
        }
// Swap the found minimum element with the first element
        temp = pData[i];
        pData[i] = pData[min_idx];
        pData[min_idx] = temp;}
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
    FILE* inFile = fopen(inputFileName,"r");
    int dataSz = 0;
    int i, n, *data;
    *ppData = NULL;

    if (inFile)
    {
        fscanf(inFile,"%d\n",&dataSz);
        *ppData = (int *)malloc(sizeof(int) * dataSz);
        // Implement parse data block
        if (*ppData == NULL)
        {
            printf("Cannot allocate memory\n");
            exit(-1);
        }
        for (i=0;i<dataSz;++i)
        {
            fscanf(inFile, "%d ",&n);
            data = *ppData + i;
            *data = n;
        }

        fclose(inFile);
    }

    return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
    int i, sz = (dataSz > 100 ? dataSz - 100 : 0);
    int firstHundred = (dataSz < 100 ? dataSz : 100);
    printf("\tData:\n\t");
    for (i=0;i<firstHundred;++i)
    {
        printf("%d ",pData[i]);
    }
    printf("\n\t");

    for (i=sz;i<dataSz;++i)
    {
        printf("%d ",pData[i]);
    }
    printf("\n\n");
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

int main(void) { // Start of Driver Code...
    clock_t start, end;
    int i;
    double cpu_time_used;
    char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};

    // for loop cycles the following operations for the 3 input files...
    for (i=0;i<3 ;++i) {
        int *pDataSrc, *pDataCopy;
        int dataSz = parseData(fileNames[i], &pDataSrc);

        if (dataSz <= 0)
            continue;

        pDataCopy = (int *)Alloc(sizeof(int)*dataSz);

        printf("---------------------------\n");
        printf("Dataset Size : %d\n",dataSz);
        printf("---------------------------\n");

        printf("Selection Sort:\n");
        memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
        extraMemoryAllocated = 0;
        start = clock();
        selectionSort(pDataCopy, dataSz);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
        printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
        printArray(pDataCopy, dataSz);

        printf("Insertion Sort:\n");
        memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
        extraMemoryAllocated = 0;
        start = clock();
        insertionSort(pDataCopy, dataSz);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
        printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
        printArray(pDataCopy, dataSz);

        printf("Bubble Sort:\n");
        memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
        extraMemoryAllocated = 0;
        start = clock();
        bubbleSort(pDataCopy, dataSz);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
        printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
        printArray(pDataCopy, dataSz);

        printf("Merge Sort:\n");
        memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
        extraMemoryAllocated = 0;
        start = clock();
        mergeSort(pDataCopy, 0, dataSz - 1);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
        printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
        printArray(pDataCopy, dataSz);

        printf("Heap Sort:\n");
        memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
        extraMemoryAllocated = 0;
        start = clock();
        heapSort(pDataCopy, dataSz);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
        printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
        printArray(pDataCopy, dataSz);

     }

} // End of Driver Code...