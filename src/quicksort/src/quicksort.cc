#include "src/quicksort/src/quicksort.h"

/**
 * Swaps two numbers given by their indices i and j
 *
 * @param i,j indices of the numbers
 */
void swap(int arr[], int i, int j) {
    int tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

int partition(int arr[], int a, int b, int splitter_index) {
    /* Swap the splitter to the end of the list */
    swap(arr, splitter_index, b);
    /* Index of the partition */
    int i = a - 1;
    for (int j = a; j < b; j++) {
        if (arr[j] < arr[b]) {
            swap(arr, ++i, j);
        }
    }
    /* Place the splitter element into the correct position */
    swap(arr, ++i, b);
    return i;
}

void _quick(int arr[], int a, int b) {
    if (a < b) {
        /* Select the last element as the splitter element */
        int splitter_index = b;

        /* Divide Step: Partition the array into left and right
         * and return the index of the splitter element. Left
         * part should be smaller and right part should be greater
         * than the splitter element. */
        int splitter_position = partition(arr, a, b, splitter_index);

        _quick(arr, a, splitter_position - 1);
        _quick(arr, splitter_position + 1, b);
    }
}

void quicksort(int arr[], int n) {
    _quick(arr, 0, n-1);
}
