#include "quicksort.hpp"

void swap(int arr[], int i, int j) {
    int tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

void quick_(int arr[], int a, int b) {
    if (b >= a) {
        /* Select the last element as the splitter element */

        /* Partition the array into left and right */
    }
}


void quicksort(int arr[], int n) {
    quick_(arr, 0, n-1);
}

