
#include <iostream>

#include "quicksort.hpp"

void print_array(int arr[], int n) {
    std::cout << "[";
    for (int i = 0; i < (n-1); i++) {
        std::cout << arr[i] << ", ";
    }
    std::cout << arr[n-1]
              << "]" << std::endl;
}

int main(int argc, char **argv) {
    int arr[] = {2, 5, 4, 5,3,2,4,5,7,12,4,6,34,65,7,55434,6,4,2,5,6,234,65356,234,6,234,6,25325,3546,23242,5,2345,432,5,25,25,234,543,52,5,46,6,65436,354,735,46};
    int n = sizeof(arr) / sizeof(arr[0]);

    std::cout << "Input: ";
    print_array(arr, n);

    quicksort(arr, n);

    std::cout << "Output: ";
    print_array(arr, n);
    return 0;
}

