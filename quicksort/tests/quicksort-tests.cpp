#include <gtest/gtest.h>
#include <algorithm>

#include "../src/quicksort.hpp"

/* NOTE: tests include only integer type arrays since the function doesn't support any other type */

#define SIZE(array) (sizeof(array) / sizeof(array[0]))

void ArrayComparison(int *arr1, int *arr2, int size) {
  for (int i = 0; i < size; i++)
    EXPECT_EQ(arr1[i], arr2[i]) << "Arrays differ at index " << i;
}

TEST(QuicksortTest, HandlesEmptyArray) {
  // init the array
  int *arr = {};
  int arrSize = 0;

  // copy the array to sort
  int arrToSort[0] = {};
  std::copy(arr, arr + arrSize, arrToSort);

  // sort the array with our implementation
  quicksort(arrToSort, 0);

  ArrayComparison(arr, arrToSort, 0);
}

TEST(QuicksortTest, HandlesSingleElement) {
  // init the array
  int arr_1[1] = {15};
  int arr_2[1] = {1243214};
  int arr_3[1] = {0};
  int arr_4[1] = {-900};
  int arrSize = SIZE(arr_1);

  // copy the array to sort
  int arrToSort_1[1] = {};
  int arrToSort_2[1] = {};
  int arrToSort_3[1] = {};
  int arrToSort_4[1] = {};

  std::copy(arr_1, arr_1 + arrSize, arrToSort_1);
  std::copy(arr_2, arr_2 + arrSize, arrToSort_2);
  std::copy(arr_3, arr_3 + arrSize, arrToSort_3);
  std::copy(arr_4, arr_4 + arrSize, arrToSort_4);

  // sort the array with our implementation
  quicksort(arrToSort_1, 1);
  quicksort(arrToSort_2, 1);
  quicksort(arrToSort_3, 1);
  quicksort(arrToSort_4, 1);

  ArrayComparison(arr_1, arrToSort_1, 1);
  ArrayComparison(arr_2, arrToSort_2, 1);
  ArrayComparison(arr_3, arrToSort_3, 1);
  ArrayComparison(arr_4, arrToSort_4, 1);
}

TEST(QuicksortTest, HandlesArraysOfEvenLength) {
  // init the array
  int arr_1[2] = {12,42};
  int arr_2[4] = {-12,42,0,5413};
  int arr_3[6] = {1,2,3,4,5,6};
  int arr_4[8] = {534,91239,3,-123242,123,-4,214,0};
  int arr_5[10] = {1,-2,5,10,3,-42,23,1,-5,23};
  int arr_6[10] = {0};

  int arrSize_1 = SIZE(arr_1);
  int arrSize_2 = SIZE(arr_2);
  int arrSize_3 = SIZE(arr_3);
  int arrSize_4 = SIZE(arr_4);
  int arrSize_5 = SIZE(arr_5);
  int arrSize_6 = SIZE(arr_6);

  // copy the array to sort
  int arrToSort_1[arrSize_1] = {};
  int arrToSort_2[arrSize_2] = {};
  int arrToSort_3[arrSize_3] = {};
  int arrToSort_4[arrSize_4] = {};
  int arrToSort_5[arrSize_5] = {};
  int arrToSort_6[arrSize_6] = {};

  std::copy(arr_1, arr_1 + arrSize_1, arrToSort_1);
  std::copy(arr_2, arr_2 + arrSize_2, arrToSort_2);
  std::copy(arr_3, arr_3 + arrSize_3, arrToSort_3);
  std::copy(arr_4, arr_4 + arrSize_4, arrToSort_4);
  std::copy(arr_5, arr_5 + arrSize_5, arrToSort_5);
  std::copy(arr_6, arr_6 + arrSize_6, arrToSort_6);

  // sort the array with our implementation
  quicksort(arrToSort_1, arrSize_1);
  quicksort(arrToSort_2, arrSize_2);
  quicksort(arrToSort_3, arrSize_3);
  quicksort(arrToSort_4, arrSize_4);
  quicksort(arrToSort_5, arrSize_5);
  quicksort(arrToSort_6, arrSize_6);

  // sort the actual arrays with lib function
  std::sort(arr_1, arr_1 + arrSize_1);
  std::sort(arr_2, arr_2 + arrSize_2);
  std::sort(arr_3, arr_3 + arrSize_3);
  std::sort(arr_4, arr_4 + arrSize_4);
  std::sort(arr_5, arr_5 + arrSize_5);
  std::sort(arr_6, arr_6 + arrSize_6);

  ArrayComparison(arr_1, arrToSort_1, arrSize_1);
  ArrayComparison(arr_2, arrToSort_2, arrSize_2);
  ArrayComparison(arr_3, arrToSort_3, arrSize_3);
  ArrayComparison(arr_4, arrToSort_4, arrSize_4);
  ArrayComparison(arr_5, arrToSort_5, arrSize_5);
  ArrayComparison(arr_6, arrToSort_6, arrSize_6);
}

TEST(QuicksortTest, HandlesArraysOfOddLength) {
  // test body
}

TEST(QuicksortTest, HandlesLargeArrays) {
  // test body
}
