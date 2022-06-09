/**
 *  gtests for quicksort implementation
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <random>
#include <limits>
#include <assert.h>

#include "../src/quicksort.hpp"

/* NOTE: tests include only integer type arrays since the function doesn't support any other type */

#define SIZE(array) (sizeof(array) / sizeof(array[0]))

static std::random_device dev;
static std::mt19937 rng(dev());
static std::uniform_int_distribution<std::mt19937::result_type> dist(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());

/**
 *  Compares the given arrays at every index
 */
void ArrayComparison(int *arr1, int *arr2, int size) {
  for (int i = 0; i < size; i++)
    EXPECT_EQ(arr1[i], arr2[i]) << "Arrays differ at index " << i;
}

/**
 *  Returns a random integer in range MIN_INT and MAX_INT
 */
int getRandomInt() {
  return dist(rng);
}

/** Test Type : arrays of even length */
class QuicksortTest :
  public testing::TestWithParam<int> {
};

/** Test Type: arrays of even length */
class QuicksortTestArraysEvenLength :
  public testing::TestWithParam<int> {
};

/** Test Type: arrays of odd length */
class QuicksortTestArraysOddLength :
  public testing::TestWithParam<int> {
};

/** Test Type: large arrays */
class QuicksortTestLargeArrays:
  public testing::TestWithParam<int> {
};

/**
 *  Test for empty array
 */
TEST(QuicksortTest, HandlesEmptyArray) {
  // init the array
  int arr[] = {};
  int arrSize = SIZE(arr);

  // copy the array to sort
  int arrToSort[] = {};
  std::copy(arr, arr + arrSize, arrToSort);

  // sort the array with our implementation
  quicksort(arrToSort, 0);

  ArrayComparison(arr, arrToSort, 0);
}

/**
 *  Test for array with a single element.
 *  Repeated 100 times.
 */
TEST(QuicksortTest, HandlesSingleElement) {
  // iteration number
  int n = 10000;

  // init the array
  int arr[1] = {0};
  int arrToSort[1] = {0};
  int arrSize = SIZE(arr);

  while(n > 0) {
    arr[0] = getRandomInt();

    // copy the array to sort
    std::copy(arr, arr + 1, arrToSort);

    // sort the array with our implementation
    quicksort(arrToSort, 1);

    // compare the arrays
    ArrayComparison(arr, arrToSort, 1);
    n--;
  }
}

/**
 *  Test even length arrays, e.g. arrays of length 2,4,6,100,1000 etc.
 *  Repeated 100 times.
 */
TEST_P(QuicksortTestArraysEvenLength, HandlesArraysOfEvenLength) {
  // iteration number
  int n = 100;

  // init the array
  int arrSize = GetParam();
  assert(arrSize > 0);
  // arr length must be even in this test case
  assert(arrSize % 2 == 0);

  int arr[arrSize] = {0};
  int arrToSort[arrSize] = {0};

  while(n > 0) {
    // initialize the array
    for (int i = 0; i < arrSize; i++) {
      int elem = getRandomInt();
      arr[i] = elem;
    }

    // copy the array to sort
    std::copy(arr, arr + arrSize, arrToSort);

    // sort the array with our implementation
    quicksort(arrToSort, arrSize);

    // sort the other with the library function
    std::sort(arr, arr + arrSize);

    // compare the arrays
    ArrayComparison(arr, arrToSort, arrSize);

    n--;
  }
}

TEST_P(QuicksortTestArraysOddLength, HandlesArraysOfOddLength) {
  // iteration number
  int n = 100;

  // init the array
  int arrSize = GetParam();
  assert(arrSize > 0);
  // arr length must be odd in this test case
  assert(arrSize % 2 == 1);

  int arr[arrSize] = {0};
  int arrToSort[arrSize] = {0};

  while(n > 0) {
    // initialize the array
    for (int i = 0; i < arrSize; i++) {
      int elem = getRandomInt();
      arr[i] = elem;
    }

    // copy the array to sort
    std::copy(arr, arr + arrSize, arrToSort);

    // sort the array with our implementation
    quicksort(arrToSort, arrSize);

    // sort the other with the library function
    std::sort(arr, arr + arrSize);

    // compare the arrays
    ArrayComparison(arr, arrToSort, arrSize);

    n--;
  }
}

TEST_P(QuicksortTestLargeArrays, HandlesLargeArrays) {
  // iteration number
  int n = 100;

  // init the array
  int arrSize = GetParam();
  assert(arrSize > 0);

  int arr[arrSize] = {0};
  int arrToSort[arrSize] = {0};

  while(n > 0) {
    // initialize the array
    for (int i = 0; i < arrSize; i++) {
      int elem = getRandomInt();
      arr[i] = elem;
    }

    // copy the array to sort
    std::copy(arr, arr + arrSize, arrToSort);

    // sort the array with our implementation
    quicksort(arrToSort, arrSize);

    // sort the other with the library function
    std::sort(arr, arr + arrSize);

    // compare the arrays
    ArrayComparison(arr, arrToSort, arrSize);

    n--;
  }
}

INSTANTIATE_TEST_SUITE_P(
    QuicksortTestSuit1,
    QuicksortTestArraysEvenLength,
    testing::Values(2,4,6,8,10,12,14,16,32,64,128,1048,2056,3000,10000,12422,100,98)
);

INSTANTIATE_TEST_SUITE_P(
    QuicksortTestSuit2,
    QuicksortTestArraysOddLength,
    testing::Values(1,3,5,7,9,11,13,15,31,63,127,1047,2055,325,10001,12429,101,99)
);

INSTANTIATE_TEST_SUITE_P(
    QuicksortTestSuit3,
    QuicksortTestLargeArrays,
    //testing::Values(44534, 1e6 - 1, 325123)
    testing::Values(44534, 99999, 75183)
);

