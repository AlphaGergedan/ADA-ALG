load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")

cc_library(
    name = "quicksort",
    srcs = ["//src/quicksort/src/quicksort.cc"],
    hdrs = ["src/quicksort.h"]
)

# //src/quicksort:quicksort
cc_binary(
    name = "quicksort-run",
    srcs = ["src/main.cc"],
    deps = [":quicksort"]
)
