# Copyright (c) 2008-2016 Andrew Sutton
# All rights reserved

add_library(numeric
  concepts.cpp
  digit.cpp
  number.cpp
  fp_natural.cpp
  ap_natural.cpp
)

add_subdirectory(digit.test)
add_subdirectory(fp_natural.test)
