set(handlefasta_lib_srcs FileOps.cpp FastxParser.cpp)
add_library(hfasta STATIC ${handlefasta_lib_srcs})

add_executable(handlefasta MakeFastq.cpp)
target_link_libraries(handlefasta z hfasta)
