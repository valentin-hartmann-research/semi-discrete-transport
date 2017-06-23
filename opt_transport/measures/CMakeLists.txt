add_library (Source Source.cpp)
add_library (Target Target.cpp)

install (TARGETS Source Target DESTINATION ${PROJECT_SOURCE_DIR}/bin)
install (FILES Source.h Target.h DESTINATION ${PROJECT_SOURCE_DIR}/include)
