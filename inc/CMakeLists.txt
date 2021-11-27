set(api_files heinz/Complex.h heinz/Vectors3D.h)

install(
    FILES ${api_files}
    DESTINATION include/heinz
    COMPONENT Headers)
