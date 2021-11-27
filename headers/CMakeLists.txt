set(api_files Complex.h Vectors3D.h)

install(
    FILES ${api_files}
    DESTINATION include/heinz
    COMPONENT Headers)
