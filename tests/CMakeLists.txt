add_executable(test_image "test_image.cpp")
target_link_libraries(test_image PRIVATE pm)

add_executable(test_camera "test_camera.cpp")
target_link_libraries(test_camera PRIVATE pm)

add_executable(test_sphere "test_sphere.cpp")
target_link_libraries(test_sphere PRIVATE pm)

add_executable(test_plane "test_plane.cpp")
target_link_libraries(test_plane PRIVATE pm)

add_executable(test_intersector "test_intersector.cpp")
target_link_libraries(test_intersector PRIVATE pm)

add_executable(visualize_photon_map "visualize_photon_map.cpp")
target_link_libraries(visualize_photon_map PRIVATE pm)

add_executable(test_path_tracing "test_path_tracing.cpp")
target_link_libraries(test_path_tracing PRIVATE pm)