# 最低限必要なCMakeのバージョンを指定
cmake_minimum_required(VERSION 3.22)

# プロジェクトの設定
project(KEKAR2023 VERSION 1.0)

# 使用するC++の標準を指定
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ビルド成果物の出力ディレクトリを`.build`以下に設定
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj)

# ヘッダーファイル用のインクルードディレクトリ
include_directories(${CMAKE_SOURCE_DIR}/include)

# ROOTを探す
find_package(ROOT REQUIRED COMPONENTS Spectrum MathMore EG Minuit Gui)
if(ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "ROOT not found. Please ensure ROOT is installed and accessible.")
endif()

# 共有コード（sharedディレクトリ）を収集
file(GLOB SHARED_SRC "${CMAKE_SOURCE_DIR}/src/shared/*.cpp")

# 各ソースファイルに対して実行ファイルを作成
file(GLOB EXECUTABLE_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
foreach(EXECUTABLE_SOURCE ${EXECUTABLE_SOURCES})
    # 実行ファイル名を取得
    get_filename_component(EXECUTABLE_NAME ${EXECUTABLE_SOURCE} NAME_WE)

    # 実行ファイルを作成し、共有コードをリンク
    add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SOURCE} ${SHARED_SRC})

    # 実行ファイルをbinディレクトリに出力
    set_target_properties(${EXECUTABLE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj
    )

    # ROOTライブラリをリンク
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${ROOT_LIBRARIES})

    message(STATUS "Executable created: ${EXECUTABLE_NAME}")
endforeach()

# # クリーンターゲットを追加
# add_custom_target(clean-all
#     COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
#     COMMAND ${CMAKE_COMMAND} -E echo "[success] Cleaned build directory."
# )
