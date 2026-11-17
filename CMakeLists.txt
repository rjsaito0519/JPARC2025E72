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
include_directories(${CMAKE_SOURCE_DIR}/common/include)

# ROOTを探す
find_package(ROOT REQUIRED COMPONENTS Spectrum MathMore EG Minuit Gui Core MathCore GenVector)
if(ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "ROOT not found. Please ensure ROOT is installed and accessible.")
endif()

# std::vector<TParticle> の辞書（G4 の PRM を読むツール用。実行時の GenerateDictionary は使わない）
set(G4PARTICLE_DICT_TOOLS lh2eloss_g4_samples)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
ROOT_GENERATE_DICTIONARY(G__g4particle TParticle.h LINKDEF ${CMAKE_SOURCE_DIR}/common/include/g4particle_linkdef.h NOINSTALL)

# 共有コード（common/src）を収集
file(GLOB SHARED_SRC "${CMAKE_SOURCE_DIR}/common/src/*.cpp")

# 各ソースファイルに対して実行ファイルを作成
# calibration/src と analysis/src の両方からソースを探す
file(GLOB EXECUTABLE_SOURCES 
    "${CMAKE_SOURCE_DIR}/calibration/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/analysis/src/*.cpp"
)
# "unused" という名前が含まれるファイルやディレクトリを除外
list(FILTER EXECUTABLE_SOURCES EXCLUDE REGEX "unused")
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

    # 辞書を使うツール: 辞書をリンクし、pcm を実行ファイルの隣に置く
    if(EXECUTABLE_NAME IN_LIST G4PARTICLE_DICT_TOOLS)
        target_sources(${EXECUTABLE_NAME} PRIVATE ${CMAKE_BINARY_DIR}/G__g4particle.cxx)
        add_custom_command(TARGET ${EXECUTABLE_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libg4particle_rdict.pcm ${CMAKE_SOURCE_DIR}/bin/)
    endif()

    message(STATUS "Executable created: ${EXECUTABLE_NAME}")
endforeach()

# # クリーンターゲットを追加
# add_custom_target(clean-all
#     COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
#     COMMAND ${CMAKE_COMMAND} -E echo "[success] Cleaned build directory."
# )
