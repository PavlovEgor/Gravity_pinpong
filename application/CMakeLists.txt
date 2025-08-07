
cmake_minimum_required(VERSION 2.8)	 # Проверка версии CMake.
										# Если версия установленной программы
										# старее указаной, произойдёт аварийный выход.

project(pinpong_in_plane)			# Название проекта

set(SOURCE_EXE pinpong_in_plane_2D.c)		# Установка переменной со списком исходников для исполняемого файла

set(SOURCE_LIB 	../src/general_ellipse_properties.c
				../src/initialize.c
				../src/step.c
				../src/saving_data.c)			# Тоже самое, но для библиотеки

add_library(PlanarPinPong STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем foo

include_directories(../include)

add_executable(pinpong_in_plane_2D ${SOURCE_EXE})	# Создает исполняемый файл с именем main

set(MYCOMPILE_FLAGS "-lm")

add_definitions(${MYCOMPILE_FLAGS})

target_link_libraries(pinpong_in_plane_2D PlanarPinPong m)		# Линковка программы с библиотекой