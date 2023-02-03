# **************************************************************************** #
#                                                                              #
#                                                                              #
#    Makefile                                          Personal Website        #
#                                                     ##################       #
#    By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||       #
#                                                     ##################       #
#    Created: 2023/02/02 14:53:12 by Zian Huang                                #
#                                                                              #
# **************************************************************************** #

objects = obj/flux_func.o obj/ghost_fluid_utilities.o obj/vec_transform.o obj/gfm_2d_euler_solver.o  obj/exec_interface.o
source_files = src/flux_func.cc src/ghost_fluid_utilities.cc src/vec_transform.cc src/gfm_2d_euler_solver.cc  src/exec_interface.cc 
inline_source_files = src/inline/cell_operation.hh src/inline/debug_tools.hh src/inline/primitive_tran.hh 

run_simulation: $(objects)
	g++ $(objects) -o bin/run_simulation

obj/exec_interface.o: src/exec_interface.cc $(inline_source_files)
	g++ -c -O3 src/exec_interface.cc -o obj/exec_interface.o

obj/gfm_2d_euler_solver.o: src/gfm_2d_euler_solver.cc $(inline_source_files)
	g++ -c -O3 src/gfm_2d_euler_solver.cc -o obj/gfm_2d_euler_solver.o

obj/vec_transform.o: src/vec_transform.cc $(inline_source_files)
	g++ -c -O3 src/vec_transform.cc -o obj/vec_transform.o

obj/flux_func.o: src/flux_func.cc $(inline_source_files)
	g++ -c -O3 src/flux_func.cc -o obj/flux_func.o

obj/ghost_fluid_utilities.o: src/ghost_fluid_utilities.cc $(inline_source_files)
	g++ -c -O3 src/ghost_fluid_utilities.cc -o obj/ghost_fluid_utilities.o

clean:
	rm bin/run_simulation $(objects)

clean_data_and_log:
	rm -r data/*.dat
	rm -r data/gif/*
	rm -r debug/*
