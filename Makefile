# **************************************************************************** #
#                                                                              #
#                                                                              #
#    Makefile                                          Personal Website        #
#                                                     ##################       #
#    By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||       #
#                                                     ##################       #
#    Created: 2023/01/20 19:44:03 by Zian Huang                                #
#    Updated: 2023/01/24 13:06:41 by Zian Huang                                #
#                                                                              #
# **************************************************************************** #

src_dir = src
obj_dir = obj
bin_dir = bin

objects = $(obj_dir)/flux_func.o $(obj_dir)/gfm_2d_euler_solver.o $(obj_dir)/vec_transform.o $(obj_dir)/exec_interface.o $(obj_dir)/ghost_fluid_utilities.o

all: flux_func.o gfm_2d_euler_solver.o vec_transform.o exec_interface.o ghost_fluid_utilities.o
	g++ -o $(bin_dir)/gfm_2d_euler_solver $(objects)

exec_interface.o: $(src_dir)/exec_interface.cc
	g++ -c -O3 $(src_dir)/exec_interface.cc -o $(obj_dir)/exec_interface.o

gfm_2d_euler_solver.o: $(src_dir)/gfm_2d_euler_solver.cc
	g++ -c -O3 $(src_dir)/gfm_2d_euler_solver.cc -o $(obj_dir)/gfm_2d_euler_solver.o

vec_transform.o: $(src_dir)/vec_transform.cc
	g++ -c -O3 $(src_dir)/vec_transform.cc -o $(obj_dir)/vec_transform.o

flux_func.o: $(src_dir)/flux_func.cc
	g++ -c -O3 $(src_dir)/flux_func.cc -o $(obj_dir)/flux_func.o

ghost_fluid_utilities.o: $(src_dir)/ghost_fluid_utilities.cc
	g++ -c -O3 $(src_dir)/ghost_fluid_utilities.cc -o $(obj_dir)/ghost_fluid_utilities.o

clean:
	rm $(bin_dir)/gfm_2d_euler_solver $(objects)
