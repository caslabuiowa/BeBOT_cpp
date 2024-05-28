1) To compile: g++ -o lgl_ps_example2_diffflat ~/dev/optimization/BeBOT_cpp/examples/lgl_ps/pardiso/example_2/lgl_ps_example2_diffflat.cpp ~/dev/optimization/BeBOT_cpp/bebot/lgl_ps.cpp ~/dev/optimization/BeBOT_cpp/bebot/lagrangepoly.cpp -I~/dev/optimization/BeBOT/include -I/usr/local/include -L/usr/local/lib -lipopt -L/usr/local/lib/pardiso/panua-pardiso-20230908-linux/lib -lpardiso -ldl -lm -lstdc++

2) To include library: export LD_LIBRARY_PATH=/usr/local/lib/pardiso/panua-pardiso-20230908-linux/lib:$LD_LIBRARY_PATH
3) ./lgl_ps_example2_diffflat

