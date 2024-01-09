1) To compile: g++ -o lgl_ps_example2_diffflat ~/dev/optimization/BeBOT_cpp/examples/lgl_ps/ma_57/example_2/lgl_ps_example2_diffflat.cpp ~/dev/optimization/BeBOT_cpp/bebot/lgl_ps.cpp ~/dev/optimization/BeBOT_cpp/bebot/lagrangepoly.cpp -I~/dev/optimization/BeBOT/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -ldl -lm -lstdc++

2) To include library: export LD_LIBRARY_PATH=/usr/local/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
3) ./lgl_ps_example2_diffflat

