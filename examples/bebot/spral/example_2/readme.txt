1) To compile: g++ -o bebot_example2 ~/dev/optimization/BeBOT_cpp/examples/bebot/spral/example_2/bebot_example2.cpp ~/dev/optimization/BeBOT_cpp/bebot/bebot.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -lspral -lm -lstdc++

2) To include library: export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
3) To include library: export OMP_CANCELLATION=TRUE

4) ./bebot_example2 


