1) To compile: g++ -o pwbebot_example2 ~/dev/optimization/BeBOT_cpp/examples/pwbebot/pardiso/example_2/pwbebot_example2.cpp ~/dev/optimization/BeBOT_cpp/bebot/piecewisebebot.cpp ~/dev/optimization/BeBOT_cpp/bebot/piecewisebernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -ldl -lm -lstdc++

2) To include library: export LD_LIBRARY_PATH=/usr/local/lib/pardiso/panua-pardiso-20230908-linux/lib:$LD_LIBRARY_PATH
3) ./pwbebot_example2 


