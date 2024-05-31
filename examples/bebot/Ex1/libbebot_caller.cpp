#include <iostream>
#include <dlfcn.h>
#include <libgen.h>
#include <limits.h>
#include <unistd.h>
#include <cstring>

typedef void PointSetProblem;
typedef PointSetProblem* (*CreateProblemFunc)(int, double, double, double, double, double, double, double);
typedef void (*SolveProblemFunc)(PointSetProblem*);
typedef void (*GetSolutionFunc)(PointSetProblem*, double*, int);
typedef double (*GetFinalObjectiveValueFunc)(PointSetProblem*);
typedef void (*DestroyProblemFunc)(PointSetProblem*);

void* load_symbol(void* handle, const char* symbol_name) {
    void* symbol = dlsym(handle, symbol_name);
    if (!symbol) {
        std::cerr << "Cannot load symbol '" << symbol_name << "': " << dlerror() << std::endl;
    }
    return symbol;
}

void* load_library(const char* lib_name) {
    char exe_path[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
    if (len == -1) {
        std::cerr << "Cannot read the executable path: " << strerror(errno) << std::endl;
        return nullptr;
    }
    exe_path[len] = '\0';

    // Get the directory of the executable
    char* dir_path = dirname(exe_path);

    // Construct the full path to the library
    std::string lib_path = std::string(dir_path) + "/" + lib_name;

    // Load the shared library
    void* handle = dlopen(lib_path.c_str(), RTLD_LAZY);
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << std::endl;
        return nullptr;
    }

    return handle;
}

int main() {
    // Load the shared library
    void* handle = load_library("libbebot.so");
    if (!handle) {
        return 1;
    }

    // Load the symbols
    CreateProblemFunc create_problem = (CreateProblemFunc)load_symbol(handle, "create_point_set_problem");
    SolveProblemFunc solve_problem = (SolveProblemFunc)load_symbol(handle, "solve_point_set_problem");
    GetSolutionFunc get_solution = (GetSolutionFunc)load_symbol(handle, "get_solution");
    GetFinalObjectiveValueFunc get_final_objective_value = (GetFinalObjectiveValueFunc)load_symbol(handle, "get_final_objective_value");
    DestroyProblemFunc destroy_problem = (DestroyProblemFunc)load_symbol(handle, "destroy_point_set_problem");

    if (!create_problem || !solve_problem || !get_solution || !get_final_objective_value || !destroy_problem) {
        dlclose(handle);
        return 1;
    }

    // Create the problem
    int N = 7;
    double tf = 10.0;
    double amax = 5.0;
    double amin = -5.0;
    double x10 = -3.0;
    double x1f = 0.0;
    double x20 = 0.0;
    double x2f = 0.0;

    std::cout << "Sending parameters to libbebot.so:" << std::endl;
    std::cout << "N: " << N << std::endl;
    std::cout << "tf: " << tf << std::endl;
    std::cout << "amax: " << amax << std::endl;
    std::cout << "amin: " << amin << std::endl;
    std::cout << "x10: " << x10 << std::endl;
    std::cout << "x1f: " << x1f << std::endl;
    std::cout << "x20: " << x20 << std::endl;
    std::cout << "x2f: " << x2f << std::endl;

    PointSetProblem* problem = create_problem(N, tf, amax, amin, x10, x1f, x20, x2f);

    // Solve the problem
    solve_problem(problem);

    // Retrieve the solution
    double solution[N + 2];
    get_solution(problem, solution, N + 2);

    // Clean up
    dlclose(handle);

    return 0;
}