from python.MPOLindbladSolver import MPOLindbladSolver
# crate the parameters dictionary
dict_val = {'tau': 0.2, 't_final': 2,
            'input_file': "C:/Users/HaggaiLanda/Box/Haggai/python/MPO/MPO.txt",
            'output_file': "C:/Users/HaggaiLanda/Box/Haggai/python/MPO/MPO"}
# initialize class - parameters, cygwin_path, simulator_path
sim = MPOLindbladSolver(dict_val, "C:/cygwin64/bin/bash.exe",
                        "/cygdrive/c/Users/HaggaiLanda/gitprojects/Lindbladian-MPO-simulator/lindblad.exe")

# execute simulator
sim.solve()


'''
sim.build_input_file(dict_val)
sim.execute("C:/cygwin64/bin/bash.exe --login -i -c",
                      "/cygdrive/c/Users/galvz/AppData/Roaming/SPB_Data/Lindbladian-MPO-simulator",
                      "C:/Users/galvz/PycharmProjects/sim_func/newfile.txt")
'''
