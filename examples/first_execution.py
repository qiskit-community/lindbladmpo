from Lindbladian-MPO-simulator
'''
sim.build_input_file(dict_val)
sim.execute("C:/cygwin64/bin/bash.exe --login -i -c",
                      "/cygdrive/c/Users/galvz/AppData/Roaming/SPB_Data/Lindbladian-MPO-simulator",
                      "C:/Users/galvz/PycharmProjects/sim_func/newfile.txt")
'''
dict_val = {'tau': 0.2, 't_final': 2, 'input_file': "C:/Users/galvz/PycharmProjects/sim_func/newfile.txt", 'output_file': "C:/Users/galvz/AppData/Roaming/SPB_Data/Lindbladian-MPO-simulator/out"}
sim = MPOLindbladSolver(dict_val, "C:/cygwin64/bin/bash.exe --login -i -c",
                        "/cygdrive/c/Users/galvz/AppData/Roaming/SPB_Data/Lindbladian-MPO-simulator")
sim.solve()
