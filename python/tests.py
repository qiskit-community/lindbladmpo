#######################################################
##################### TEST AREA #######################
#######################################################
#from 'C:\Users\eldor\PycharmProjects\Lindbladian-MPO-simulator\python\MPOLindbladSolver.py' import MPOLindbladSolver

from python import MPOLindbladSolver

def arg_check_1():
    # this check should pass (no errors)
    input_dict = {}
    input_dict['n_qubits'] = 5
    input_dict['time'] = 20
    input_dict['time_step'] = 0.1

    input_dict['h_x'] = 1.111
    input_dict['h_y'] = [2.222, 4, 2e-2, 777777, -0.3]
    input_dict['h_z'] = 3.333

    input_dict['J_z'] = [[2, 5, 6, 7, 9],
                         [4.55, -4.1, 12, -33, 10],
                         [4.55, -1.1, 17, 0, 10],
                         [4.55, -4.1, 61, -33, 10],
                         [4.55, -1.1, 11, -33, 10]]

    input_dict['J'] = [[2, 5, 6, 7, 9],
                       [4.55, -4.1, 11, -33, 11e-9],
                       [4.55, -2.1, 11, 0, 10],
                       [6.25, -3.1, 2e-5, -33, 10],
                       [4.15, -4.1, 11, -33, 11]]

    input_dict['g_0'] = [2.222, 4, 2e-2, 777777, -0.3]
    input_dict['g_1'] = -99
    input_dict['g_2'] = -4e-10
    input_dict['init_pure_state'] = "+a"
    #    input_dict['l_x'] =
    #    input_dict['l_y'] =
    input_dict['b_periodic_x'] = True
    input_dict['b_periodic_y'] = False
    input_dict['trotter_order'] = 3
    input_dict['max_dim'] = 100
    input_dict['max_dim_rho'] = 200
    input_dict['cut_off'] = 1e-10
    input_dict['cut_off_rho'] = 1e-11
    input_dict['b_force_rho_trace'] = False
    input_dict['b_force_rho_hermitian'] = True
    input_dict['output_step'] = 1
    input_dict['save_state_file'] = "C:\\output_of_my_simulation"
    input_dict['1q_components'] = ["x", "y", "z"]
    input_dict['1q_sites'] = [1, 3, 5]
    input_dict['2q_components'] = ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]
    input_dict['2q_sites'] = [(1, 2), (3, 5), (2, 4)]

    warng_out = MPOLindbladSolver.check_argument_correctness(input_dict)
    if (warng_out == ""):
        print("Check 1 Passed")
    else:
        print("Check 1 Failed")
        print('output:\n', warng_out)


def arg_check_2():
    input_dict = {}
    input_dict['n_qubits'] = 5
    input_dict['time'] = ""  # this should fail
    input_dict['time_step'] = 0.1

    input_dict['h_x'] = 1.111
    input_dict['h_y'] = [2.222, 4, 2e-2, 777777, -0.3]
    input_dict['h_z'] = 3.333

    input_dict['J_z'] = [[2, 5, 6, 7, 9],
                         [4.55, -4.1, 12, -33, 10],
                         [4.55, -1.1, 17, 0],
                         [4.55, -4.1, 61, -33, 10],
                         [4.55, -1.1, 11, -33, 10]]  # this should fail

    input_dict['J'] = [[2, 5, 6, 7, 9],
                       [4.55, -4.1, 11, -33, 11e-9],
                       [4.55, -2.1, 11, 0, 10],
                       [6.25, -3.1, 2e-5, -33, 10]]  # this should fail

    input_dict['g_0'] = [2.222, 4, 777777, -0.3]  # this should fail
    input_dict['g_1'] = -99
    input_dict['g_2'] = "aaa"  # this should fail
    input_dict['init_pure_state'] = "+a"
    #    input_dict['l_x'] =
    #    input_dict['l_y'] =
    input_dict['b_periodic_x'] = 0  # this should fail
    input_dict['b_periodic_y'] = 1  # this should fail
    input_dict['trotter_order'] = 3
    input_dict['max_dim'] = 100
    input_dict['max_dim_rho'] = 200
    input_dict['cut_off'] = 1e-10
    input_dict['cut_off_rho'] = 1e-11
    input_dict['b_force_rho_trace'] = "False"  # this should fail
    input_dict['b_force_rho_hermitian'] = True
    input_dict['output_step'] = 0.6  # this should fail
    input_dict['save_state_file'] = "C:\\output_of_my_simulation"
    input_dict['1q_components'] = ["x", "y", "j"]  # this should fail
    input_dict['1q_sites'] = [1, 3, 5]
    input_dict['2q_components'] = ["XX", "XY", "A", "YY", "YZ", "ZZ"]  # this should fail
    input_dict['2q_sites'] = [(1, 2), (3, 5, 4), (2, 4)]  # this should fail

    warng_out = MPOLindbladSolver.check_argument_correctness(input_dict)
    if (warng_out.count("Error") == 13):
        print("Check 2 Passed")
    else:
        print("Check 2 Failed")
        print('output:\n' + warng_out)


def arg_check_3():
    input_dict = {}
    input_dict['n_qubits'] = 9
    input_dict['time'] = "22"
    input_dict['time_step'] = "aaa"

    input_dict['h_x'] = "aa1.111"
    input_dict['h_y'] = "[2.222,       4,   2e-2, 777777,   -0.3]"
    input_dict['h_z'] = "aa3.333"

    input_dict['J_z'] = [[2, 5, 6, '7', 9],
                         [4.55, -4.1, 12, -33, 10],
                         [4.55, -1.1, 17, 0],
                         [4.55, -4.1, 61, -33, 10],
                         [4.55, -1.1, 11, -33, 10]]

    input_dict['J'] = [[2, 5, '6', 7, 9],
                       [4.55, -4.1, 11, -33, 11e-9],
                       [4.55, -2.1, 11, 0, 10],
                       [6.25, -3.1, 2e-5, -33, 10]]

    input_dict['g_0'] = [2.222, 4]
    input_dict['g_1'] = (-11, 22)
    input_dict['g_2'] = "aaa"
    input_dict['init_pure_state'] = "+a2"

    input_dict['b_periodic_x'] = 0
    input_dict['b_periodic_y'] = 1
    input_dict['trotter_order'] = (3.99, 99)
    input_dict['max_dim'] = 100.99
    input_dict['max_dim_rho'] = 200.99
    input_dict['cut_off'] = "1e-10"
    input_dict['cut_off_rho'] = "1e-11"
    input_dict['b_force_rho_trace'] = "False"
    input_dict['b_force_rho_hermitian'] = 222
    input_dict['output_step'] = 0.6
    input_dict['save_state_file'] = "C:\\output_of_my_simulation"
    input_dict['1q_components'] = ["x", "y", "j"]
    input_dict['1q_sites'] = [1, 3, 1, 5]
    input_dict['2q_components'] = ["XX", "XY", "A", "YY", "YZ", "ZZ"]
    input_dict['2q_sites'] = [(1, 2), (3, 5, 4), (2, 4)]

    warng_out = MPOLindbladSolver.check_argument_correctness(input_dict)
    if (warng_out.count("Error") == 25):
        print("Check 3 Passed")
    else:
        print("Check 3 Failed")
        print('output:\n' + warng_out)


def arg_check_4():
    # this check should pass (no errors)
    input_dict = {}
    #    input_dict['n_qubits'] =
    input_dict['time'] = 20
    input_dict['time_step'] = 0.1

    input_dict['h_x'] = 1.111
    input_dict['h_y'] = [2.222, 4, 2e-2, 777777, -0.3]
    input_dict['h_z'] = 3.333

    input_dict['J_z'] = [[2, 5, 6, 7, 9],
                         [4.55, -4.1, 12, -33, 10],
                         [4.55, -1.1, 17, 0, 10],
                         [4.55, -4.1, 61, -33, 10],
                         [4.55, -1.1, 11, -33, 10]]

    input_dict['J'] = [[2, 5, 6, 7, 9],
                       [4.55, -4.1, 11, -33, 11e-9],
                       [4.55, -2.1, 11, 0, 10],
                       [6.25, -3.1, 2e-5, -33, 10],
                       [4.15, -4.1, 11, -33, 11]]

    input_dict['g_0'] = [2.222, 4, 2e-2, 777777, -0.3]
    input_dict['g_1'] = -99
    input_dict['g_2'] = -4e-10
    input_dict['init_pure_state'] = "+a"
    input_dict['l_x'] = 1
    input_dict['l_y'] = 5
    input_dict['b_periodic_x'] = True
    input_dict['b_periodic_y'] = False
    input_dict['trotter_order'] = 3
    input_dict['max_dim'] = 100
    input_dict['max_dim_rho'] = 200
    input_dict['cut_off'] = 1e-10
    input_dict['cut_off_rho'] = 1e-11
    input_dict['b_force_rho_trace'] = False
    input_dict['b_force_rho_hermitian'] = True
    input_dict['output_step'] = 1
    input_dict['save_state_file'] = "C:\\output_of_my_simulation"
    input_dict['1q_components'] = ["x", "y", "z"]
    input_dict['1q_sites'] = [1, 3, 5]
    input_dict['2q_components'] = ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]
    input_dict['2q_sites'] = [(1, 2), (3, 5), (2, 4)]

    warng_out = MPOLindbladSolver.check_argument_correctness(input_dict)
    if (warng_out == ""):
        print("Check 4 Passed")
    else:
        print("Check 4 Failed")
        print('output:\n', warng_out)


#
#   def second_check():
#       input_dict = {}
#       input_dict['n_qubits'] =
#       input_dict['time'] =
#       input_dict['time_step'] =
#       input_dict['h_x'] =
#       input_dict['h_y'] =
#       input_dict['h_z'] =
#       input_dict['J_z'] =
#       input_dict['J'] =
#       input_dict['g_0'] =
#       input_dict['g_1'] =
#       input_dict['g_2'] =
#       input_dict['init_pure_state'] =
#       input_dict['l_x'] =
#       input_dict['l_y'] =
#       input_dict['b_periodic_x'] =
#       input_dict['b_periodic_y'] =
#       input_dict['trotter_order'] =
#       input_dict['max_dim'] =
#       input_dict['max_dim_rho'] =
#       input_dict['cut_off'] =
#       input_dict['cut_off_rho'] =
#       input_dict['b_force_rho_trace'] =
#       input_dict['b_force_rho_hermitian'] =
#       input_dict['output_step'] =
#       input_dict['save_state_file'] =
#       input_dict['1q_components'] =
#       input_dict['1q_sites'] =
#       input_dict['2q_components'] =
#       input_dict['2q_sites'] =
#
#       warng_out = MPOLindbladSolver.check_argument_correctness(input_dict)
#       print ('output:',warng_out)
#
#   def second_check():
#       input_dict = {}
#       input_dict['n_qubits'] =
#       input_dict['time'] =
#       input_dict['time_step'] =
#       input_dict['h_x'] =
#       input_dict['h_y'] =
#       input_dict['h_z'] =
#       input_dict['J_z'] =
#       input_dict['J'] =
#       input_dict['g_0'] =
#       input_dict['g_1'] =
#       input_dict['g_2'] =
#       input_dict['init_pure_state'] =
#       input_dict['l_x'] =
#       input_dict['l_y'] =
#       input_dict['b_periodic_x'] =
#       input_dict['b_periodic_y'] =
#       input_dict['trotter_order'] =
#       input_dict['max_dim'] =
#       input_dict['max_dim_rho'] =
#       input_dict['cut_off'] =
#       input_dict['cut_off_rho'] =
#       input_dict['b_force_rho_trace'] =
#       input_dict['b_force_rho_hermitian'] =
#       input_dict['output_step'] =
#       input_dict['save_state_file'] =
#       input_dict['1q_components'] =
#       input_dict['1q_sites'] =
#       input_dict['2q_components'] =
#       input_dict['2q_sites'] =
#
#       warng_out = MPOLindbladSolver.check_argument_correctness(input_dict)
#       print ('output:',warng_out)
#

arg_check_1()
arg_check_2()
arg_check_3()
arg_check_4()

