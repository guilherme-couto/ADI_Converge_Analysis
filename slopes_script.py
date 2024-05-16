import numpy as np
import matplotlib.pyplot as plt
import os

real_type = 'double'
model = 'AFHN'

def run_all_simulations(method, dts, dxs, thetas):
    
    # Compile (sm_80 for A100-Ampere; sm_86 for RTX3050-Ampere; sm_89 for RTX 4070-Ada)
    os.system('nvcc -Xcompiler -fopenmp -lpthread -lcusparse convergence.cu -o convergence -O3 -w')

    for i in range(len(dts)):
        dx = dxs[i]
        dt = dts[i]
        tts = []
        
        if method != 'theta-ADI':
            tts = ['0.00']
        else:
            tts = thetas 
        for theta in tts:
            simulation_line = f'./convergence {method} {dt} {dx} {theta}'
            print(f'Executing {simulation_line}...')
            os.system(simulation_line)

def read_errors(method, dts, dxs, thetas):
    errors = []
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        tts = []
        
        if method != 'theta-ADI':
            with open(f'./simulation-files/{real_type}/{model}/{method}/infos-{dt}-{dx}.txt', 'r') as file:
                for line in file:
                    if 'Norm-2 Error' in line:
                        line = line.split('=')
                        errors.append(float(line[-1]))
        else:
            tts = thetas    
    return errors

def calculate_slopes(errors, dts):
    slopes = []
    slopes.append('-----')
    for i in range(1, len(errors)):
        slope = (np.log10(errors[i])-np.log10(errors[0])) / (np.log10(float(dts[i]))-np.log10(float(dts[0])))
        slopes.append(f'{(slope):.3f}')
    return slopes

def run_script(alpha):
    # 1st order approx (dt = a*dx²)
    # 2nd order approx (dt = a*dx)
    thetas = ['0.50']
    methods = ['ADI']

    # Create directories
    if not os.path.exists(f'./simulation-files/simulation-graphs'):
        os.makedirs(f'./simulation-files/simulation-graphs')
    if not os.path.exists(f'./simulation-files/simulation-analysis'):
        os.makedirs(f'./simulation-files/simulation-analysis')

    # dts = [0.00005, 0.0001, 0.0002, 0.0004, 0.0005, 0.000625]
    dts = [0.0001, 0.0002, 0.0004, 0.0005, 0.000625, 0.00125]
    # dts = [0.025, 0.0125, 0.00625, 0.00313, 0.00156]
    dts.sort()

    dts = [f'{dt:.8f}' for dt in dts]
    dxs = [f'{(float(dt) / alpha):.6f}' for dt in dts]

    analysis_path = f'./simulation-files/simulation-analysis/analysis_{alpha}.txt'
    analysis_file = open(analysis_path, 'w')

    for method in methods:
        run_all_simulations(method, dts, dxs, thetas)
        errors = read_errors(method, dts, dxs, thetas)
        slopes = calculate_slopes(errors, dts)

        analysis_file.write(f'For method {method}\n')
        analysis_file.write(f'dt \t|\t dx \t|\t N-2 Error \t|\t slope\n')
        print(f'For method {method}')
        print(f'dt \t|\t dx \t|\t N-2 Error \t|\t slope')
        for i in range(len(errors)):
            analysis_file.write(f'{dts[i]}\t|\t{dxs[i]}\t|\t{(errors[i]):.6f}\t|\t{slopes[i]}\n')
            print(f'{dts[i]}\t|\t{dxs[i]}\t|\t{(errors[i]):.6f}\t|\t{slopes[i]}')
        analysis_file.write('\n\n')
        print('\n')
        
        plt.loglog([float(dt) for dt in dts], errors, '-o', label=f'{method}')
    
    plt.xlabel('dt')
    plt.ylabel('Error')
    plt.title(f'Convergence Analysis - 2nd Order (a = {(alpha):.3f})')
    plt.legend()
    plt.savefig(f'./simulation-files/simulation-graphs/convergence-analysis_{alpha}.png')
    plt.close()


alphas = [0.01]
for a in alphas:
    run_script(a)

#dt_dx_analysis############################################################################################
# os.system('nvcc -Xcompiler -fopenmp -lpthread -lcusparse convergence.cu -o convergence -O3 -w')

# terminal_outputs = []
# dt_dx_file = open('dt_dx_analysis.txt', 'w')
# # dt_dx_file = open('dt_dx_analysis_exp.txt', 'w')
# dts = [0.01, 0.005, 0.001, 0.0005, 0.00025, 0.0001, 0.00008, 0.00005, 0.00001]

# for dt in [f'{value:.8f}' for value in dts]:
#     dxs = [0.01, 0.008, 0.00625, 0.004, 0.002, 0.001, 0.0008, 0.000625, 0.0005, 0.0004, 0.0002, 0.0001, 0.00008, 0.00005]

#     for dx in [f'{value:.6f}' for value in dxs]:
#         # simulation_line = f'./convergence theta-ADI {dt} {dx} 0.50'
#         simulation_line = f'./convergence SSI-ADI {dt} {dx} 0'
#         print(f'Executing {simulation_line}...')
#         os.system(simulation_line)
#         print('Simulation finished!\n')

#         # Save in the terminal output the value of the first element of the output file
#         # output_file = open(f'./simulation-files/double/{model}/theta-ADI/0.50/last-{dt}-{dx}.txt', 'r')
#         output_file = open(f'./simulation-files/double/{model}/SSI-ADI/last-{dt}-{dx}.txt', 'r')
#         first_element = output_file.readline().split()[0]
#         output = f'For dt = {dt} and dx = {dx}, the first element is {first_element}'
#         terminal_outputs.append(output)
#         print(output)
#         dt_dx_file.write(f'{output}\n')
#         output_file.close()
        
#         # os.system('rm -f ./simulation-files/double/{model}/theta-ADI/0.50/*.txt')
#         # os.system('rm -f ./simulation-files/double/{model}/FE/*.txt')


# # Print the terminal outputs
# for output in terminal_outputs:
#     print(output)

# exit() 
##########################################################################################################