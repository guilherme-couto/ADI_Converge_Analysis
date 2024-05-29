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
        # slope = (np.log10(errors[i])-np.log10(errors[0])) / (np.log10(float(dts[i]))-np.log10(float(dts[0])))
        slope = (np.log10(errors[i-1])-np.log10(errors[i])) / (np.log10(float(dts[i-1]))-np.log10(float(dts[i])))
        slopes.append(f'{(slope):.3f}')
    return slopes

def plot_last_frame_and_exact(method, dts, dxs, alpha):
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        
        if method != 'theta-ADI':

            # Read data from the text file
            data_last = np.genfromtxt(f'./simulation-files/{real_type}/{model}/{method}/last-{dt}-{dx}.txt', dtype=float)

            # Read data from the text file
            data_exact = np.genfromtxt(f'./simulation-files/{real_type}/{model}/{method}/exact-{dt}-{dx}.txt', dtype=float)

            # Get the greater value to be the vmax
            max_value = data_last.max()
            if data_exact.max() > max_value:
                max_value = data_exact.max()

            # Get the minimum value to be the vmin
            min_value = data_last.min()
            if data_exact.min() > min_value:
                min_value = data_exact.min()

            # Plot the last
            plt.figure()
            plt.imshow(data_last, cmap='viridis', vmin=-1.0, vmax=1.0)
            plt.colorbar(label='Value')
            plt.xticks([])
            plt.yticks([])
            plt.title(f'Last dt={dt} dx={dx}')
            plt.savefig(f'./simulation-files/simulation-graphs/last-{dt}-{dx}_{alpha}.png')
            plt.close()

            # Plot the exact
            plt.figure()
            plt.imshow(data_exact, cmap='viridis', vmin=-1.0, vmax=1.0)
            plt.colorbar(label='Value')
            plt.xticks([])
            plt.yticks([])
            plt.title(f'Exact dt={dt} dx={dx}')
            plt.savefig(f'./simulation-files/simulation-graphs/exact-{dt}-{dx}_{alpha}.png')
            plt.close()

def plot_exact(method, dts, dxs, alpha):
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        
        if method != 'theta-ADI':

            # Read data from the text file
            data = np.genfromtxt(f'./simulation-files/{real_type}/{model}/{method}/exact-{dt}-{dx}.txt', dtype=float)

            # Plot the data
            plt.figure()
            plt.imshow(data, cmap='viridis', vmin=0, vmax=5)
            plt.colorbar(label='Value')
            plt.xticks([])
            plt.yticks([])
            plt.title(f'Exact dt={dt} dx={dx}')
            plt.savefig(f'./simulation-files/simulation-graphs/exact-{dt}-{dx}_{alpha}.png')
            plt.close()

def plot_errors(method, dts, dxs, alpha):
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        
        if method != 'theta-ADI':

            # Read data from the text file
            data = np.genfromtxt(f'./simulation-files/{real_type}/{model}/{method}/errors-{dt}-{dx}.txt', dtype=float)

            # Plot the data
            plt.figure()
            plt.imshow(data, cmap='viridis')
            plt.colorbar(label='Value')
            plt.xticks([])
            plt.yticks([])
            plt.title(f'Errors dt={dt} dx={dx} (max_error={data.max()})')
            plt.savefig(f'./simulation-files/simulation-graphs/errors-{dt}-{dx}_{alpha}.png')
            plt.close()
    

def run_script(alpha):
    # 1st order approx (dt = a*dx²)
    # 2nd order approx (dt = a*dx)
    thetas = ['0.50']
    methods = ['SSI-ADI']

    # Create directories
    if not os.path.exists(f'./simulation-files/simulation-graphs'):
        os.makedirs(f'./simulation-files/simulation-graphs')
    if not os.path.exists(f'./simulation-files/simulation-analysis'):
        os.makedirs(f'./simulation-files/simulation-analysis')

    # dts = [0.025, 0.0125, 0.00625] # Works for DIFF
    dts = [0.01, 0.005, 0.0025, 0.00125] # Works for LINMONO
    dts[::-1].sort()

    dts = [f'{dt:.8f}' for dt in dts]
    dxs = [f'{(float(dt) / alpha):.6f}' for dt in dts]

    analysis_path = f'./simulation-files/simulation-analysis/analysis_{alpha}.txt'
    analysis_file = open(analysis_path, 'w')

    for method in methods:
        run_all_simulations(method, dts, dxs, thetas)
        errors = read_errors(method, dts, dxs, thetas)
        # slopes = calculate_slopes(errors, dts)
        slopes = calculate_slopes(errors, dts)

        analysis_file.write(f'For method {method}\n')
        analysis_file.write(f'dt \t\t\t|\tdx \t\t\t|\tN-2 Error \t|\tslope\n')
        analysis_file.write('---------------------------------------------------------\n')
        print(f'For method {method}')
        print(f'dt \t\t|\tdx \t\t|\tN-2 Error \t|\tslope')
        print('------------------------------------------------------------------------------------')
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

    for method in methods:
        plot_last_frame_and_exact(method, dts, dxs, alpha)
        # plot_exact(method, dts, dxs, alpha)
        plot_errors(method, dts, dxs, alpha)

def main():
    # alphas = [0.5] # Works for DIFF
    alphas = [0.1]
    for a in alphas:
        run_script(a)

if __name__ == "__main__":
    main()

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

