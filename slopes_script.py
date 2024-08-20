import numpy as np
import matplotlib.pyplot as plt
import os

model = 'AFHN'

def run_all_simulations(method, dts, dxs, thetas, real_type, serial_or_gpu, problem):
    
    if real_type == 'float':
        double_or_float = 'USE_FLOAT'
    elif real_type == 'double':
        double_or_float = 'USE_DOUBLE'
    else:
        raise ValueError('Invalid real type')
    
    # Compile (sm_80 for A100-Ampere; sm_86 for RTX3050-Ampere; sm_89 for RTX 4070-Ada)
    compile_command = f'nvcc -Xcompiler -fopenmp -lpthread -lcusparse convergence.cu -o convergence -O3 -w -D{double_or_float} -D{serial_or_gpu} -D{problem}'
    print(f'Compiling {compile_command}...')
    os.system(compile_command)

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

def read_errors(method, dts, dxs, theta, real_type):
    errors = []
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        
        infos_path = f'./simulation_files/{real_type}/{model}/{method}/infos_{dt}_{dx}.txt'
        if method == 'theta-ADI':
            infos_path = f'./simulation_files/{real_type}/{model}/{method}/{theta}/infos_{dt}_{dx}.txt'
        
        with open(infos_path, 'r') as file:
            for line in file:
                if 'Norm-2 Error' in line:
                    line = line.split('=')
                    errors.append(float(line[-1]))
    return errors

def calculate_slopes(errors, dts):
    slopes = []
    slopes.append('-----')
    for i in range(1, len(errors)):
        slope = (np.log10(errors[i-1])-np.log10(errors[i])) / (np.log10(float(dts[i-1]))-np.log10(float(dts[i])))
        slopes.append(f'{(slope):.3f}')
    return slopes

def plot_last_frame_and_exact(method, dts, dxs, alpha, real_type, serial_or_gpu, problem):
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        
        if method != 'theta-ADI':

            # Read data from the text file
            data_last = np.genfromtxt(f'./simulation_files/{real_type}/{model}/{method}/last_{dt}_{dx}.txt', dtype=float)

            # Read data from the text file
            data_exact = np.genfromtxt(f'./simulation_files/{real_type}/{model}/{method}/exact_{dt}_{dx}.txt', dtype=float)

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
            plt.savefig(f'./simulation_files/simulation_graphs/last_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}_alpha{alpha}.png')
            plt.close()

            # Plot the exact
            plt.figure()
            plt.imshow(data_exact, cmap='viridis', vmin=-1.0, vmax=1.0)
            plt.colorbar(label='Value')
            plt.xticks([])
            plt.yticks([])
            plt.title(f'Exact dt={dt} dx={dx}')
            plt.savefig(f'./simulation_files/simulation_graphs/exact_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}_alpha{alpha}.png')
            plt.close()

def plot_exact(method, dts, dxs, alpha, real_type, serial_or_gpu, problem):
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        
        if method != 'theta-ADI':

            # Read data from the text file
            data = np.genfromtxt(f'./simulation_files/{real_type}/{model}/{method}/exact_{dt}_{dx}.txt', dtype=float)

            # Plot the data
            plt.figure()
            plt.imshow(data, cmap='viridis', vmin=0, vmax=5)
            plt.colorbar(label='Value')
            plt.xticks([])
            plt.yticks([])
            plt.title(f'Exact dt={dt} dx={dx}')
            plt.savefig(f'./simulation_files/simulation_graphs/exact_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}_alpha{alpha}.png')
            plt.close()

def plot_errors(method, dts, dxs, alpha, real_type, serial_or_gpu, problem):
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        
        if method != 'theta-ADI':

            # Read data from the text file
            data = np.genfromtxt(f'./simulation_files/{real_type}/{model}/{method}/errors_{dt}_{dx}.txt', dtype=float)

            # Plot the data
            plt.figure()
            plt.imshow(data, cmap='viridis')
            plt.colorbar(label='Value')
            plt.xticks([])
            plt.yticks([])
            plt.title(f'Errors dt={dt} dx={dx} (max_error={data.max()})')
            plt.savefig(f'./simulation_files/simulation_graphs/errors_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}_alpha{alpha}.png')
            plt.close()
    
def run_script(alpha, thetas, methods, dts, dxs, real_type="float", serial_or_gpu="SERIAL", problem="MONOAFHN"):
    
    # Create directories
    if not os.path.exists(f'./simulation_files/simulation_graphs'):
        os.makedirs(f'./simulation_files/simulation_graphs')
    if not os.path.exists(f'./simulation_files/simulation_analysis'):
        os.makedirs(f'./simulation_files/simulation_analysis')

    analysis_path = f'./simulation_files/simulation_analysis/analysis_{real_type}_{serial_or_gpu}_{problem}_alpha{alpha}.txt'
    analysis_file = open(analysis_path, 'w')

    plt.figure()
    for method in methods:
        run_all_simulations(method, dts, dxs, thetas, real_type, serial_or_gpu, problem)

        if method == 'theta-ADI':
            for theta in thetas:
                errors = read_errors(method, dts, dxs, theta, real_type)
                slopes = calculate_slopes(errors, dts)

                analysis_file.write(f'For method {method} with theta = {theta}\n')
                analysis_file.write(f'dt \t\t\t|\tdx \t\t\t|\tN-2 Error \t|\tslope\n')
                analysis_file.write('---------------------------------------------------------\n')
                print(f'For method {method} with theta = {theta}')
                print(f'dt \t\t|\tdx \t\t|\tN-2 Error \t|\tslope')
                print('------------------------------------------------------------------------------------')
                for i in range(len(errors)):
                    analysis_file.write(f'{dts[i]}\t|\t{dxs[i]}\t|\t{(errors[i]):.6f}\t|\t{slopes[i]}\n')
                    print(f'{dts[i]}\t|\t{dxs[i]}\t|\t{(errors[i]):.6f}\t|\t{slopes[i]}')
                analysis_file.write('\n\n')
                print('\n')
                
                plt.loglog([float(dt) for dt in dts], errors, '-o', label=f'{method} ({theta})')

        else:
            errors = read_errors(method, dts, dxs, '0', real_type)
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
    plt.savefig(f'./simulation_files/simulation_graphs/convergence_analysis_{real_type}_{serial_or_gpu}_{problem}_alpha{alpha}.png')
    plt.close()

    # for method in methods:
    #     plot_last_frame_and_exact(method, dts, dxs, alpha, real_type, serial_or_gpu, problem)
    #     plot_errors(method, dts, dxs, alpha, real_type, serial_or_gpu, problem)

def main():

    thetas = ['0.50', '0.66', '1.00']
    methods = ['SSI-ADI', 'theta-ADI']
    real_type = 'double'
    problem = 'MONOAFHN' # DIFF, MONOAFHN, LINMONO
    serial_or_gpu = 'GPU' # SERIAL, GPU
    serial_or_gpu = serial_or_gpu.upper()

    dts = [0.025, 0.0125, 0.00625] # Works for DIFF and MONOAFHN
    # dts = [0.01, 0.005, 0.0025, 0.00125] # Works for LINMONO
    dts[::-1].sort()
    dts = [f'{dt:.8f}' for dt in dts]
    
    alphas = [0.5] # Works for DIFF and MONOAFHN
    # alphas = [0.1]

    for a in alphas:
        # 1st order approx (dt = a*dx²)
        # 2nd order approx (dt = a*dx)
        
        dxs = [f'{(float(dt) / a):.6f}' for dt in dts]
        run_script(a, thetas, methods, dts, dxs, real_type, serial_or_gpu, problem)

if __name__ == "__main__":
    main()



