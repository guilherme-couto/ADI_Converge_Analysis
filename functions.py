import numpy as np
import matplotlib.pyplot as plt
import os
import imageio.v2

cell_model = 'AFHN'

def get_gpu_architecture():
    try:
        # Run `nvidia-smi` and get the output
        stream = os.popen('nvidia-smi --query-gpu=compute_cap --format=csv,noheader')
        output = stream.read().strip()
        
        if not output:
            raise Exception("No output from nvidia-smi")

        # Get the major and minor version of the GPU architecture
        major, minor = output.split('.')
        architecture = f"sm_{major}{minor}"
        
        return architecture

    except Exception as e:
        print(f"Failed to determine GPU architecture: {e}")
        return None

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
        
        infos_path = f'./simulation_files/{real_type}/{cell_model}/{method}/infos_{dt}_{dx}.txt'
        if method == 'theta-ADI':
            infos_path = f'./simulation_files/{real_type}/{cell_model}/{method}/{theta}/infos_{dt}_{dx}.txt'
        
        if not os.path.exists(infos_path):
            raise FileNotFoundError(f'File {infos_path} not found')                                             
            
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

def plot_last_frame(method, dt, dx, real_type, serial_or_gpu, problem, theta='0.00'):
    save_dir = './simulation_files/simulation_figures'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
    
    max_value = 0.0
    min_value = 0.0
    
    if problem == 'MONOAFHN':
        max_value = 100.0
        min_value = 0.0
    
    if method != 'theta-ADI':
        # Read data from the text file
        file_path = f'./simulation_files/{real_type}/{cell_model}/{method}/last_{dt}_{dx}.txt'
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        data_last = np.genfromtxt(file_path, dtype=float)

        if max_value == 0.0 and min_value == 0.0:
            # Get the greater value to be the vmax
            max_value = data_last.max()
            # Get the minimum value to be the vmin
            min_value = data_last.min()
            
        # Plot the last
        plt.figure(figsize=(6, 6))
        plt.imshow(data_last, cmap='plasma', vmin=min_value, vmax=max_value, origin='lower')
        plt.colorbar(label='Value', fraction=0.04, pad=0.04)
        plt.xticks([])
        plt.yticks([])
        plt.title(f'Last frame {method} dt={dt} dx={dx}')
        plt.tight_layout()
        plt.savefig(f'{save_dir}/last_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}.png')
        plt.close()
        
        print(f'Last frame saved to {save_dir}/last_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}.png')
    
    else:
        # Read data from the text file
        file_path = f'./simulation_files/{real_type}/{cell_model}/{method}/{theta}/last_{dt}_{dx}.txt'
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        data_last = np.genfromtxt(file_path, dtype=float)

        # Plot the last
        plt.figure(figsize=(6, 6))
        plt.imshow(data_last, cmap='plasma', vmin=min_value, vmax=max_value, origin='lower')
        plt.colorbar(label='Value', fraction=0.04, pad=0.04)
        plt.xticks([])
        plt.yticks([])
        plt.title(f'Last frame {method} ({theta}) dt={dt} dx={dx}')
        plt.tight_layout()
        plt.savefig(f'{save_dir}/last_{dt}_{dx}_{theta}_{real_type}_{serial_or_gpu}_{problem}.png')
        plt.close()
        
        print(f'Last frame saved to {save_dir}/last_{dt}_{dx}_{theta}_{real_type}_{serial_or_gpu}_{problem}.png')

def plot_last_frame_and_exact(method, dt, dx, real_type, serial_or_gpu, problem, theta='0.00'):
    save_dir = './simulation_files/simulation_figures'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
        
    if method != 'theta-ADI':

        # Read data from the text file
        file_path = f'./simulation_files/{real_type}/{cell_model}/{method}/last_{dt}_{dx}.txt'
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        data_last = np.genfromtxt(file_path, dtype=float)

        # Read data from the text file
        file_path = f'./simulation_files/{real_type}/{cell_model}/{method}/exact_{dt}_{dx}.txt'
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        data_exact = np.genfromtxt(file_path, dtype=float)

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
        plt.imshow(data_last, cmap='viridis', vmin=-1.0, vmax=1.0, origin='lower')
        plt.colorbar(label='Value', fraction=0.04, pad=0.04)
        plt.xticks([])
        plt.yticks([])
        plt.title(f'Last dt={dt} dx={dx}')
        plt.savefig(f'{save_dir}/last_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}.png')
        plt.close()

        # Plot the exact
        plt.figure()
        plt.imshow(data_exact, cmap='viridis', vmin=-1.0, vmax=1.0, origin='lower')
        plt.colorbar(label='Value', fraction=0.04, pad=0.04)
        plt.xticks([])
        plt.yticks([])
        plt.title(f'Exact dt={dt} dx={dx}')
        plt.savefig(f'{save_dir}/exact_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}.png')
        plt.close()
        
    else:
        print("Dont have support for theta-ADI method in plot_last_frame_and_exact function yet")

def plot_exact(method, dt, dx, real_type, serial_or_gpu, problem, theta='0.00'):
    save_dir = './simulation_files/simulation_figures'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    if method != 'theta-ADI':

        # Read data from the text file
        file_path = f'./simulation_files/{real_type}/{cell_model}/{method}/exact_{dt}_{dx}.txt'
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        data = np.genfromtxt(file_path, dtype=float)

        # Plot the data
        plt.figure()
        plt.imshow(data, cmap='viridis', vmin=0, vmax=5, origin='lower')
        plt.colorbar(label='Value', fraction=0.04, pad=0.04)
        plt.xticks([])
        plt.yticks([])
        plt.title(f'Exact dt={dt} dx={dx}')
        plt.savefig(f'{save_dir}/exact_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}.png')
        plt.close()
        
    else:
        print("Dont have support for theta-ADI method in plot_exact function yet")

def plot_errors(method, dt, dx, real_type, serial_or_gpu, problem, theta='0.00'):
    save_dir = './simulation_files/simulation_figures'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    if method != 'theta-ADI':

        # Read data from the text file
        file_path = f'./simulation_files/{real_type}/{cell_model}/{method}/errors_{dt}_{dx}.txt'
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        data = np.genfromtxt(file_path, dtype=float)

        # Plot the data
        plt.figure()
        plt.imshow(data, cmap='viridis', origin='lower')
        plt.colorbar(label='Value', fraction=0.04, pad=0.04)
        plt.xticks([])
        plt.yticks([])
        plt.title(f'Errors dt={dt} dx={dx} (max_error={data.max()})')
        plt.savefig(f'{save_dir}/errors_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}.png')
        plt.close()
    
    else:
        print("Dont have support for theta-ADI method in plot_errors function yet")

def create_gif(method, dt, dx, real_type, serial_or_gpu, problem, theta='0.00'):
    # Create gif directory
    save_dir = './simulation_files/simulation_gifs'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    times = []
    frame = []
    frames = []

    frames_file = f'./simulation_files/{real_type}/{cell_model}/{method}/frames_{dt}_{dx}.txt'
    if method == 'theta-ADI':
        frames_file = f'./simulation_files/{real_type}/{cell_model}/{method}/{theta}/frames_{dt}_{dx}.txt'
    
    if not os.path.exists(frames_file):
        raise FileNotFoundError(f'File {frames_file} not found')
        
    f = open(frames_file, 'r')
    
    line = f.readline()
    line = line.split()
    frame_count = 0
    
    while True:
        if not line:
            break
        if len(line) == 1:
            times.append(float(line[0]))
            frame = []
            line = f.readline()
            if not line:
                break
            line = line.split()
        else:
            while len(line) > 1:
                # line = line.split()
                frame.append(line)
                line = f.readline()
                if not line:
                    break
                line = line.split()
                
            frame_name = f'{save_dir}/frame_{frame_count}.png'
            frames.append(frame_name)
            
            for i in range(len(frame)):
                frame[i] = [float(x) for x in frame[i]]
            
            plt.figure(figsize=(6, 6))
            if cell_model == 'AFHN':
                plt.imshow(frame, cmap='plasma', vmin=0.0, vmax=100, origin='lower')
            elif cell_model == 'TT2':
                plt.imshow(frame, cmap='plasma', vmin=-90.0, vmax=50, origin='lower')
            plt.colorbar(label='V (mV)', fraction=0.04, pad=0.04)
            if method != 'theta-ADI':
                plt.title(f'{serial_or_gpu} {problem} {cell_model} {method} dt={dt} dx={dx} ({times[frame_count]:.2f} ms)')
            else:
                plt.title(f'{serial_or_gpu} {problem} {cell_model} {method} ({theta}) dt={dt} dx={dx} ({times[frame_count]:.2f} ms)')
            plt.xticks([])
            plt.yticks([])
            plt.tight_layout()
            plt.savefig(frame_name)
            plt.close()
            
            frame_count += 1

    # Build gif
    gif_path = f'{save_dir}/gif_{dt}_{dx}_{real_type}_{serial_or_gpu}_{problem}.gif'
    if method == 'theta-ADI':
        gif_path = f'{save_dir}/gif_{dt}_{dx}_{theta}_{real_type}_{serial_or_gpu}_{problem}.gif'
    with imageio.v2.get_writer(gif_path, mode='I') as writer:
        for frame in frames:
            image = imageio.v2.imread(frame)
            writer.append_data(image)

    # Remove files
    for png in set(frames):
        if png.find('lastframe') == -1:
            os.remove(png)
            
    print(f'Gif saved to {gif_path}')
    
def run_script(alpha, thetas, methods, dts, dxs, real_type="float", serial_or_gpu="SERIAL", problem="MONOAFHN"):
    
    # Create directories
    graph_dir = f'./simulation_files/simulation_graphs'
    if not os.path.exists(graph_dir):
        os.makedirs(graph_dir)
    convergence_analysis_dir = f'./simulation_files/simulation_analysis'
    if not os.path.exists(convergence_analysis_dir):
        os.makedirs(convergence_analysis_dir)

    analysis_path = f'{convergence_analysis_dir}/convergence_analysis_{real_type}_{serial_or_gpu}_{problem}.txt'
    analysis_file = open(analysis_path, 'w')

    plt.figure()
    for method in methods:
        run_all_simulations(method, dts, dxs, thetas, real_type, serial_or_gpu, problem)

        if method == 'theta-ADI':
            for theta in thetas:
                errors = read_errors(method, dts, dxs, theta, real_type)
                slopes = calculate_slopes(errors, dts)

                analysis_file.write(f'For method {method} with theta = {theta} and alpha = {alpha}\n')
                analysis_file.write(f'dt \t\t\t|\tdx \t\t\t|\tN-2 Error \t|\tslope\n')
                analysis_file.write('---------------------------------------------------------\n')
                print(f'For method {method} with theta = {theta} and alpha = {alpha}')
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

            analysis_file.write(f'For method {method} and alpha = {alpha}\n')
            analysis_file.write(f'dt \t\t\t|\tdx \t\t\t|\tN-2 Error \t|\tslope\n')
            analysis_file.write('---------------------------------------------------------\n')
            print(f'For method {method} and alpha = {alpha}')
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
    plt.savefig(f'{graph_dir}/convergence_analysis_{real_type}_{serial_or_gpu}_{problem}.png')
    plt.close()
    
def read_values_with_rate(filename, rate):
    selected_values = []
    
    with open(filename, 'r') as file:
        
        # Read the file line by line
        for line_index, line in enumerate(file):
            
            # If the line index is multiple of the rate, select the line
            if line_index % rate == 0:
            
                values = line.split()
                values = [float(value) for value in values]
                
                # Select columns with the rate
                selected_values.extend(values[::rate])
    
    return selected_values