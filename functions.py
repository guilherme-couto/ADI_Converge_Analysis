import numpy as np
import matplotlib.pyplot as plt
import os
import imageio.v2

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

def run_all_simulations_for_convergence_analysis(serial_or_gpu, real_type, problem, cell_model, method, dts, dxs, thetas):
    
    if real_type == 'float':
        double_or_float = 'USE_FLOAT'
    elif real_type == 'double':
        double_or_float = 'USE_DOUBLE'
    else:
        raise ValueError('Invalid real type')
    
    # Compile (sm_80 for A100-Ampere; sm_86 for RTX3050-Ampere; sm_89 for RTX 4070-Ada)
    compile_command = f'nvcc -Xcompiler -fopenmp -lpthread -lcusparse convergence.cu -o convergence -O3 -w -arch={get_gpu_architecture()} -DCONVERGENCE_ANALYSIS -D{problem} -D{serial_or_gpu} -D{double_or_float} -D{cell_model}'
    print(f'Compiling {compile_command}...')
    os.system(compile_command)

    for i in range(len(dts)):
        dx = dxs[i]
        dt = dts[i]
        tts = []
        
        if 'theta' not in method:
            tts = ['0.00']
        else:
            tts = thetas 
        for theta in tts:
            simulation_line = f'./convergence {method} {dt} {dx} {theta}'
            print(f'Executing {simulation_line}...')
            os.system(simulation_line)

def read_errors(serial_or_gpu, real_type, problem, cell_model, method, dts, dxs, theta='0.00'):
    errors = []
    for i in range(len(dts)):
        dt = dts[i]
        dx = dxs[i]
        
        infos_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/infos/infos_{dt}_{dx}.txt'
        if 'theta' in method:
            infos_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/infos/infos_{dt}_{dx}.txt'
        
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

def plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/figures/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
    
    max_value = 0.0
    min_value = 0.0
    
    if problem == 'MONODOMAIN' and cell_model == 'AFHN':
        max_value = 100.0
        min_value = 0.0
    
    if cell_model == 'AFHN':
        variables = ['V', 'W']
    elif cell_model == 'TT2':
        variables = ['V', 'X_r1', 'X_r2', 'X_s', 'm', 'h', 'j', 'd', 'f', 'f2', 'fCaSS', 's', 'r', 'Ca_i', 'Ca_SR', 'Ca_SS', 'R_prime', 'Na_i', 'K_i']
    
    for variable in variables:
        if 'theta' not in method:
            file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe/last{variable}_{dt}_{dx}.txt'
            title = f'Last frame {variable} dt={dt} dx={dx}'
        else:
            file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe/last{variable}_{dt}_{dx}.txt'
            title = f'Last frame {variable} ({theta}) dt={dt} dx={dx}'
            
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        # Read data from the text file
        data_last = np.genfromtxt(file_path, dtype=float)

        max_value = data_last.max()
        max_value = max_value + 0.2 * abs(max_value)
        min_value = data_last.min()
        min_value = min_value - 0.2 * abs(min_value)

        if problem == 'CABLEEQ':
            plt.figure()
            # Criando o vetor x
            x = np.arange(0, (len(data_last)-1) * float(dx), float(dx))
            x = np.append(x, (len(data_last)-1) * float(dx))
            
            if cell_model == 'AFHN':
                plt.plot(x, data_last)
            elif cell_model == 'TT2':
                plt.plot(x, data_last)
            plt.ylim(min_value, max_value)
            plt.title(f'{title}')
            plt.ylabel('V (mV)')
            plt.xlabel('L (cm)')
        
        else:
            # Plot the last
            plt.figure(figsize=(6, 6))
            plt.imshow(data_last, cmap='plasma', vmin=min_value, vmax=max_value, origin='lower')
            plt.colorbar(label='Value', fraction=0.04, pad=0.04)
            plt.xticks([])
            plt.yticks([])
            plt.title(f'{title}')
            plt.tight_layout()
        
        plt.savefig(f'{save_dir}/last{variable}_{dt}_{dx}.png')
        plt.close()
        
        print(f'Last frame of {variable} saved to {save_dir}/last{variable}_{dt}_{dx}.png')

def plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/figures/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
    
    max_value = 0.0
    min_value = 0.0
    
    if cell_model == 'AFHN':
        max_value = 100.0
        min_value = 0.0
    elif cell_model == 'TT2':
        max_value = 100.0
        min_value = -90.0
    
    if 'theta' not in method:
        file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe/last_{dt}_{dx}.txt'
        title = f'Last frame dt={dt} dx={dx}'
    else:
        file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe/last_{dt}_{dx}.txt'
        title = f'Last frame ({theta}) dt={dt} dx={dx}'
        
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'File {file_path} not found')
    
    # Read data from the text file
    data_last = np.genfromtxt(file_path, dtype=float)

    if max_value == 0.0 and min_value == 0.0:
        # Get the greater value to be the vmax
        max_value = data_last.max()
        # Get the minimum value to be the vmin
        min_value = data_last.min()

    if problem == 'CABLEEQ':
        plt.figure()
        # Criando o vetor x
        x = np.arange(0, (len(data_last)-1) * float(dx), float(dx))
        x = np.append(x, (len(data_last)-1) * float(dx))
        
        plt.plot(x, data_last)
        plt.ylim(min_value, max_value)
        plt.title(f'{title}')
        plt.ylabel('V (mV)')
        plt.xlabel('L (cm)')
        plt.savefig(f'{save_dir}/last_{dt}_{dx}.png')
        plt.close()

        # Plot AP too
        if 'theta' not in method:
            file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/AP/AP_{dt}_{dx}.txt'
            title = f'AP of cell 0.5 cm dt={dt} dx={dx}'
        else:
            file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/AP/AP_{dt}_{dx}.txt'
            title = f'AP of cell 0.5 cm ({theta}) dt={dt} dx={dx}'
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        # Read data from the text file
        data_AP = np.genfromtxt(file_path, dtype=float)

        plt.figure()
        t = [float(dt)*i for i in range(0,len(data_AP))]
        plt.plot(t, data_AP)
        plt.ylim(min_value, max_value)
        plt.title(title)
        plt.ylabel('V (mV)')
        plt.xlabel('t (ms)')
        plt.savefig(f'{save_dir}/AP_{dt}_{dx}.png')
        plt.close()
    
    else:
        # Plot the last
        plt.figure(figsize=(6, 6))
        plt.imshow(data_last, cmap='plasma', vmin=min_value, vmax=max_value, origin='lower')
        plt.colorbar(label='Value', fraction=0.04, pad=0.04)
        plt.xticks([])
        plt.yticks([])
        plt.title(f'{title}')
        plt.tight_layout()
        plt.savefig(f'{save_dir}/last_{dt}_{dx}.png')
        plt.close()
    
    print(f'Last frame saved to {save_dir}/last_{dt}_{dx}.png')
    if problem == 'CABLEEQ':
        os.remove(file_path)
        print(f'Action potential saved to {save_dir}/AP_{dt}_{dx}.png')

def plot_last_frame_and_exact(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/figures/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
        
    if 'theta' not in method:
        file_path_last = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe/last_{dt}_{dx}.txt'
        file_path_exact = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/exact/exact_{dt}_{dx}.txt'
        title_last = f'Last frame dt={dt} dx={dx}'
        title_exact = f'Exact dt={dt} dx={dx}'
    else:
        file_path_last = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe/last_{dt}_{dx}.txt'
        file_path_exact = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/exact/exact_{dt}_{dx}.txt'
        title_last = f'Last frame ({theta}) dt={dt} dx={dx}'
        title_exact = f'Exact ({theta}) dt={dt} dx={dx}'
        
    if not os.path.exists(file_path_last):
        raise FileNotFoundError(f'File {file_path_last} not found')
    if not os.path.exists(file_path_exact):
        raise FileNotFoundError(f'File {file_path_exact} not found')
    
    # Read data from the text file
    data_last = np.genfromtxt(file_path_last, dtype=float)    
    data_exact = np.genfromtxt(file_path_exact, dtype=float)

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
    plt.title(f'{title_last}')
    plt.savefig(f'{save_dir}/last_{dt}_{dx}.png')
    plt.close()
    
    print(f'Last frame saved to {save_dir}/last_{dt}_{dx}.png')

    # Plot the exact
    plt.figure()
    plt.imshow(data_exact, cmap='viridis', vmin=-1.0, vmax=1.0, origin='lower')
    plt.colorbar(label='Value', fraction=0.04, pad=0.04)
    plt.xticks([])
    plt.yticks([])
    plt.title(f'{title_exact}')
    plt.savefig(f'{save_dir}/exact_{dt}_{dx}.png')
    plt.close()
    
    print(f'Exact saved to {save_dir}/exact_{dt}_{dx}.png')

def plot_exact(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/figures/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    if 'theta' not in method:
        file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/exact/exact_{dt}_{dx}.txt'
        title = f'Exact dt={dt} dx={dx}'
    else:
        file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/exact/exact_{dt}_{dx}.txt'
        title = f'Exact ({theta}) dt={dt} dx={dx}'   
        
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'File {file_path} not found')
    
    # Read data from the text file
    data = np.genfromtxt(file_path, dtype=float)

    # Plot the data
    plt.figure()
    plt.imshow(data, cmap='viridis', vmin=0, vmax=5, origin='lower')
    plt.colorbar(label='Value', fraction=0.04, pad=0.04)
    plt.xticks([])
    plt.yticks([])
    plt.title(f'{title}')
    plt.savefig(f'{save_dir}/exact_{dt}_{dx}.png')
    plt.close()
    
    print(f'Exact saved to {save_dir}/exact_{dt}_{dx}.png')
    
def plot_errors(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/errors_figures/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    if 'theta' not in method:
        file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/errors/errors_{dt}_{dx}.txt'
        title = f'Errors dt={dt} dx={dx}'
    else:
        file_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/errors/errors_{dt}_{dx}.txt'
        title = f'Errors ({theta}) dt={dt} dx={dx}'
        
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'File {file_path} not found')
    
    # Read data from the text file
    data = np.genfromtxt(file_path, dtype=float)

    # Plot the data
    plt.figure()
    plt.imshow(data, cmap='viridis', origin='lower')
    plt.colorbar(label='Value', fraction=0.04, pad=0.04)
    plt.xticks([])
    plt.yticks([])
    plt.title(f'{title} (max_error={data.max()})')
    plt.savefig(f'{save_dir}/errors_{dt}_{dx}.png')
    plt.close()
    
    print(f'Errors saved to {save_dir}/errors_{dt}_{dx}.png')

def plot_difference_map_from_data(data, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/difference_maps/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
    
    side_length = int(np.sqrt(len(data)))
    data = abs(data)
    data = data.reshape((side_length, side_length))
    
    if 'theta' not in method:
        title = f'DiffMap dt={dt} dx={dx}'
    else:
        title = f'DiffMap ({theta}) dt={dt} dx={dx}'
    
    # Plot the data
    plt.figure(figsize=(6, 6))
    plt.imshow(data, cmap='viridis', origin='lower', vmin=0, vmax=25)
    plt.colorbar(label='Value', fraction=0.04, pad=0.04)
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    plt.title(f'{title} (max_error={(data.max()):.2f})')
    plt.savefig(f'{save_dir}/diffmap_{dt}_{dx}.png')
    plt.close()
    
    print(f'Difference map saved to {save_dir}/diffmap_{dt}_{dx}.png')

def plot_difference_vector_from_data(data, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/difference_vector/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
    
    # Absolute value
    data = abs(data)
    
    if 'theta' not in method:
        title = f'DiffVec dt={dt} dx={dx}'
    else:
        title = f'DiffVec ({theta}) dt={dt} dx={dx}'
    
    # x-axis
    xaxis = [i * float(dx) for i in range(len(data))]
    
    # Plot the data
    plt.figure()
    plt.plot(xaxis, data, linestyle='-', color='b')
    plt.xlabel(f'L (cm)')
    plt.ylabel('|Diff| (mV)')
    plt.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    plt.grid()
    plt.title(f'{title} (max_error={(data.max()):.2f})')
    plt.savefig(f'{save_dir}/diff_{dt}_{dx}.png')
    plt.close()
    
    print(f'Difference vector saved to {save_dir}/diff_{dt}_{dx}.png')

def create_gif(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    # Create gif directory
    save_dir = f'./simulation_files/gifs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    times = []
    frame = []
    frames = []

    if 'theta' not in method:
        frames_file = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/frames/frames_{dt}_{dx}.txt'
        title = f'Simulation dt={dt} dx={dx}'
    else:
        frames_file = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/frames/frames_{dt}_{dx}.txt'
        title = f'Simulation ({theta}) dt={dt} dx={dx}'
    
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
            
            if problem == 'CABLEEQ':
                plt.figure()
                # Criando o vetor x
                x = np.arange(0, (len(frame[0])-1) * float(dx), float(dx))
                x = np.append(x, (len(frame[0])-1) * float(dx))
                
                if cell_model == 'AFHN':
                    # plt.imshow(frame, cmap='plasma', vmin=0.0, vmax=100, origin='lower')
                    plt.plot(x, frame[0])
                    plt.ylim(0.0, 100.0)
                elif cell_model == 'TT2':
                    plt.plot(x, frame[0])
                    plt.ylim(-90.0, 100.0)
                plt.title(f'{title} ({times[frame_count]:.2f} ms)')
                # plt.tight_layout()
                plt.ylabel('V (mV)')
                plt.xlabel('L (cm)')
                plt.savefig(frame_name)
                plt.close()
            
            else:
                plt.figure(figsize=(6, 6))
                if cell_model == 'AFHN':
                    plt.imshow(frame, cmap='plasma', vmin=0.0, vmax=100, origin='lower')
                elif cell_model == 'TT2':
                    plt.imshow(frame, cmap='plasma', vmin=-90.0, vmax=100, origin='lower')
                plt.colorbar(label='V (mV)', fraction=0.04, pad=0.04)
                plt.title(f'{title} ({times[frame_count]:.2f} ms)')
                plt.xticks([])
                plt.yticks([])
                plt.tight_layout()
                plt.savefig(frame_name)
                plt.close()
            
            frame_count += 1

    # Build gif
    gif_path = f'{save_dir}/gif_{dt}_{dx}.gif'
    with imageio.v2.get_writer(gif_path, mode='I') as writer:
        for frame in frames:
            image = imageio.v2.imread(frame)
            writer.append_data(image)

    # Remove files
    for png in set(frames):
        if png.find('lastframe') == -1:
            os.remove(png)
    
    os.remove(frames_file)
    print(f'Gif saved to {gif_path}')
      
def run_script_for_convergence_analysis(alpha, serial_or_gpu, real_type, problem, cell_model, methods, dts, dxs, thetas):
    graph_dir = f'./simulation_files/graphs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}'
    if not os.path.exists(graph_dir):
        os.makedirs(graph_dir)
        
    plt.figure()
    for method in methods:
        # Create directories
        convergence_analysis_dir = f'./simulation_files/analysis/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
        if not os.path.exists(convergence_analysis_dir):
            os.makedirs(convergence_analysis_dir)

        analysis_path = f'{convergence_analysis_dir}/convergence_analysis.txt'
        analysis_file = open(analysis_path, 'w')
        
        run_all_simulations_for_convergence_analysis(serial_or_gpu, real_type, problem, cell_model, method, dts, dxs, thetas)

        if 'theta' in method:
            for theta in thetas:
                errors = read_errors(serial_or_gpu, real_type, problem, cell_model, method, dts, dxs, theta)
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
            errors = read_errors(serial_or_gpu, real_type, problem, cell_model, method, dts, dxs)
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
    plt.title(f'Convergence Analysis with 2nd Order (a = {(alpha):.3f})')
    plt.legend()
    plt.savefig(f'{graph_dir}/convergence_analysis.png')
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