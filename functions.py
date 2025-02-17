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

def compile(real_type, serial_or_gpu, problem, cell_model, init, shift_state, frames, save_last_frame, save_last_state, measure_velocity, method, theta=None):
    compile_command = f'nvcc -Xcompiler -fopenmp -lpthread -lcusparse main.cu -o {method} -O3 -arch={get_gpu_architecture()} -w '
    
    if real_type == 'double':
        compile_command += '-DUSE_DOUBLE '
    elif real_type == 'float':
        compile_command += '-DUSE_FLOAT '

    if serial_or_gpu == 'GPU':
        compile_command += '-DGPU '
    elif serial_or_gpu == 'SERIAL':
        compile_command += '-DSERIAL '

    if problem == 'MONODOMAIN':
        compile_command += '-DMONODOMAIN '
    elif problem == 'CABLEEQ':
        compile_command += '-DCABLEEQ '

    if cell_model == 'AFHN':
        compile_command += '-DAFHN '
    elif cell_model == 'TT2':
        compile_command += '-DTT2 -DENDO '
    elif cell_model == 'MV':
        compile_command += '-DMV -DEPI '

    if init == 'restore_state':
        compile_command += '-DRESTORE_STATE '

    if shift_state:
        compile_command += '-DSHIFT_STATE '

    if frames:
        compile_command += '-DSAVE_FRAMES '

    if save_last_frame:
        compile_command += '-DSAVE_LAST_FRAME '

    if save_last_state:
        compile_command += '-DSAVE_LAST_STATE '
        
    if measure_velocity:
        compile_command += '-DMEASURE_VELOCITY '

    if method == 'ADI':
        compile_command += '-DADI '
    elif method == 'OS-ADI':
        compile_command += '-DOSADI '
    elif method == 'SSI-ADI':
        compile_command += '-DSSIADI '
    elif method == 'theta-SSI-ADI':
        compile_command += '-DTHETASSIADI '
        compile_command += f'-DTHETA={theta} '
    elif method == 'theta-SSI-RK2':
        compile_command += '-DTHETASSIRK2 '
        compile_command += f'-DTHETA={theta} '
    elif method == 'FE':
        compile_command += '-DFE '
    
    print(f'Compiling {compile_command}...\n')
    os.system(compile_command)

def run_all_simulations_for_convergence_analysis(serial_or_gpu, real_type, problem, cell_model, method, dts, dxs, thetas):
    
    if real_type == 'float':
        double_or_float = 'USE_FLOAT'
    elif real_type == 'double':
        double_or_float = 'USE_DOUBLE'
    else:
        raise ValueError('Invalid real type')
    
    # Compile (sm_80 for A100-Ampere; sm_86 for RTX3050-Ampere; sm_89 for RTX 4070-Ada)
    compile_command = f'nvcc -Xcompiler -fopenmp -lpthread -lcusparse convergence.cu -o convergence -O3 -w -arch={get_gpu_architecture()} -DCONVERGENCE_ANALYSIS_FORCING_TERM -D{problem} -D{serial_or_gpu} -D{double_or_float} -D{cell_model}'
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
        
        infos_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/infos/infos.txt'
        if 'theta' in method:
            infos_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/infos/infos.txt'
        
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

def plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta='0.00'):
    save_dir = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if problem == 'CABLEEQ':
        save_dir = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
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
        variables = ['Vm', 'W']
    elif cell_model == 'TT2':
        variables = ['Vm', 'X_r1', 'X_r2', 'X_s', 'm', 'h', 'j', 'd', 'f', 'f2', 'fCaSS', 's', 'r', 'Ca_i', 'Ca_SR', 'Ca_SS', 'R_prime', 'Na_i', 'K_i']
    elif cell_model == 'MV':
        variables = ['Vm', 'v', 'w', 's']
    
    for variable in variables:
        if 'theta' not in method:
            file_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe{variable}.txt'
            title = f'Last frame {variable} {method} dt={dt} dx={dx} dy={dy}'
            if problem == 'CABLEEQ':
                file_path = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe{variable}.txt'
                title = f'Last frame {variable} {method} dt={dt} dx={dx}'
        else:
            file_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe{variable}.txt'
            title = f'Last frame {variable} {method} theta={theta} dt={dt} dx={dx} dy={dy}'
            if problem == 'CABLEEQ':
                file_path = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe{variable}.txt'
                title = f'Last frame {variable} {method} theta={theta} dt={dt} dx={dx}'
            
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
            
            plt.ylim(min_value, max_value)
            plt.title(f'{title}')
            plt.ylabel('Vm (mV)')
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
        
        plt.savefig(f'{save_dir}/lastframe{variable}.png')
        plt.close()
        
        print(f'Plot of {variable} last frame saved to {save_dir}/lastframe{variable}.png')

def plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta='0.00'):
    save_dir = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if problem == 'CABLEEQ':
        save_dir = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
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
    elif cell_model == 'MV':
        max_value = 50.0
        min_value = -90.0
    
    if 'theta' not in method:
        file_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe.txt'
        title = f'Last frame {method} dt={dt} dx={dx} dy={dy}'
        if problem == 'CABLEEQ':
            file_path = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe.txt'
            title = f'Last frame {method} dt={dt} dx={dx}'
    else:
        file_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe.txt'
        title = f'Last frame {method} theta={theta} dt={dt} dx={dx} dy={dy}'
        if problem == 'CABLEEQ':
            file_path = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe.txt'
            title = f'Last frame {method} theta={theta} dt={dt} dx={dx}'
        
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
        plt.ylabel('Vm (mV)')
        plt.xlabel('L (cm)')
        plt.savefig(f'{save_dir}/lastframe.png')
        plt.close()

        # Plot AP too
        if 'theta' not in method:
            file_path = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/AP.txt'
            title = f'AP of cell 0.5 cm {method} dt={dt} dx={dx}'
        else:
            file_path = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/AP.txt'
            title = f'AP of cell 0.5 cm {method} theta={theta} dt={dt} dx={dx}'
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f'File {file_path} not found')
        
        # Read data from the text file
        data_AP = np.genfromtxt(file_path, dtype=float)

        plt.figure()
        t = [float(dt)*i for i in range(0, len(data_AP))]
        plt.plot(t, data_AP)
        plt.ylim(min_value, max_value)
        plt.title(title)
        plt.ylabel('Vm (mV)')
        plt.xlabel('t (ms)')
        plt.savefig(f'{save_dir}/AP.png')
        plt.close()
        
        os.remove(file_path)
        print(f'Action potential saved to {save_dir}/AP.png')
    
    else:
        # Plot the last
        plt.figure(figsize=(6, 6))
        plt.imshow(data_last, cmap='plasma', vmin=min_value, vmax=max_value, origin='lower')
        plt.colorbar(label='Value', fraction=0.04, pad=0.04)
        plt.xticks([])
        plt.yticks([])
        plt.title(f'{title}')
        plt.tight_layout()
        plt.savefig(f'{save_dir}/lastframe.png')
        plt.close()
    
    print(f'Plot of last frame saved to {save_dir}/lastframe.png')

def plot_last_frame_and_exact(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
        
    if 'theta' not in method:
        file_path_last = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe.txt'
        file_path_exact = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/exact/exact.txt'
        title_last = f'Last frame dt={dt} dx={dx}'
        title_exact = f'Exact dt={dt} dx={dx}'
    else:
        file_path_last = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe.txt'
        file_path_exact = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/exact/exact.txt'
        title_last = f'Last frame theta={theta} dt={dt} dx={dx}'
        title_exact = f'Exact theta={theta} dt={dt} dx={dx}'
        
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
    plt.savefig(f'{save_dir}/last.png')
    plt.close()
    
    print(f'Plot of last frame saved to {save_dir}/last.png')

    # Plot the exact
    plt.figure()
    plt.imshow(data_exact, cmap='viridis', vmin=-1.0, vmax=1.0, origin='lower')
    plt.colorbar(label='Value', fraction=0.04, pad=0.04)
    plt.xticks([])
    plt.yticks([])
    plt.title(f'{title_exact}')
    plt.savefig(f'{save_dir}/exact.png')
    plt.close()
    
    print(f'Plot of exact saved to {save_dir}/exact.png')

def plot_exact(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta='0.00'):
    save_dir = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    if 'theta' not in method:
        file_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/exact/exact.txt'
        title = f'Exact dt={dt} dx={dx}'
    else:
        file_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/exact/exact.txt'
        title = f'Exact theta={theta} dt={dt} dx={dx}'   
        
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
    plt.savefig(f'{save_dir}/exact.png')
    plt.close()
    
    print(f'Plot of exact saved to {save_dir}/exact.png')
    
def plot_errors(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta='0.00'):
    save_dir = f'./simulation_files/errors_dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    if 'theta' not in method:
        file_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/errors/errors.txt'
        title = f'Errors dt={dt} dx={dx}'
    else:
        file_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/errors/errors.txt'
        title = f'Errors theta={theta} dt={dt} dx={dx}'
        
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
    plt.savefig(f'{save_dir}/errors.png')
    plt.close()
    
    print(f'Plot of errors saved to {save_dir}/errors.png')

def plot_difference_map_from_data(data, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, Nx, Ny, theta='0.00'):
    save_dir = f'./simulation_files/difference_maps/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
    
    data = abs(data)
    data = data.reshape((Ny, Nx))
    
    if 'theta' not in method:
        title = f'DiffMap dt={dt} dx={dx} dy={dy}'
    else:
        title = f'DiffMap theta={theta} dt={dt} dx={dx} dy={dy}'
    
    # Plot the data
    plt.figure()
    plt.imshow(data, cmap='viridis', origin='lower', vmin=0, vmax=100)
    plt.colorbar(label='Value', fraction=0.04, pad=0.04)
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    plt.title(f'{title} (max_error={(data.max()):.2f})')
    plt.savefig(f'{save_dir}/diffmap.png')
    plt.close()
    
    print(f'Plot of difference map saved to {save_dir}/diffmap.png')

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
        title = f'DiffVec theta={theta} dt={dt} dx={dx}'
    
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
    plt.savefig(f'{save_dir}/diff.png')
    plt.close()
    
    print(f'Plot of difference vector saved to {save_dir}/diff.png')

def create_gif(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta='0.00'):
    # Create gif directory
    save_dir = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if problem == 'CABLEEQ':
        save_dir = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
    if 'theta' in method:
        save_dir += f'/{theta}'
    if not os.path.exists(f'{save_dir}'):
        os.makedirs(f'{save_dir}')
        
    times = []
    frame = []
    frames = []

    if 'theta' not in method:
        frames_file = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/frames.txt'
        title = f'Simulation {method} dt={dt} dx={dx} dy={dy}'
        if problem == 'CABLEEQ':
            frames_file = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/frames.txt'
            title = f'Simulation {method} dt={dt} dx={dx}'
    else:
        frames_file = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/frames.txt'
        title = f'Simulation {method} theta={theta} dt={dt} dx={dx} dy={dy}'
        if problem == 'CABLEEQ':
            frames_file = f'./simulation_files/dt_{dt}_dx_{dx}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/frames.txt'
            title = f'Simulation {method} theta={theta} dt={dt} dx={dx}'
    
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
                elif cell_model == 'MV':
                    plt.plot(x, frame[0])
                    plt.ylim(-90.0, 50.0)
                plt.title(f'{title} ({times[frame_count]:.2f} ms)')
                # plt.tight_layout()
                plt.ylabel('Vm (mV)')
                plt.xlabel('L (cm)')
                plt.savefig(frame_name)
                plt.close()
            
            else:
                plt.figure()
                if cell_model == 'AFHN':
                    plt.imshow(frame, cmap='plasma', vmin=0.0, vmax=100, origin='lower')
                elif cell_model == 'TT2':
                    plt.imshow(frame, cmap='plasma', vmin=-90.0, vmax=100, origin='lower')
                elif cell_model == 'MV':
                    plt.imshow(frame, cmap='plasma', vmin=-90.0, vmax=50.0, origin='lower')
                plt.colorbar(label='Vm (mV)', fraction=0.04, pad=0.04, orientation='horizontal')
                plt.title(f'{title} ({times[frame_count]:.2f} ms)')
                plt.xticks([])
                plt.yticks([])
                plt.tight_layout()
                plt.savefig(frame_name)
                plt.close()
            
            frame_count += 1

    # Build gif
    gif_path = f'{save_dir}/gif.gif'
    images = []
    for frame in frames:
        image = imageio.v2.imread(frame)
        images.append(image)
    imageio.v2.mimsave(
        gif_path, 
        images,
        duration=0.9,
        loop=0  # Set infinite loop
    )

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
                
                plt.loglog([float(dt) for dt in dts], errors, '-', marker='x', label=f'{method} theta={theta}')

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
            
            plt.loglog([float(dt) for dt in dts], errors, '-', marker='o', label=f'{method}')
    
    plt.xlabel('dt')
    plt.ylabel('Error')
    plt.title(f'Convergence Analysis with 2nd Order (a = {(alpha):.3f})')
    plt.legend()
    plt.savefig(f'{graph_dir}/convergence_analysis.png')
    plt.close()
    
def read_values_with_rate(filename, rate_x, rate_y):
    selected_values = []
    Nx = 0
    Ny = 0
    
    with open(filename, 'r') as file:
        
        # Read the file line by line
        for line_index, line in enumerate(file):
            
            # If the line index is multiple of the rate, select the line
            if line_index % rate_y == 0:
            
                values = line.split()
                values = [float(value) for value in values]
                
                # Select columns with the rate
                selected_values.extend(values[::rate_x])

                # Get the number of columns
                Nx = len(values[::rate_x])

                # Get the number of lines
                Ny += 1
    
    return selected_values, Nx, Ny