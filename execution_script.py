from functions import *

def main():
    dts = [0.0001, 0.005, 0.01, 0.02, 0.04] # Dont work for MONODOMAIN  with dx=0.0005, but work for CABLEEQ
    # dts = ['0.00050', '0.00080', '0.00100', '0.00160']
    methods = ['SSI-ADI', 'theta-ADI'] #'SSI-ADI', 'theta-ADI', 'theta-RK2' (CABLEEQ)
    thetas = ['0.50', '0.66', '1.00']

    # refs
    dx = 0.0005
    dy = 0.0005
    dx = 0.005
    dy = 0.005
    dts = [0.0001]
    dts = [0.005]
    methods = ['theta-ADI']
    thetas = ['0.50']

    real_type = 'double'
    serial_or_gpu = 'SERIAL'
    problem = 'MONODOMAIN'
    cell_model = 'AFHN' # 'AFHN', 'TT2'
    init = 'initial_conditions' # 'initial_conditions', 'restore_state'
    shift_state = False
    frames = True
    save_last_state = True
    
    # Compile (arch=sm_80 for A100-Ampere; arch=sm_86 for RTX3050-Ampere; arch=sm_89 for RTX 4070-Ada)
    compile_command = f'nvcc -Xcompiler -fopenmp -lpthread -lcusparse convergence.cu -o convergence -O3 -arch={get_gpu_architecture()} -w '
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
    if init == 'restore_state':
        compile_command += '-DRESTORE_STATE '
    if shift_state:
        compile_command += '-DSHIFT_STATE '
    if frames:
        compile_command += '-DSAVE_FRAMES '
    if save_last_state:
        compile_command += '-DSAVE_LAST_STATE '
    print(f'Compiling {compile_command}...')
    os.system(compile_command)
    
    for method in methods:
        for i in range(len(dts)):
            dt = dts[i]
                
            if 'theta' not in method:
                execution_args = f'{method} {dt} {dx} {dy}'
            
                if problem == 'CABLEEQ':
                    execution_args = f'{method} {dt} {dx}'
                    
                os.system(f'./convergence {execution_args}')
                
                plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                if save_last_state:
                    plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                if frames:
                    create_gif(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                    
            else:
                for theta in thetas:
                    execution_args = f'{method} {dt} {dx} {dy} {theta}' 
            
                    if problem == 'CABLEEQ':
                        execution_args = f'{method} {dt} {dx} {theta}'
                    
                    if dt == 0.0001 and theta != '0.50':
                        continue
                    
                    os.system(f'./convergence {execution_args}')
                    
                    plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                    if save_last_state:
                        plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                    if frames:
                        create_gif(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)

if __name__ == '__main__':
    main()