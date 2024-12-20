from functions import *

def main():
    dts = ['0.00010', '0.00500', '0.01000', '0.02000', '0.04000'] # Dont work for MONODOMAIN  with dx=0.0005, but work for CABLEEQ
    # dts = ['0.00050', '0.00080', '0.00100', '0.00160']
    methods = ['SSI-ADI', 'theta-ADI'] #'SSI-ADI', 'theta-ADI', 'theta-RK2' (CABLEEQ)
    thetas = ['0.50', '0.66', '1.00']

    # refs
    dx = '0.00050'
    dts = ['0.00010']
    # dts = ['0.01000']
    methods = ['theta-ADI']
    thetas = ['0.50']

    real_type = 'double'
    serial_or_gpu = 'SERIAL'
    problem = 'MONODOMAIN'
    cell_model = 'TT2' # 'AFHN', 'TT2'
    init = 'initial_conditions' # 'spiral', 'initial_conditions', 'restore_and_shift'
    frames = True
    save_last_state = True
    
    # Compile (arch=sm_80 for A100-Ampere; arch=sm_86 for RTX3050-Ampere; arch=sm_89 for RTX 4070-Ada)
    compile_command = f'nvcc -Xcompiler -fopenmp -lpthread convergence.cu -o convergence -O3 -arch={get_gpu_architecture()} -w '
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
    if init == 'spiral':
        compile_command += '-DINIT_WITH_SPIRAL '
    elif init == 'restore_and_shift':
        compile_command += '-DRESTORE_STATE_AND_SHIFT '
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
                tts = ['0.00']
            else:
                tts = thetas

            for theta in tts:

                if dt == '0.00010' and theta != '0.50':
                    continue

                os.system(f'./convergence {method} {dt} {dx} {theta}')
                plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta)
                if save_last_state:
                    plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta)
                if frames:
                    create_gif(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta)

if __name__ == '__main__':
    main()