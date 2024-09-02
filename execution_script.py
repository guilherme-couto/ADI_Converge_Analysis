from functions import *

def main():
    dts = ['0.0050', '0.0100', '0.0200', '0.0300', '0.0400', '0.0500']
    dxs = ['0.0050', '0.0100', '0.0200']
    methods = ['SSI-ADI', 'theta-ADI']
    thetas = ['0.50', '0.66', '1.00']
    
    real_type = 'float'
    serial_or_gpu = 'GPU'
    problem = 'MONODOMAIN'
    cell_model = 'AFHN'
    init = 'spiral'
    frames = True
    
    # Compile (arch=sm_80 for A100-Ampere; arch=sm_86 for RTX3050-Ampere; arch=sm_89 for RTX 4070-Ada)
    compile_command = f'nvcc -Xcompiler -fopenmp -lpthread -lcusparse convergence.cu -o convergence -O3 -arch={get_gpu_architecture()} -w '
    if real_type == 'double':
        compile_command += '-DUSE_DOUBLE '
    elif real_type == 'float':
        compile_command += '-DUSE_FLOAT '
    if serial_or_gpu == 'GPU':
        compile_command += '-DGPU '
    elif serial_or_gpu == 'CPU':
        compile_command += '-DSERIAL '
    if problem == 'MONODOMAIN':
        compile_command += '-DMONODOMAIN '
    if cell_model == 'AFHN':
        compile_command += '-DAFHN '
    if init == 'spiral':
        compile_command += '-DINIT_WITH_SPIRAL '
    if frames:
        compile_command += '-DSAVE_FRAMES '
    print(f'Compiling {compile_command}...')
    os.system(compile_command)
    
    for method in methods:
        for i in range(len(dts)):
            dt = dts[i]
            # dx = dxs[i]
            dx = '0.0100'
            
            if method == 'SSI-ADI':
                tts = ['0.00']
            else:
                tts = thetas
            for theta in tts:
                os.system(f'./convergence {method} {dt} {dx} {theta}')
                plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta)
                create_gif(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta)

if __name__ == '__main__':
    main()