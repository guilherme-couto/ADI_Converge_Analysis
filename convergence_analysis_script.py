from functions import *

def main():

    thetas = ['0.50', '0.66', '1.00']
    methods = ['SSI-ADI', 'theta-ADI']
    real_type = 'double'
    problem = 'MONODOMAIN' # DIFF, LINMONO
    cell_model = 'AFHN' # AFHN
    serial_or_gpu = 'GPU' # SERIAL, GPU
    serial_or_gpu = serial_or_gpu.upper()

    dts = [0.025, 0.0125, 0.00625] # Works for DIFF and MONODOMAIN (with AFHN)
    # dts = [0.01, 0.005, 0.0025, 0.00125] # Works for LINMONO
    dts[::-1].sort()
    dts = [f'{dt:.5f}' for dt in dts]
    
    alphas = [0.5] # Works for DIFF and MONODOMAIN (with AFHN)
    # alphas = [0.1]
    
    # Compile (arch=sm_80 for A100-Ampere; arch=sm_86 for RTX3050-Ampere; arch=sm_89 for RTX 4070-Ada)
    compile_command = f'nvcc -Xcompiler -fopenmp -lpthread -lcusparse convergence.cu -o convergence -O3 -arch={get_gpu_architecture()} -DCONVERGENCE_ANALYSIS -w '
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
    print(f'Compiling {compile_command}...')
    os.system(compile_command)

    for a in alphas:
        # 1st order approx (dt = a*dx²)
        # 2nd order approx (dt = a*dx)
        
        dxs = [f'{(float(dt) / a):.5f}' for dt in dts]
        run_script_for_convergence_analysis(a, serial_or_gpu, real_type, problem, cell_model, methods, dts, dxs, thetas)

if __name__ == "__main__":
    main()



