from functions import *

def main():
    dts = [0.005, 0.01, 0.02, 0.04] # Dont work for MONODOMAIN  with dx=0.0005, but work for CABLEEQ
    dts = [0.001, 0.002, 0.004, 0.005]
    methods = ['SSI-ADI', 'OS-ADI', 'FE'] #'SSI-ADI', 'theta-SSI-ADI', 'theta-RK2' (CABLEEQ), 'FE', 'OS-ADI'
    thetas = ['0.50', '0.66', '1.00']

    # refs - convergence analysis
    # dx = 0.0005
    # dy = 0.0005
    # dts = [0.0001]
    # methods = ['SSI-ADI']

    # refs - error analysis
    dx = 0.005
    dy = 0.005
    dx = 0.02 # 200 um, sugestão Berg de uso no MonoAlg
    dy = 0.02 # 200 um, sugestão Berg de uso no MonoAlg
    dts = [0.02]
    methods = ['SSI-ADI']
    
    real_type = 'double'
    serial_or_gpu = 'GPU'
    problem = 'MONODOMAIN'
    cell_model = 'MV' # 'AFHN', 'TT2', 'MV' (only in GPU by now)
    init = 'initial_conditions' # 'initial_conditions', 'restore_state'
    shift_state = False
    frames = True
    save_last_frame = True
    save_last_state = False
    measure_velocity = True
    
    for method in methods:
        for i in range(len(dts)):
            dt = dts[i]

            execution_args = f'{dt} {dx} {dy}'
            if problem == 'CABLEEQ':
                execution_args = f'{dt} {dx}'
                
            if 'theta' not in method:
                compile(real_type, serial_or_gpu, problem, cell_model, init, shift_state, frames, save_last_frame, save_last_state, measure_velocity, method)

                os.system(f'./{method} {execution_args}')
                
                plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                if save_last_state:
                    plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                if frames:
                    create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                    
            else:
                for theta in thetas:
                    compile(real_type, serial_or_gpu, problem, cell_model, init, shift_state, frames, save_last_frame, save_last_state, measure_velocity, method, theta)
                    
                    if dt == 0.0001 and theta != '0.50':
                        continue
                    
                    # os.system(f'valgrind --leak-check=full --track-origins=yes ./convergence {execution_args}')
                    os.system(f'./{method} {execution_args}')
                    
                    plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                    if save_last_state:
                        plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                    if frames:
                        create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)

if __name__ == '__main__':
    main()
