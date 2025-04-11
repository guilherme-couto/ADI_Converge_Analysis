from functions import *

def main():
    
    # Parameters
    real_type = 'double'
    serial_or_gpu = 'OPENMP'
    problem = 'MONODOMAIN'
    cell_model = 'MV' # 'AFHN', 'TT2', 'MV' (only in GPU by now)
    init = 'initial_conditions' # 'initial_conditions', 'restore_state'
    shift_state = False
    save_frames = False
    save_last_frame = True
    save_last_state = False
    measure_velocity = False

    # option 0 - vary dt | option 1 - vary dx and dy
    option = 0
    
    if option == 0:

        # dts = [0.005, 0.01, 0.02, 0.04] # Dont work for MONODOMAIN  with dx=0.0005, but work for CABLEEQ
        # thetas = ['0.50', '0.66', '1.00']

        # refs - convergence analysis
        # dx = 0.0005
        # dy = 0.0005
        # dts = [0.0001]
        # methods = ['SSI-ADI']

        # refs - error analysis
        # Referece solution config - running in cluster (A100)
        # dx = 0.002 # Ref solution
        # dy = 0.002 # Ref solution
        # dts = [0.0001] # Ref solution
        # methods = ['SSI-ADI'] # Ref solution

        dts = [0.01]
        dx = 0.01
        dy = dx
        methods = ['OS-ADI', 'SSI-ADI', 'FE']

        for method in methods:
            for i in range(len(dts)):
                dt = dts[i]

                execution_args = f'{dt} {dx} {dy}'
                if problem == 'CABLEEQ':
                    execution_args = f'{dt} {dx}'
                    
                if 'theta' not in method:
                    compile(real_type, serial_or_gpu, problem, cell_model, method, init, shift_state, save_frames, save_last_frame, save_last_state, measure_velocity)

                    os.system(f'./{method} {execution_args}')
                    
                    plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                    if save_last_state:
                        plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                    if save_frames:
                        create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                        
                else:
                    for theta in thetas:
                        compile(real_type, serial_or_gpu, problem, cell_model, method, init, shift_state, save_frames, save_last_frame, save_last_state, measure_velocity, theta)
                        
                        # os.system(f'valgrind --leak-check=full --track-origins=yes ./convergence {execution_args}')
                        os.system(f'./{method} {execution_args}')
                        
                        plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                        if save_last_state:
                            plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                        if save_frames:
                            create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
    
    # Run fixing dt and varying dx and dy
    elif option == 1:

        dxs = [0.032, 0.02, 0.016, 0.008, 0.004]
        dys = dxs
        dt = 0.001
        methods = ['OS-ADI', 'SSI-ADI', 'FE']

        for method in methods:
            for i in range(len(dxs)):
                dx = dxs[i]
                dy = dys[i]

                execution_args = f'{dt} {dx} {dy}'
                if problem == 'CABLEEQ':
                    execution_args = f'{dt} {dx}'
                    
                if 'theta' not in method:
                    compile(real_type, serial_or_gpu, problem, cell_model, method, init, shift_state, save_frames, save_last_frame, save_last_state, measure_velocity)

                    os.system(f'./{method} {execution_args}')
                    
                    plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                    if save_last_state:
                        plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                    if save_frames:
                        create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                        
                else:
                    for theta in thetas:
                        compile(real_type, serial_or_gpu, problem, cell_model, method, init, shift_state, save_frames, save_last_frame, save_last_state, measure_velocity, theta)
                        
                        os.system(f'./{method} {execution_args}')
                        
                        plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                        if save_last_state:
                            plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                        if save_frames:
                            create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)

if __name__ == '__main__':
    main()
