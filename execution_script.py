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
    # Referece solution config - running in cluster (A100)
    dx = 0.002 # Ref solution
    dy = 0.002 # Ref solution
    dts = [0.0001] # Ref solution
    methods = ['OS-ADI', 'FE'] # Ref solution

    # Berg suggestion for normal simulations
    # dx = 0.02 # 200 um, sugestão Berg de uso no MonoAlg
    # dy = 0.02 # 200 um, sugestão Berg de uso no MonoAlg
    # dts = [0.001] # Ref solution
    # methods = ['SSI-ADI'] # Ref solution

    # Simulations
    # dxs = [0.005, 0.01, 0.02, 0.04, 0.05]
    # dxs = [0.04, 0.05]
    # dys = dxs
    dts = [0.02]
    dx = 0.02
    dy = dx
    # dt = 0.01
    methods = ['FE']
    
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
                    
                    # os.system(f'valgrind --leak-check=full --track-origins=yes ./convergence {execution_args}')
                    os.system(f'./{method} {execution_args}')
                    
                    plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                    if save_last_state:
                        plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
                    if frames:
                        create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
    
    # Run fixing dx and dy and varying dt
    # for method in methods:
    #     for i in range(len(dts)):
    #         dt = dts[i]
    #         dx = 0.01
    #         dy = 0.01

    #         execution_args = f'{dt} {dx} {dy}'
    #         if problem == 'CABLEEQ':
    #             execution_args = f'{dt} {dx}'
                
    #         if 'theta' not in method:
    #             compile(real_type, serial_or_gpu, problem, cell_model, init, shift_state, frames, save_last_frame, save_last_state, measure_velocity, method)

    #             os.system(f'./{method} {execution_args}')
                
    #             plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
    #             if save_last_state:
    #                 plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
    #             if frames:
    #                 create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                    
    #         else:
    #             for theta in thetas:
    #                 compile(real_type, serial_or_gpu, problem, cell_model, init, shift_state, frames, save_last_frame, save_last_state, measure_velocity, method, theta)
                    
    #                 # os.system(f'valgrind --leak-check=full --track-origins=yes ./convergence {execution_args}')
    #                 os.system(f'./{method} {execution_args}')
                    
    #                 plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
    #                 if save_last_state:
    #                     plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
    #                 if frames:
    #                     create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
    
    # Run fixing dt and varying dx and dy
    # for method in methods:
    #     for i in range(len(dxs)):
    #         # dt = 0.005
    #         dx = dxs[i]
    #         dy = dys[i]

    #         execution_args = f'{dt} {dx} {dy}'
    #         if problem == 'CABLEEQ':
    #             execution_args = f'{dt} {dx}'
                
    #         if 'theta' not in method:
    #             compile(real_type, serial_or_gpu, problem, cell_model, init, shift_state, frames, save_last_frame, save_last_state, measure_velocity, method)

    #             os.system(f'./{method} {execution_args}')
                
    #             plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
    #             if save_last_state:
    #                 plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
    #             if frames:
    #                 create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy)
                    
    #         else:
    #             for theta in thetas:
    #                 compile(real_type, serial_or_gpu, problem, cell_model, init, shift_state, frames, save_last_frame, save_last_state, measure_velocity, method, theta)
                    
    #                 os.system(f'./{method} {execution_args}')
                    
    #                 plot_last_frame(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
    #                 if save_last_state:
    #                     plot_last_frame_state_variables(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)
    #                 if frames:
    #                     create_GIF(serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, theta)

if __name__ == '__main__':
    main()
