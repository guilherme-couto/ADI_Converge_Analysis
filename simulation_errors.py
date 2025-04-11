import numpy as np
import os
from functions import *

def main():

    methods = ['FE', 'OS-ADI', 'SSI-ADI'] #'SSI-ADI', 'theta-SSI-ADI', 'SSI-CN' (CABLEEQ), 'theta-RK2' (CABLEEQ), 'FE', 'OS-ADI'
    thetas = ['0.50', '0.66', '1.00']

    real_type = 'double'
    serial_or_gpu = 'GPU'
    problem = 'MONODOMAIN'
    cell_model = 'MV' # 'AFHN', 'TT2', 'MV'
    
    # Find the largest dx to use as base for the read rate
    # base_dx = 0
    # for dx in dxs:
    #     if float(dx) > base_dx:
    #         base_dx = float(dx)
    # base_dx = 0.01

    # 0 for varying dt, 1 for varying dx and dy
    option = 1

    if option == 0:

        # dts = [0.005, 0.01, 0.02] # Dont work for MONODOMAIN  with dx=0.0005, but work for CABLEEQ
        dts = [0.002, 0.003, 0.004, 0.005, 0.01, 0.02, 0.04, 0.05, 0.08, 0.16, 0.2, 0.32, 0.5, 0.64, 0.8]

        dx = reference_dx
        dy = reference_dy
        
        # Read reference solution
        reference_dx = 0.002
        reference_dy = 0.002
        reference_solution_path = f'./reference_solutions/{real_type}/{problem}/{cell_model}/lastframe.txt'
        
        if not os.path.exists(reference_solution_path):
            raise FileNotFoundError(f'Reference solution not found at {reference_solution_path}')
        
        print(f'Reading reference solution from {reference_solution_path}')
        reference_data, Nx, Ny = read_values_with_rate(reference_solution_path, 1, 1)
        print(f'Reference solution read successfully. Total size: {len(reference_data)}')
        print()

        for method in methods:
            # Create error analysis file
            error_analysis_dir = f'./simulation_files/analysis_dt/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
            if not os.path.exists(error_analysis_dir):
                os.makedirs(error_analysis_dir)
                
            error_analysis_path = f'{error_analysis_dir}/error_analysis.txt'
            ea_file = open(error_analysis_path, 'w')
            
            # Prepare error analysis file
            ea_file.write(f'For method {method}\n')
            
            if 'theta' not in method:
                ea_file.write(f'dt (ms)\t| dx (cm)\t| dy (cm)\t| N-2 Error\t\t| slope\t\t| REP (%)\n')
                ea_file.write('---------------------------------------------------------------------------------\n')
                
                # Comparative plot
                plt.figure()
                plt.plot(reference_data, label='ref')

                # Initialize lists to store dt and n2_error values
                n2_errors = []

                # Iterate over the dts
                prev_error = 0
                for i in range(len(dts)):
                    dt = dts[i]
                                        
                    simulation_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe.txt'
                    if not os.path.exists(simulation_path):
                        raise FileNotFoundError(f'Simulation file not found at {simulation_path}')
                    
                    # Read simulation file
                    print(f'Reading files with rate 1 for x and 1 for y')
                    print(f'Reading simulation file from {simulation_path}')
                    simulation_data, Nx, Ny = read_values_with_rate(simulation_path, 1, 1)
                    print(f'Simulation file read successfully. Total size: {len(simulation_data)}')
                    
                    # Plot data to compare
                    plt.plot(simulation_data, label=f'dt={dt}')
                    
                    # Plot difference map between reference and simulation
                    difference = np.array(simulation_data) - np.array(reference_data)
                    if problem != 'CABLEEQ':
                        plot_difference_map_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, Nx, Ny)
                    else:
                        plot_difference_vector_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx)
                    
                    # Calculate the N-2 error
                    n2_error = np.linalg.norm(difference, 2) * np.sqrt(1/(len(difference)))
                    print(f'N-2 error: {n2_error}')
                    
                    # Calculate the convergence rate (slope)
                    slope = '-----'
                    if i > 0:
                        convergence_rate = (np.log10(n2_error)-np.log10(prev_error)) / (np.log10(dt)-np.log10(dts[i-1]))
                        slope = f'{(convergence_rate):.3f}'
                    print(f'Convergence rate: {slope}')
                    print()

                    # Calculate the relative error percentage
                    sum_abs_diff = np.sum(np.abs(difference))
                    sum_abs_ref = np.sum(np.abs(reference_data))
                    relative_error = sum_abs_diff / sum_abs_ref
                    relative_error_percentage = relative_error * 100
                    print(f'Relative error percentage: {relative_error_percentage:.2f}%')
                    
                    # Write to error analysis file
                    ea_file.write(f'{dt} \t| {dx} \t| {dy} \t| {(n2_error):.3e} \t| {slope} \t| {relative_error_percentage:.2f}\n')
                    
                    prev_error = n2_error

                    # Store dt and n2_error for later plotting
                    n2_errors.append(n2_error)
                    
                ea_file.write('\n')

                plt.grid()
                plt.legend()
                plt.title(f'Last Frame Comparative for {method}')
                plt.savefig(f'{error_analysis_dir}/lastframe_comparative.png')
                plt.close()

                # Convert dt_values and n2_errors to numpy arrays
                dt_values = np.array(dts)
                n2_errors = np.array(n2_errors)

                # Fit a line to the log-log data (least squares method)
                if np.isnan(n2_errors).all():
                    ea_file.write(f'Slope of the fitted line (least squares in log-log): NaN\n\n')
                else:
                    # Calculate the logarithms of dt_values and n2_errors
                    log_dt_values = np.log10(dt_values)
                    log_n2_errors = np.log10(n2_errors)

                    try:
                        coefficients = np.polyfit(log_dt_values, log_n2_errors, 1)  # 1 indicates a linear fit
                        line_slope = coefficients[0]

                        # Create the linear fit function in log-log space
                        linear_fit = np.poly1d(coefficients)

                        # Plot the data of N2-error vs dt in log-log scale
                        plt.figure()
                        plt.loglog(dt_values, n2_errors, 'o', color='blue', label='N-2 Error')
                        plt.loglog(dt_values, 10**linear_fit(log_dt_values), color='red', label='Linear Fit', linestyle='--')
                        plt.xlabel('dt')
                        plt.ylabel('N-2 Error')
                        plt.title(f'N-2 Error vs dt for {method}')
                        plt.legend()
                        plt.grid()
                        plt.savefig(f'{error_analysis_dir}/n2_error_vs_dt_loglog.png')
                        plt.close()

                        # Write the slope of the fitted line to the file
                        ea_file.write(f'Slope of the fitted line (least squares in log-log): {line_slope:.6f}\n\n')
                    
                    except Exception as e:
                        print(f'Error during linear fit: {e}')
                        ea_file.write(f'Error during linear fit: {e}\n\n')
            
            else:
                thetas_linear_fits = []
                for theta in thetas:
                    ea_file.write(f'with theta={theta}\n')
                    ea_file.write(f'dt (ms)\t| dx (cm)\t| dy (cm)\t| N-2 Error\t\t| slope\t\t| REP (%)\n')
                    ea_file.write('-----------------------------------------------------------------------------------\n')
                    
                    # Comparative plot
                    if problem != 'MONODOMAIN':
                        plt.figure()
                        plt.plot(reference_data, label='ref')

                    # Initialize lists to store dt and n2_error values
                    n2_errors = []

                    # Iterate over the dts and dxs
                    prev_error = 0
                    for i in range(len(dts)):
                        dt = dts[i]
                                            
                        simulation_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe.txt'
                        if not os.path.exists(simulation_path):
                            raise FileNotFoundError(f'Simulation file not found at {simulation_path}')
                        
                        # Read simulation file
                        print(f'Reading simulation file from {simulation_path}')
                        simulation_data, Nx, Ny = read_values_with_rate(simulation_path, 1, 1)
                        print(f'Simulation file read successfully. Total size: {len(simulation_data)}')                        

                        # Plot difference map between reference and simulation
                        difference = np.array(simulation_data) - np.array(reference_data)
                        if problem != 'CABLEEQ':
                            plot_difference_map_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, Nx, Ny, theta)
                        else:
                            # Plot data to compare
                            plt.plot(simulation_data, label=f'dt={dt}')
                            plot_difference_vector_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta)               
                        
                        # Calculate the N-2 error
                        n2_error = np.linalg.norm(difference, 2) * np.sqrt(1/(len(difference)))
                        print(f'N-2 error: {n2_error}')
                        
                        # Calculate the convergence rate (slope)
                        slope = '-----'
                        if i > 0:
                            convergence_rate = (np.log10(n2_error)-np.log10(prev_error)) / (np.log10(dt)-np.log10(dts[i-1]))
                            slope = f'{(convergence_rate):.3f}'
                        print(f'Convergence rate: {slope}')
                        print()

                        # Calculate the relative error percentage
                        sum_abs_diff = np.sum(np.abs(difference))
                        sum_abs_ref = np.sum(np.abs(reference_data))
                        relative_error = sum_abs_diff / sum_abs_ref
                        relative_error_percentage = relative_error * 100
                        print(f'Relative error percentage: {relative_error_percentage:.2f}%')
                        
                        # Write to error analysis file
                        ea_file.write(f'{dt} \t| {dx} \t| {dy} \t| {(n2_error):.3e} \t| {slope} \t| {relative_error_percentage:.2f}\n')
                        
                        prev_error = n2_error

                        # Store dt and n2_error for later plotting
                        n2_errors.append(n2_error)
                        
                    ea_file.write('\n')

                    if problem != 'MONODOMAIN':
                        plt.grid()
                        plt.legend()
                        plt.title(f'Last Frame Comparative for {method} (theta={theta})')
                        plt.savefig(f'{error_analysis_dir}/lastframe_comparative_{theta}.png')
                        plt.close()

                    # Convert dt_values and n2_errors to numpy arrays
                    dt_values = np.array(dts)
                    n2_errors = np.array(n2_errors)

                    # Fit a line to the log-log data (least squares method)
                    if np.isnan(n2_errors).all():
                        ea_file.write(f'Slope of the fitted line (least squares in log-log): NaN\n\n')
                    else:
                        # Calculate the logarithms of dt_values and n2_errors
                        log_dt_values = np.log10(dt_values)
                        log_n2_errors = np.log10(n2_errors)
                        
                        try:

                            # Fit a line to the log-log data (least squares method)
                            coefficients = np.polyfit(log_dt_values, log_n2_errors, 1)  # 1 indicates a linear fit
                            line_slope = coefficients[0]

                            # Create the linear fit function in log-log space
                            linear_fit = np.poly1d(coefficients)

                            # Plot the data of N2-error vs dt in log-log scale
                            plt.figure()
                            plt.loglog(dt_values, n2_errors, 'o', color='blue', label='N-2 Error')
                            plt.loglog(dt_values, 10**linear_fit(log_dt_values), color='red', label='Linear Fit', linestyle='--')
                            plt.xlabel('dt')
                            plt.ylabel('N-2 Error')
                            plt.title(f'N-2 Error vs dt for {method} (theta={theta})')
                            plt.legend()
                            plt.grid()
                            plt.savefig(f'{error_analysis_dir}/n2_error_vs_dt_loglog_{theta}.png')
                            plt.close()

                            thetas_linear_fits.append(10**linear_fit(log_dt_values))

                            # Write the slope of the fitted line to the file
                            ea_file.write(f'Slope of the fitted line (least squares in log-log): {line_slope:.6f}\n\n')
                            
                        except Exception as e:
                            print(f'Error during linear fit: {e}')
                            ea_file.write(f'Error during linear fit: {e}\n\n')
                
                # plt.figure()
                # for t in range(len(thetas)):
                #     plt.loglog(dt_values, thetas_linear_fits[t], label=f'theta={thetas[t]}', linestyle='-')
                # plt.xlabel('dt')
                # plt.ylabel('N-2 Error')
                # plt.title(f'N-2 Error vs dt for {method}')
                # plt.legend()
                # plt.grid()
                # plt.savefig(f'{error_analysis_dir}/n2_error_vs_dt_loglog.png')
                # plt.close()

    elif option == 1:

        dxs = [0.004, 0.008, 0.016, 0.02]
        dys = dxs

        reference_dx = 0.002
        reference_dy = 0.002

        dt = 0.001

        for method in methods:
        
            # Create error analysis file
            error_analysis_dir = f'./simulation_files/analysis_dx/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
            if not os.path.exists(error_analysis_dir):
                os.makedirs(error_analysis_dir)
                
            error_analysis_path = f'{error_analysis_dir}/error_analysis.txt'
            ea_file = open(error_analysis_path, 'w')
            
            # Prepare error analysis file
            ea_file.write(f'For method {method}\n') 
            ea_file.write(f'dt (ms)\t| dx (cm)\t| dy (cm)\t| N-2 Error\t\t| slope\t\t| REP (%)\n')
            ea_file.write('---------------------------------------------------------------------------------\n')

            # Initialize lists to store dt and n2_error values
            n2_errors = []

            # Iterate over the dts
            prev_error = 0
            for i in range(len(dxs)):
                dx = dxs[i]
                dy = dys[i]

                rate_x = int(dx / reference_dx)
                rate_y = int(dy / reference_dy)

                # Read reference solution
                reference_dx = 0.002
                reference_dy = 0.002
                reference_solution_path = f'./reference_solutions/{real_type}/{problem}/{cell_model}/lastframe.txt'
                
                if not os.path.exists(reference_solution_path):
                    raise FileNotFoundError(f'Reference solution not found at {reference_solution_path}')
                
                print(f'Reading reference solution from {reference_solution_path} with rate {rate_x} for x and {rate_y} for y')
                reference_data, Nx, Ny = read_values_with_rate(reference_solution_path, rate_x, rate_y)
                print(f'Reference solution read successfully. Total size: {len(reference_data)}')
                print()
                
                print(f'Reading simulation files with rate 1 for x and 1 for y')

                simulation_path = f'./simulation_files/dt_{dt}_dx_{dx}_dy_{dy}/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe.txt'
                if not os.path.exists(simulation_path):
                    raise FileNotFoundError(f'Simulation file not found at {simulation_path}')
                
                # Read simulation file
                print(f'Reading simulation file from {simulation_path}')
                simulation_data, Nx, Ny = read_values_with_rate(simulation_path, 1, 1)
                print(f'Simulation file read successfully. Total size: {len(simulation_data)}')
                
                # Plot data to compare
                plt.plot(simulation_data, label=f'dx={dx} dy={dy}')
                
                # Plot difference map between reference and simulation
                difference = np.array(simulation_data) - np.array(reference_data)
                if problem != 'CABLEEQ':
                    plot_difference_map_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, dy, Nx, Ny)
                else:
                    plot_difference_vector_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx)

                # Calculate the N-2 error
                n2_error = np.linalg.norm(difference, 2) * np.sqrt(1/(len(difference)))
                print(f'N-2 error: {n2_error}')

                # Calculate the convergence rate (slope)
                slope = '-----'
                if i > 0:
                    convergence_rate = (np.log10(n2_error)-np.log10(prev_error)) / (np.log10(dx)-np.log10(dxs[i-1]))
                    slope = f'{(convergence_rate):.3f}'
                print(f'Convergence rate: {slope}')
                print()

                # Calculate the relative error percentage
                sum_abs_diff = np.sum(np.abs(difference))
                sum_abs_ref = np.sum(np.abs(reference_data))
                relative_error = sum_abs_diff / sum_abs_ref
                relative_error_percentage = relative_error * 100
                print(f'Relative error percentage: {relative_error_percentage:.2f}%')

                # Write to error analysis file
                ea_file.write(f'{dt} \t| {dx} \t| {dy} \t| {(n2_error):.3e} \t| {slope} \t| {relative_error_percentage:.2f}\n')
                prev_error = n2_error

                # Store dt and n2_error for later plotting
                n2_errors.append(n2_error)

            # Convert dt_values and n2_errors to numpy arrays
            dx_values = np.array(dxs)
            n2_errors = np.array(n2_errors)

            # Fit a line to the log-log data (least squares method)
            if np.isnan(n2_errors).all():
                ea_file.write(f'Slope of the fitted line (least squares in log-log): NaN\n\n')

            else:

                # Calculate the logarithms of dt_values and n2_errors
                log_dx_values = np.log10(dx_values)
                log_n2_errors = np.log10(n2_errors)

                try:
                    coefficients = np.polyfit(log_dx_values, log_n2_errors, 1)  # 1 indicates a linear fit
                    line_slope = coefficients[0]

                    # Create the linear fit function in log-log space
                    linear_fit = np.poly1d(coefficients)

                    # Plot the data of N2-error vs dt in log-log scale
                    plt.figure()
                    plt.loglog(dx_values, n2_errors, 'o', color='blue', label='N-2 Error')
                    plt.loglog(dx_values, 10**linear_fit(log_dx_values), color='red', label='Linear Fit', linestyle='--')
                    plt.xlabel('dx')
                    plt.ylabel('N-2 Error')
                    plt.title(f'N-2 Error vs dx for {method}')
                    plt.legend()
                    plt.grid()
                    plt.savefig(f'{error_analysis_dir}/n2_error_vs_dx_loglog.png')
                    plt.close()

                    # Write the slope of the fitted line to the file
                    ea_file.write(f'Slope of the fitted line (least squares in log-log): {line_slope:.6f}\n\n')
                
                except Exception as e:
                    print(f'Error during linear fit: {e}')
                    ea_file.write(f'Error during linear fit: {e}\n\n')


if __name__ == '__main__':
    main()