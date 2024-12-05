import numpy as np
import os
from functions import *

def main():
    # dts = ['0.00500', '0.01000', '0.02000', '0.04000', '0.08000', '0.10000']
    dts = ['0.00500', '0.01000', '0.02000', '0.04000'] # Dont work for MONODOMAIN  with dx=0.0005, but work for CABLEEQ
    # dts = ['0.00050', '0.00080', '0.00100', '0.00160']
    methods = ['theta-ADI'] #'SSI-ADI', 'theta-ADI', 'SSI-CN' (CABLEEQ), 'theta-RK2' (CABLEEQ)
    thetas = ['0.50', '0.66', '1.00']

    real_type = 'double'
    serial_or_gpu = 'SERIAL'
    problem = 'MONODOMAIN'
    cell_model = 'AFHN' # 'AFHN', 'TT2'
    
    # Find the largest dx to use as base for the read rate
    # base_dx = 0
    # for dx in dxs:
    #     if float(dx) > base_dx:
    #         base_dx = float(dx)
    # base_dx = 0.01
    
    # Read reference solution
    reference_dt = '0.00010'
    reference_dx = '0.00050'
    reference_solution_path = f'./reference_solutions/{real_type}/{problem}/{cell_model}/last_{reference_dt}_{reference_dx}.txt'
    
    if not os.path.exists(reference_solution_path):
        raise FileNotFoundError(f'Reference solution not found at {reference_solution_path}')
    
    base_dx = float(reference_dx)
    print(f'Reading reference solution from {reference_solution_path}')
    reference_data = read_values_with_rate(reference_solution_path, int(base_dx/float(reference_dx)))
    print(f'Reference solution read successfully. Total size: {len(reference_data)}')
    print()
    
    dx = reference_dx
    rate = int(base_dx/float(dx))
    print(f'Reading files with rate {rate}')

    for method in methods:
        # Create error analysis file
        error_analysis_dir = f'./simulation_files/analysis/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
        if not os.path.exists(error_analysis_dir):
            os.makedirs(error_analysis_dir)
            
        error_analysis_path = f'{error_analysis_dir}/error_analysis.txt'
        ea_file = open(error_analysis_path, 'w')
        
        # Prepare error analysis file
        ea_file.write(f'For method {method}\n')
        
        if 'theta' not in method:
            ea_file.write(f'dt \t\t| dx \t\t| N-2 Error \t| slope\n')
            ea_file.write('---------------------------------------------------------\n')
            
            # Comparative plot
            plt.figure()
            plt.plot(reference_data, label='ref')

            # Initialize lists to store dt and n2_error values
            n2_errors = []

            # Iterate over the dts and dxs
            prev_error = 0
            for i in range(len(dts)):
                dt = dts[i]
                                    
                simulation_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe/last_{dt}_{dx}.txt'
                if not os.path.exists(simulation_path):
                    raise FileNotFoundError(f'Simulation file not found at {simulation_path}')
                
                # Read simulation file
                print(f'Reading simulation file from {simulation_path}')
                simulation_data = read_values_with_rate(simulation_path, rate)
                print(f'Simulation file read successfully. Total size: {len(simulation_data)}')
                
                # Plot data to compare
                plt.plot(simulation_data, label=f'dt={dt}')
                
                # Plot difference map between reference and simulation
                difference = np.array(simulation_data) - np.array(reference_data)
                if problem != 'CABLEEQ':
                    plot_difference_map_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx)
                else:
                    plot_difference_vector_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx)
                
                # Calculate the N-2 error
                n2_error = np.linalg.norm(difference, 2) * np.sqrt(1/(len(difference)))
                print(f'N-2 error: {n2_error}')
                
                # Calculate the convergence rate (slope)
                slope = '-----'
                if i > 0:
                    convergence_rate = (np.log10(n2_error)-np.log10(prev_error)) / (np.log10(float(dts[i]))-np.log10(float(dts[i-1])))
                    slope = f'{(convergence_rate):.3f}'
                print(f'Convergence rate: {slope}')
                print()
                
                # Write to error analysis file
                ea_file.write(f'{dt}\t| {dx}\t| {(n2_error):.6f} \t| {slope}\n')
                
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
            dt_values = np.array(dts, dtype=float)
            n2_errors = np.array(n2_errors)

            # Calculate the logarithms of dt_values and n2_errors
            log_dt_values = np.log10(dt_values)
            log_n2_errors = np.log10(n2_errors)

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
            plt.title(f'N-2 Error vs dt for {method}')
            plt.legend()
            plt.grid()
            plt.savefig(f'{error_analysis_dir}/n2_error_vs_dt_loglog.png')
            plt.close()

            # Write the slope of the fitted line to the file
            ea_file.write(f'Slope of the fitted line (least squares in log-log): {line_slope:.6f}\n\n')
        
        else:
            thetas_linear_fits = []
            for theta in thetas:
                ea_file.write(f'with theta={theta}\n')
                ea_file.write(f'dt \t\t| dx \t\t| N-2 Error \t| slope\n')
                ea_file.write('---------------------------------------------------------\n')
                
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
                                        
                    simulation_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/{theta}/lastframe/last_{dt}_{dx}.txt'
                    if not os.path.exists(simulation_path):
                        raise FileNotFoundError(f'Simulation file not found at {simulation_path}')
                    
                    # Read simulation file
                    print(f'Reading simulation file from {simulation_path}')
                    simulation_data = read_values_with_rate(simulation_path, rate)
                    print(f'Simulation file read successfully. Total size: {len(simulation_data)}')                        

                    # Plot difference map between reference and simulation
                    difference = np.array(simulation_data) - np.array(reference_data)
                    if problem != 'CABLEEQ':
                        plot_difference_map_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta)
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
                        convergence_rate = (np.log10(n2_error)-np.log10(prev_error)) / (np.log10(float(dts[i]))-np.log10(float(dts[i-1])))
                        slope = f'{(convergence_rate):.3f}'
                    print(f'Convergence rate: {slope}')
                    print()
                    
                    # Write to error analysis file
                    ea_file.write(f'{dt}\t| {dx}\t| {(n2_error):.6f} \t| {slope}\n')
                    
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
                dt_values = np.array(dts, dtype=float)
                n2_errors = np.array(n2_errors)

                # Calculate the logarithms of dt_values and n2_errors
                log_dt_values = np.log10(dt_values)
                log_n2_errors = np.log10(n2_errors)

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
            
            plt.figure()
            for t in range(len(thetas)):
                plt.loglog(dt_values, thetas_linear_fits[t], label=f'theta={thetas[t]}', linestyle='-')
            plt.xlabel('dt')
            plt.ylabel('N-2 Error')
            plt.title(f'N-2 Error vs dt for {method}')
            plt.legend()
            plt.grid()
            plt.savefig(f'{error_analysis_dir}/n2_error_vs_dt_loglog.png')
            plt.close()

if __name__ == '__main__':
    main()