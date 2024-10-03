import numpy as np
import os
from functions import *

def main():
    dts = ['0.00500', '0.01000', '0.02000', '0.04000', '0.08000', '0.16000']
    dxs = ['0.00500', '0.01000', '0.02000']
    methods = ['SSI-ADI', 'theta-ADI'] #'SSI-ADI', 'theta-ADI', 'SSI-CN' (CABLEEQ)
    thetas = ['0.50', '0.66', '1.00']

    # dts = ['0.00050']
    methods = ['SSI-CN']

    real_type = 'double'
    serial_or_gpu = 'SERIAL'
    problem = 'CABLEEQ'
    cell_model = 'TT2' # 'AFHN', 'TT2'
    
    # Find the largest dx to use as base for the read rate
    base_dx = 0
    for dx in dxs:
        if float(dx) > base_dx:
            base_dx = float(dx)
    base_dx = 0.01
    
    # Read reference solution
    reference_dt = '0.00050'
    reference_dx = '0.01000'
    reference_solution_path = f'./reference_solutions/{real_type}/{problem}/{cell_model}/last_{reference_dt}_{reference_dx}.txt'
    
    if not os.path.exists(reference_solution_path):
        raise FileNotFoundError(f'Reference solution not found at {reference_solution_path}')
    
    print(f'Reading reference solution from {reference_solution_path}')
    reference_data = read_values_with_rate(reference_solution_path, int(base_dx/float(reference_dx)))
    print(f'Reference solution read successfully. Total size: {len(reference_data)}')
    print(reference_data)
    print()
    
    dx = '0.01000'
    rate = int(base_dx/float(dx))
    print(f'Reading files with rate {rate}')
    plt.figure()
    plt.plot(reference_data, label='ref')
    for method in methods:
        # Create error analysis file
        error_analysis_dir = f'./simulation_files/analysis/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}'
        if not os.path.exists(error_analysis_dir):
            os.makedirs(error_analysis_dir)
            
        error_analysis_path = f'{error_analysis_dir}/error_analysis.txt'
        ea_file = open(error_analysis_path, 'w')
        
        # Prepare error analysis file
        ea_file.write(f'For method {method}\n')
        
        if method != 'theta-ADI':
            ea_file.write(f'dt \t\t| dx \t\t| N-2 Error \t| slope\n')
            ea_file.write('---------------------------------------------------------\n')
                
            # Iterate over the dts and dxs
            prev_error = 0
            for i in range(len(dts)):
                dt = dts[i]
                # dx = dxs[i]
                dx = '0.01000'
                                    
                simulation_path = f'./simulation_files/outputs/{serial_or_gpu}/{real_type}/{problem}/{cell_model}/{method}/lastframe/last_{dt}_{dx}.txt'
                if not os.path.exists(simulation_path):
                    raise FileNotFoundError(f'Simulation file not found at {simulation_path}')
                
                # Read simulation file
                print(f'Reading simulation file from {simulation_path}')
                simulation_data = read_values_with_rate(simulation_path, rate)
                print(f'Simulation file read successfully. Total size: {len(simulation_data)}')
                plt.plot(simulation_data, label=f'sim{dt}')
                
                # Plot difference map between reference and simulation
                difference = np.array(simulation_data) - np.array(reference_data)
                if problem != 'CABLEEQ':
                    plot_difference_map_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx)
                else:
                    plot_difference_vector_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx)
                
                # Calculate the N-2 error
                # n2_error = np.linalg.norm(difference, 2) * float(dx)
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
                
            ea_file.write('\n\n')
        
        else:
            for theta in thetas:
                ea_file.write(f'with theta={theta}\n')
                ea_file.write(f'dt \t\t| dx \t\t| N-2 Error \t| slope\n')
                ea_file.write('---------------------------------------------------------\n')
                
                # Iterate over the dts and dxs
                prev_error = 0
                for i in range(len(dts)):
                    dt = dts[i]
                    # dx = dxs[i]
                    dx = '0.01000'
                                        
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
                        plot_difference_vector_from_data(difference, serial_or_gpu, real_type, problem, cell_model, method, dt, dx, theta)               
                    
                    # Calculate the N-2 error
                    n2_error = np.linalg.norm(difference, 2) * float(dx)
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
                    
                ea_file.write('\n\n')
    
    plt.grid()
    plt.legend()
    plt.savefig('comp_cable.png')
    plt.close()

if __name__ == '__main__':
    main()