import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

def read_error_file(file_path):
    """Read error analysis file with more flexible parsing"""
    dts, errors, slopes = [], [], []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        # Skip header lines (looking for the line with dashes)
        start_idx = 0
        for i, line in enumerate(lines):
            if '----' in line:
                start_idx = i + 1
                break

        # Process data lines
        for line in lines[start_idx:]:
            if not line.strip():
                continue
            # Split by whitespace and filter out empty strings
            parts = [p for p in line.replace('\t', ' ').replace('|', '').split(' ') if p]

            try:
                if len(parts) >= 4:  # Need at least dt and error
                    dt = float(parts[0])
                    error = float(parts[3])
                    slope = float(parts[4]) if parts[4].isnumeric() else np.nan

                    if not np.isnan(error):
                        
                        dts.append(dt)
                        errors.append(error)
                        slopes.append(slope)
                    
            except (ValueError, IndexError):
                if parts[0].startswith('Slope'):
                    fit_slope = float(parts[-1])
                continue
    return np.array(dts), np.array(errors), np.array(slopes), fit_slope

# Method directories and their display names
methods = {
    'FE': 'FE',
    'OS-ADI': 'OS-ADI',
    'SSI-ADI': 'SSI-ADI'
}

# Colors and markers for each method
style = {
    'FE': {'color': 'red'},
    'OS-ADI': {'color': 'orange'},
    'SSI-ADI': {'color': 'royalblue'}
}

plt.figure(figsize=(12, 8))

for method in methods:
    file_path = os.path.join(method, 'error_analysis.txt')
    
    if not os.path.exists(file_path):
        print(f"Warning: File not found for method {method} - {file_path}")
        continue
    
    try:
        dt, error, slopes, fit_slope = read_error_file(file_path)
        
        if len(dt) > 0:
            # Plot the convergence points
            # plt.loglog(dt, error, 'o', **style[method], label=f'{methods[method]}')
            
            # Plot linear fit
            log_dt = np.log10(dt)
            log_error = np.log10(error)
            fit = np.polyfit(log_dt, log_error, 1)
            plt.loglog(dt, 10**np.polyval(fit, log_dt), 'o',
                      label=f'{methods[method]} (slope: {fit[0]:.3f})',
                      **style[method], linestyle='-', linewidth=2, markersize=8)
            
            # Add slope annotation for each point
            # for i in range(0, len(dt)):
            #     if not np.isnan(slopes[i]):
            #         plt.annotate(f'{slopes[i]:.2f}', 
            #                     (dt[i], error[i]),
            #                     textcoords="offset points", 
            #                     xytext=(10,-5 if i%2 else 5),
            #                     ha='center', fontsize=9)
    except Exception as e:
        print(f"Error processing {method}: {str(e)}")

# Plot configuration
plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))


plt.title('Convergence Analysis', fontsize=20, fontweight='bold')
plt.xlabel('Time step (ms) ($\\times$ 1e-3)', fontsize=20)
plt.ylabel('L2 Error Norm', fontsize=20)
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend(fontsize=20)
dts_values = np.array([1, 2, 3, 4, 5])  # Valores que serão mostrados nos ticks
dts_ticks = dts_values * 1e-3  # Valores reais usados no gráfico
plt.xticks(dts_ticks, dts_values)
plt.tick_params(axis='both', which='major', labelsize=18)

plt.tight_layout()
plt.savefig('convergence_analysis.png')
plt.close()

print("Plot saved as 'convergence_analysis.png'")