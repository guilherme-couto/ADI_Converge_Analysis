from functions import *

def main():

    thetas = ['0.50', '0.66', '1.00']
    methods = ['SSI-ADI', 'theta-ADI']
    real_type = 'double'
    problem = 'MONOAFHN' # DIFF, MONOAFHN, LINMONO
    serial_or_gpu = 'GPU' # SERIAL, GPU
    serial_or_gpu = serial_or_gpu.upper()

    dts = [0.025, 0.0125, 0.00625] # Works for DIFF and MONOAFHN
    # dts = [0.01, 0.005, 0.0025, 0.00125] # Works for LINMONO
    dts[::-1].sort()
    dts = [f'{dt:.8f}' for dt in dts]
    
    alphas = [0.5] # Works for DIFF and MONOAFHN
    # alphas = [0.1]

    for a in alphas:
        # 1st order approx (dt = a*dx²)
        # 2nd order approx (dt = a*dx)
        
        dxs = [f'{(float(dt) / a):.6f}' for dt in dts]
        run_script(a, thetas, methods, dts, dxs, real_type, serial_or_gpu, problem)

if __name__ == "__main__":
    main()



