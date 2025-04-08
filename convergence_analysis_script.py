from functions import *

def main():

    thetas = ['0.50', '0.66', '1.00']
    methods = ['SSI-ADI', 'theta-ADI']
    real_type = 'double'
    problem = 'MONODOMAIN' # DIFF, LINMONO
    cell_model = 'AFHN' # AFHN
    serial_or_gpu = 'SERIAL' # SERIAL, GPU

    dts = [0.025, 0.0125, 0.00625] # Works for DIFF and MONODOMAIN (with AFHN)
    # dts = [0.01, 0.005, 0.0025, 0.00125] # Works for LINMONO
    dts[::-1].sort()
    dts = [f'{dt}' for dt in dts]
    
    alphas = [0.5] # Works for DIFF and MONODOMAIN (with AFHN)
    # alphas = [0.1]

    for a in alphas:
        # 1st order approx (dt = a*dx²)
        # 2nd order approx (dt = a*dx)
        
        dxs = [f'{(float(dt) / a):.5f}' for dt in dts]
        run_script_for_convergence_analysis(a, serial_or_gpu, real_type, problem, cell_model, methods, dts, dxs, thetas)

if __name__ == "__main__":
    main()



