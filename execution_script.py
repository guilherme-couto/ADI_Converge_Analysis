from functions import *

def main():
    dts = ['0.0050', '0.0100', '0.0200']
    dxs = ['0.0050', '0.0100', '0.0200']
    methods = ['SSI-ADI', 'theta-ADI']
    thetas = ['0.50', '0.66', '1.00']

    for method in methods:
        for i in range(len(dts)):
            dt = dts[i]
            dx = dxs[i]
            
            if method == 'SSI-ADI':
                tts = ['0.00']
            else:
                tts = thetas
            for theta in tts:
                os.system(f'./convergence {method} {dt} {dx} {theta}')
                plot_last_frame(method, dt, dx, 'float', 'GPU', 'MONOAFHN', theta)
                create_gif(method, dt, dx, 'float', 'GPU', 'MONOAFHN', theta)

if __name__ == '__main__':
    main()