import argparse
import functions

# Dict to map the option to the correct function
options = {
    "gif": functions.create_gif,
    "lastframe": functions.plot_last_frame,
    "exact": functions.plot_exact,
    "errors": functions.plot_errors,
}

def main():
    parser = argparse.ArgumentParser(description="Choose an option to execute a specific function with parameters.")

    # Main argument to choose the function
    parser.add_argument("--option", "--o", choices=options.keys(), help="Option to choose the function to execute.")

    # Arguments for the functions
    parser.add_argument("--serial_or_gpu", "--e", type=str, help="Execution mode (SERIAL or GPU)")
    parser.add_argument("--real_type", "--rt", type=str, help="The real type (e.g., float or double)")
    parser.add_argument("--problem", "--p", type=str, help="Problem identifier")
    parser.add_argument("--cell_model", "--cm", type=str, help="Cell model identifier")
    parser.add_argument("--method", "--m", type=str, help="The method to use")
    parser.add_argument("--dt", type=str, help="Time step size (must be > 0 and < 1) with 5 decimal places")
    parser.add_argument("--dx", type=str, help="Space step size (must be > 0 and < 1) with 5 decimal places")
    parser.add_argument("--theta", type=str, help="Theta value for the method")
    
    args = parser.parse_args()

    # Check if the parameters dt and dx are valid
    if not (0 < float(args.dt) < 1):
        raise ValueError("The parameter 'dt' must be a float greater than 0 and less than 1.")
    if not (0 < float(args.dx) < 1):
        raise ValueError("The parameter 'dx' must be a float greater than 0 and less than 1.")
    if not (len(args.dt.split('.')[1]) == 5):
        raise ValueError("The parameter 'dt' must have 5 decimal places.")
    if not (len(args.dx.split('.')[1]) == 5):
        raise ValueError("The parameter 'dx' must have 5 decimal places.")

    # Get the function to call
    func_to_call = options[args.option]

    # Call the function with the parameters
    func_to_call(args.serial_or_gpu, args.real_type, args.problem, args.cell_model, args.method, args.dt, args.dx, args.theta)

if __name__ == "__main__":
    main()
