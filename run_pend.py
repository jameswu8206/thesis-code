import argparse
import re
import subprocess
import os

DEFAULT_CONFIGS = {
    'N': [10, 15, 20, 30, 40, 50],
    'MAXITER': [3, 10,30, 50],
    'FORMULATION': [0,1, 2]
}

# --- IMPORTANT: Update this to the exact path of your C file ---
DEFAULT_C_FILE_PATH = "/Users/jameswu8206/thesis/thesis-codes/mpc_example/invert_pendulum.c"
DEFAULT_BUILD_DIR = "/Users/jameswu8206/thesis/thesis-codes/mpc_example/build"
DEFAULT_EXECUTABLE = "./pendulum"

def modify_c_file(c_file_path, formulation, n, dt, maxiter):
    with open(c_file_path, 'r') as f:
        content = f.read()
    
    # Update macros
    content = re.sub(r'#define FORMULATION\s+\d+', f'#define FORMULATION {formulation}', content)
    content = re.sub(r'#define N\s+\d+', f'#define N {n}', content)
    content = re.sub(r'#define DT\s+[\d\.]+', f'#define DT {dt:.4f}', content)
    
    # Update OSQP max_iter setting (this will replace it in all the #if / #elif blocks)
    content = re.sub(r'settings->max_iter\s*=\s*\d+;', f'settings->max_iter = {maxiter};', content)
    
    with open(c_file_path, 'w') as f:
        f.write(content)


def recompile(build_dir):
    subprocess.run(["make", "-B", "-j4"], cwd=build_dir, capture_output=True)


def run_test(build_dir, executable):
    res = subprocess.run([executable], cwd=build_dir, capture_output=True, text=True)
    return res.stdout

def parse_output(output):
    iter_time = "N/A"
    cost = "N/A"
    status_fail = False
    fail_msg = ""
    
    # Check for solver failure
    if "Solver failed!" in output:
        status_fail = True
        m = re.search(r'Solver failed! status=\S+ \((.*?)\)', output)
        if m:
            fail_msg = m.group(1)
        else:
            fail_msg = "Failed"
            
    # We are now grabbing ONLY Iteration time and Final Cost
    m_iter = re.search(r'Iteration time-> total: [\d\.]+ ms, avg: ([\d\.]+) ms', output)
    if m_iter:
        iter_time = m_iter.group(1)
        
    m_cost = re.search(r'Final End Cost \(log10\) -> ([\-\d\.]+)', output)
    if m_cost:
        cost = m_cost.group(1)
        
    if status_fail:
        return f"{fail_msg}/N/A"
    else:
        return f"{iter_time}/{cost}"

def _parse_int_list(value):
    if value is None:
        return None
    if isinstance(value, list):
        return [int(v) for v in value]
    parts = [v.strip() for v in str(value).split(',') if v.strip()]
    return [int(v) for v in parts]


def _parse_args():
    parser = argparse.ArgumentParser(description="Grid search for OSQP pendulum.")
    parser.add_argument(
        "--n-values",
        "-n",
        help="Comma-separated list of N values (e.g. 10,20,30).",
    )
    parser.add_argument(
        "--maxiter",
        "-m",
        help="Comma-separated list of OSQP max_iter values.",
    )
    parser.add_argument(
        "--formulation",
        "-f",
        help="Comma-separated list of formulation IDs.",
    )
    parser.add_argument(
        "--c-file",
        help="Path to the pendulum C source file.",
        default=DEFAULT_C_FILE_PATH,
    )
    parser.add_argument(
        "--build-dir",
        help="Build directory containing the pendulum executable.",
        default=DEFAULT_BUILD_DIR,
    )
    parser.add_argument(
        "--exe",
        help="Executable name or path relative to build directory.",
        default=DEFAULT_EXECUTABLE,
    )
    return parser.parse_args()


def _resolve_configs(args):
    configs = {
        "N": _parse_int_list(args.n_values) or DEFAULT_CONFIGS["N"],
        "MAXITER": _parse_int_list(args.maxiter) or DEFAULT_CONFIGS["MAXITER"],
        "FORMULATION": _parse_int_list(args.formulation) or DEFAULT_CONFIGS["FORMULATION"],
    }
    return configs


def main():
    args = _parse_args()
    configs = _resolve_configs(args)
    print("Starting grid search for OSQP...")
    
    for formulation in configs['FORMULATION']:
        print(f"\n### OSQP Formulation {formulation}\n")
        print("| N \\ MAXITER | " + " | ".join(map(str, configs['MAXITER'])) + " |")
        print("|---|" + "|".join(["---"] * len(configs['MAXITER'])) + "|")
        
        for n in configs['N']:
            dt = 1.0 / n
            row = [f"**{n}**"]
            
            for maxiter in configs['MAXITER']:
                modify_c_file(args.c_file, formulation, n, dt, maxiter)
                recompile(args.build_dir)
                out = run_test(args.build_dir, args.exe)
                parsed = parse_output(out)
                row.append(parsed)
                
            # Print the row dynamically to the terminal
            print("| " + " | ".join(row) + " |")

if __name__ == "__main__":
    main()
