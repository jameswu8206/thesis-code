import subprocess
import re
import os

configs = {
    'N': [10, 15, 20, 30, 40, 50],
    'MAXITER': [3, 10,30,50],
    'FORMULATION': [1,2]
}

def modify_config(formulation, n, dt, maxiter):
    h_file = "/Users/jameswu8206/thesis/thesis-codes/mpc_example/include/mpc_pipeline_config.h"
    with open(h_file, 'r') as f:
        content = f.read()
    
    content = re.sub(r'#define FORMULATION\s+\d+', f'#define FORMULATION {formulation}', content)
    content = re.sub(r'#define N\s+\d+', f'#define N {n}', content)
    content = re.sub(r'#define DT\s+[\d\.]+', f'#define DT {dt:.4f}', content)
    content = re.sub(r'#define MAXITER\s+\d+', f'#define MAXITER {maxiter}', content)
    
    with open(h_file, 'w') as f:
        f.write(content)

def recompile():
    res = subprocess.run(
        ["make", "-B", "-j4"],
        cwd="/Users/jameswu8206/thesis/thesis-codes/mpc_example/build",
        capture_output=True,
        text=True
    )
    if res.returncode != 0:
        print("COMPILE FAILED")
        print(res.stdout)
        print(res.stderr)
        raise RuntimeError("Compilation failed")

def run_test():
    res = subprocess.run(["./mpc_test"], cwd="/Users/jameswu8206/thesis/thesis-codes/mpc_example/build", capture_output=True, text=True)
    return res.stdout

def parse_output(output):
    solver_time = "N/A"
    tsim = "N/A"
    cost = "N/A"
    status_fail = False
    fail_msg = ""
    
    if "Solver failed!" in output:
        status_fail = True
        m = re.search(r'Solver failed! status=\S+ \((.*?)\)', output)
        if m:
            fail_msg = m.group(1)
        else:
            fail_msg = "Failed"
            
    m_solver = re.search(r'Solver time\s*-> total: [\d\.]+ ms, avg: ([\d\.]+) ms', output)
    if m_solver:
        solver_time = m_solver.group(1)

    m_tsim = re.search(r'Actual Tsim\s*->\s*([\d\.]+) s', output)
    if m_tsim:
        tsim = m_tsim.group(1)
        
    m_cost = re.search(r'Final End Cost \(log10\) -> ([\-\d\.]+)', output)
    if m_cost:
        cost = m_cost.group(1)
        
    if status_fail:
        return f"{fail_msg}"
    else:
        return f"{solver_time}/{cost}"

def main():
    print("Starting grid search...")
    for formulation in configs['FORMULATION']:
        print(f"\n# Formulation {formulation}\n")
        # Header row with MAXITER values
        print("| N \\ MAXITER | " + " | ".join(map(str, configs['MAXITER'])) + " |")
        print("|---|" + "|".join(["---"] * len(configs['MAXITER'])) + "|")
        
        # Data rows for each N
        for n in configs['N']:
            row = [f"**{n}**"]
            dt = 2.0 / n
            for maxiter in configs['MAXITER']:
                modify_config(formulation, n, dt, maxiter)
                recompile()
                out = run_test()
                parsed = parse_output(out)
                
                row.append(parsed)
            print("| " + " | ".join(row) + " |")

if __name__ == "__main__":
    main()
