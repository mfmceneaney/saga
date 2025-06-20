from saga.orchestrate import create_jobs, submit_jobs

# Create job submission structure
methods = {"method":["HB","LF"]}
fitvars = {"asymfitvars":[["costheta1"],["costheta2"]]}
seeds   = {"inject_seed":[2**i for i in range(4)]}

# Results file paths and config
base_dir    = "results/"
submit_path = base_dir+"submit.sh"
yaml_path   = base_dir+"args.yaml"
out_path    = base_dir+"jobs.txt"
configs = dict(
    methods,
    **fitvars,
    **seeds
)
create_jobs(configs,base_dir,submit_path,yaml_path)
submit_jobs(configs,base_dir,submit_path,out_path,dry_run=True)
