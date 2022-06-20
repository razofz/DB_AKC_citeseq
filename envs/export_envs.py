import os
import json


ENV_BASENAME = "DB_AKC"

os.system("conda env list --json > all-envs.json")
envs = json.load(open("all-envs.json", "r"))

relevant_envs = [os.path.basename(x) for x in envs["envs"] if ENV_BASENAME in x]

for env in relevant_envs:
    print(f"Exporting environment {env} to file {env}.yaml..")
    os.system(f"conda env export --no-builds --name {env} -f {env}.yaml")
