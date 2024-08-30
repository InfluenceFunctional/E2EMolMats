from mxtaltools.common.config_processing import load_yaml
import yaml
from copy import copy

base_config = load_yaml('latents_base.yaml')

melts_dict = {
    'acridine/Form2': 394.68468468468467,
    'acridine/Form3': 394.68468468468467,
    'acridine/Form4': 351.3613613613614,
    'acridine/Form6': 358.8088088088088,
    'acridine/Form7': 372.5025025025025,
    'acridine/Form8': 371.7017017017017,
    'acridine/Form9': 376.02602602602605
}
counter = 0
for k, v in melts_dict.items():
    new_config = copy(base_config)
    new_config['structure_identifier'] = k
    new_config['temperature'] = v
    new_config['run_name'] = f'acridine_latents{counter}'

    with open(str(counter) + '.yaml', 'w') as outfile:
        yaml.dump(new_config, outfile, default_flow_style=False)

    counter += 1
