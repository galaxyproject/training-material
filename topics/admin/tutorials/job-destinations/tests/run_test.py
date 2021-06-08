import os

import bioblend.galaxy

from galaxy.tool_util.verify.wait import wait_on


N_CHARS = int(os.environ.get('N_CHARS', 15))
N_THREADS = os.environ.get('N_THREADS', 1)
GALAXY_URL = os.environ.get('GALAXY_URL', 'https://your-galaxy')
GALAXY_API_KEY = os.environ.get('GALAXY_API_KEY', 'adminkey')

gi = bioblend.galaxy.GalaxyInstance(GALAXY_URL, GALAXY_API_KEY)
history_id = gi.histories.create_history('test_history')['id']
job = gi.tools.paste_content('A' * N_CHARS, history_id=history_id)
hda_id = job['outputs'][0]['id']
rval = gi.tools.run_tool(history_id, 'testing', {'input1': {'src': 'hda', 'id': hda_id}})
job_id = rval['jobs'][0]['id']
wait_on(lambda: gi.jobs.show_job(job_id)['state'] == 'ok', 'waiting for test job', timeout=60)
output_id = rval['outputs'][0]['id']
gi.datasets.download_dataset(output_id)
assert gi.datasets.download_dataset(output_id).decode() == f"Running with '{N_THREADS}' threads\n"
