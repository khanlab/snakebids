from .. import generate_inputs_config

def test_t1w():
    #create config
    bids_dir = 'snakebids/tests/data/bids_t1w'
    derivatives = False
    debug = False
    pybids_inputs = { 't1': {'filters': {'suffix': 'T1w'}, 'wildcards': ['acquisition','subject','session','run']}} 
    config = { 'bids_dir': bids_dir, 'derivatives': derivatives, 'debug': debug, 'pybids_inputs': pybids_inputs }
    generate_inputs_config(config)
    assert config['input_lists'] ==  {'t1': {'acq': ['mprage'], 'subject': ['001']}}
    assert config['input_zip_lists'] == {'t1': {'acq': ['mprage'], 'subject': ['001']}}
    assert config['input_wildcards'] == {'t1': {'acq': '{acq}', 'subject': '{subject}'}}
    assert config['subjects'] == ['001']
    assert config['sessions'] == []
    assert config['subj_wildcards'] == {'subject': '{subject}'}

