from bids import BIDSLayout, BIDSLayoutIndexer
import bids
import os
import json
import re
import yaml
import shutil

bids.config.set_option('extension_initial_dot', True)

def bids(root=None, datatype=None, prefix=None, suffix=None, subject=None, session=None,include_subject_dir=True,include_session_dir=True,**entities):
    """Helper function for generating bids paths for snakemake workflows

    File path is of the form:

    [root]/[sub-{subject}]/[ses-{session]/[prefix]_[sub-{subject}]_[ses-{session}]_[{key}-{val}_ ... ]_[suffix]

    root -- root folder to include in the path (e.g. 'results'))
    datatype -- folder to include after sub-/ses- (e.g. anat, dwi )
    prefix -- string to prepend to the file name (typically not defined, unless you want tpl-{tpl}, or a datatype)
    suffix -- bids suffix including extension (e.g. 'T1w.nii.gz')
    subject -- subject to use, for folder and filename
    session -- session to use, for folder and filename
    include_subject_dir -- whether to include the sub-{subject} folder if subject defined (default: True)
    include_session_dir -- whether to include the ses-{session} folder if session defined (default: True)
    **entities -- dictionary of bids entities (e.g. space=T1w for space-T1w)

    Returns: bids-like file path

    Example:

        Below is a rule using bids naming for input and output:

        rule proc_img:
            input: 'sub-{subject}_T1w.nii.gz'
            output: 'sub-{subject}_space-snsx32_desc-preproc_T1w.nii.gz'

        With bids() you can instead use:

         rule proc_img:
            input: bids(subject='{subject}',suffix='T1w.nii.gz')
            output: bids(subject='{subject}',space='snsx32',desc='preproc',suffix='T1w.nii.gz')

        Note that here we are not actually using "functions as inputs" in snakemake, which would require
        a function definition with wildcards as the argument, and restrict to input/params, but bids()
        is being used simply to return a string.

        Also note that space, desc and suffix are NOT wildcards here, only {subject} is.
        This makes it easy to combine wildcards and non-wildcards with bids-like naming.

        However, you can still use bids() in a lambda function, this is especially useful if your wildcards
        are named the same as bids entities (e.g. {subject}, {session}, {task} etc..):
 
        rule proc_img:
            input: lambda wildcards: bids(**wildcards,suffix='T1w.nii.gz')
            output: bids(subject='{subject}',space='snsx32',desc='preproc',suffix='T1w.nii.gz')

        Or another example where you may have many bids-like wildcards used in your workflow:
        
        rule denoise_func:
            input: lambda wildcards: bids(**wildcards, suffix='bold.nii.gz')
            output: bids(subject='{subject}',session='{session}',task='{task}',acq='{acq}',desc='denoise',suffix='bold.nii.gz')

        In this example, all the wildcards will be determined from the output and passed on to bids() for inputs.
        The output filename will have a 'desc-denoise' flag added to it.


        Also note that even if you supply entities in a different order, the
        entities will be ordered based on the OrderedDict defined here.
        If you entities not known are provided, they will be just be placed
        at the end (before the suffix), the the order provided.

        Note: For maximum flexibility all arguments are optional (if none are specified, will return empty string)

    -- some code adapted from mne-bids
      https://mne.tools/mne-bids/stable/_modules/mne_bids/utils.html


    """

    from collections import OrderedDict
    from os.path import join
    

    #replace underscores in keys (needed to that users can use reserved keywords by appending a _)
    entities = { k.replace('_', ''): v for k, v in entities.items() }
       
 
   
    #strict ordering of bids entities is specified here:
    order = OrderedDict([('task', None),
                         ('acq', None),
                         ('ce', None),
                         ('rec', None),
                         ('dir', None),
                         ('run', None),
                         ('mod', None),
                         ('echo', None),
                         ('hemi', None),
                         ('space', None),
                         ('res', None),
                         ('den', None),
                         ('label', None),
                         ('desc', None)])

    #now add in entities (this preserves ordering above)
    for key, val in entities.items():
        order[key] = val

    #initialize lists for filename and folder
    # will append to these, then '_'.join() os.path.join() respectively
    filename = []
    folder = []

    #root directory
    if isinstance(root,str):
        folder.append(root)

    #if prefix is defined, put it before other anything else
    if isinstance(prefix, str):
        filename.append(prefix)

    #if subject defined then append to file and folder
    if isinstance(subject,str):
        if include_subject_dir is True:
            folder.append(f'sub-{subject}')
        filename.append(f'sub-{subject}')

    #if session defined then append to file and folder
    if isinstance(session,str):
        if include_session_dir is True:
            folder.append(f'ses-{session}')
        filename.append(f'ses-{session}')
    
    if isinstance(datatype,str):
        folder.append(datatype)
    
    #add the entities
    for key, val in order.items():
        if val is not None:
            filename.append(f'{key}-{val}')

    #if suffix is defined, append it
    if isinstance(suffix, str):
        filename.append(suffix)


    if len(filename) == 0:    
        return ''

    #now, join up the lists:
    filename = '_'.join(filename)

    if len(folder)>0:
        filename = join(*folder,filename)
    
    return filename


def get_filtered_ziplist_index(zip_list, wildcards, subj_wildcards):
    """ Use this function when you have wildcards for a single scan instance, 
        and want to know the index of that scan, amongst that subject's scan
        instances. 

    Example:
    >>> snakebids.get_filtered_ziplist_index({'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'] ,'acq': ['98','98','98','98','99','99','99','99'], 'subject': ['01','01','02','02','01','01','02','02' ] }, {'dir': 'PA', 'acq': '99', 'subject': '01'}, { 'subject': '{subject}' })
    3

    """
    #get the subject/(session) dict:
    subj_dict = { key:wildcards[key]  for key in subj_wildcards.keys()}

    #now filter the list based on subj_wildcards
    zip_list_filtered = filter_list(zip_list,subj_dict)

    #get the index of the wildcard from this filtered list
    indices = filter_list(zip_list_filtered,wildcards, return_indices_only=True)
    if len(indices) == 1:
        return indices[0]
    else:
        return indices

      


# this function is used when you are expanding over some subset of the wildcards
#  i.e. if your output file doesn't contain all the wildcards in input_wildcards
def filter_list(zip_list, wildcards, return_indices_only=False):
    size_list = len(next(iter(zip_list)))
    keep_indices = set()
    for key,val in wildcards.items():
        #get indices where wildcard exists
        if not key in zip_list.keys():
            continue
        indices = set([i for i,v in enumerate(zip_list[key]) if v == val])
        if len(keep_indices) == 0:
            keep_indices = indices
        else:
            keep_indices = keep_indices.intersection(indices)
    #now we have the indices, so filter the lists
    if return_indices_only:
        return list(keep_indices)
    else:
        return {key: [ zip_list[key][i] for i in keep_indices  ] for key,val in zip_list.items()}



# use pybids to get the paths to input images, then create a input path and wildcards for each suffix type

def read_bids_tags(bids_json=None):
    if bids_json == None:
        bids_json = os.path.join(os.path.dirname(os.path.realpath(__file__)),'bids_tags.json')
    with open(bids_json, 'r') as infile:
        bids_tags = json.load(infile)
    return bids_tags

def get_input_config_from_bids(config, bids_layout, inputs_dict, limit_to=None, **filters ):
    """ returns: dict with input_path and input_wildcards"""

    bids_tags = read_bids_tags()

    config.update( dict({'input_path': {}, 'input_zip_lists': {}, 'input_lists': {}, 'input_wildcards': {}}))


    if limit_to == None:
        inputs_to_iterate = inputs_dict.keys()
    else:
        inputs_to_iterate = limit_to

    

    for input_name in inputs_to_iterate:

        if config['debug']==True: print(f'grabbing pybids inputs for {input_name}..')

        imgs, = [bids_layout.get(**inputs_dict[input_name]['filters'], **filters)]
        if len(imgs) == 0:
            print(f'WARNING: no images found for {input_name}')
            continue
        
        if config['debug']==True: print(f'  found {len(imgs)} images')
        paths = set()
        zip_lists = {}
        input_lists = {}
        wildcards = {}
        for img in imgs:
            path = img.path
            for wildcard_name in inputs_dict[input_name]['wildcards']:

                if wildcard_name in bids_tags:
                    tag = bids_tags[wildcard_name]
                else:
                    tag = wildcard_name  #if it's not in the bids_tags dictionary, then just use the name itself as the tag

                   
                #this changes e.g. sub-001 to sub-{subject} in the path (so snakemake can use the wildcards)
                if wildcard_name in img.get_entities():
                    
                    if config['debug']==True: print(f'    wildcard {wildcard_name} found entities for {path}')
                    ##HACK FIX FOR acq vs acquisition etc  -- should eventually update the bids() function to also use bids_tags.json, where e.g. acquisition -> acq is defined.. -- then, can use wildcard_name instead of out_name.. 
                    if wildcard_name not in ['subject', 'session']:
                        out_name = tag
                    else:
                        out_name = wildcard_name
 
                    if out_name not in zip_lists:
                        zip_lists[out_name] = []
                        input_lists[out_name] = set()
                        wildcards[out_name] = {}

                    if config['debug']==True: print(f'    wildcard {wildcard_name} found entities for {path}')
                    pattern = '{tag}-([a-zA-Z0-9]+)'.format(tag=tag)
                    replace = '{tag}-{{{replace}}}'.format(tag=tag,replace=out_name)
                    match = re.search(pattern,path)
                    replaced = re.sub(pattern,  replace , path)
                    #update the path with the {wildcards} -- uses the value from the string (not from the pybids entities), since that has issues with integer formatting (e.g. for run=01)
                    path = replaced
                    zip_lists[out_name].append(match[1])
                    input_lists[out_name].add(match[1])
                    wildcards[out_name] = f'{{{out_name}}}'
                
            paths.add(path)
        

        #now, check to see if unique
        if len(paths) > 1:
            print(f'WARNING: more than one snakemake filename for {input_name}, taking the first')
            print(f'  To correct this, use the --filter_{input_name} option to narrow the search')
            print(paths)
    
        in_path = list(paths)[0]

        #convert sets to lists
        for key,val in input_lists.items():
            input_lists[key] = list(val)
            
            
        config['input_path'][input_name] = in_path
        config['input_zip_lists'][input_name] = zip_lists
        config['input_lists'][input_name] = input_lists
        config['input_wildcards'][input_name] = wildcards


    return config



def generate_inputs_config(config,limit_to=None):
    """ returns: updated config dict; function will also write the inputs_config.yml to standard output path """
    


    if not 'search_terms' in config.keys():
        config['search_terms'] = dict()

    if 'participant_label'  and 'exclude_participant_label' in config.keys():
        if not config['participant_label'] == None and not config['exclude_participant_label'] == None:
            print('ERROR: cannot define both participant_label and exclude_participant_label at the same time')
            return None

    #add participant_label or exclude_participant_label to search terms (if defined)
    # we make the subject key in search_terms a list so we can have both participant_label and exclude_participant_label defined 
    if 'participant_label' in config.keys():
        if not config['participant_label'] == None:
            if not 'subject' in config['search_terms'].keys():
                config['search_terms']['subject'] = []
            if isinstance(config['participant_label'], list): 
                config['search_terms']['subject'] = config['search_terms']['subject'] + config['participant_label']
            else:
                config['search_terms']['subject'].append(config['participant_label'])

    if 'exclude_participant_label' in config.keys():
        if not config['exclude_participant_label'] == None:
            if not 'subject' in config['search_terms'].keys():
                config['search_terms']['subject'] = []
            if isinstance(config['exclude_participant_label'], list): # if multiple subjects to exclude, combine with with subj1|subj2|...
                exclude_string = '|'.join(config['exclude_participant_label']) 
                
            else:
                exclude_string = config['exclude_participant_label'] #if not, then string is the label itself
            config['search_terms']['regex_search'] = True
            config['search_terms']['subject'].append(f'^((?!({exclude_string})).)*$') #regex to exclude subjects
       


    #generate inputs based on config
    layout = BIDSLayout(
        config['bids_dir'],
        derivatives=config['derivatives'],
        validate=False,
        indexer=BIDSLayoutIndexer(validate=False,
                index_metadata=False)
    )

    #this will populate input_path, input_lists, input_zip_lists, and input_wildcards 
    inputs_config_dict = get_input_config_from_bids(config=config,
                                                bids_layout=layout,
                                                inputs_dict=config['pybids_inputs'],
                                                limit_to = limit_to,
                                                **config['search_terms'])

    #populate subjects, sessions and subj_wildcards in the config
    inputs_config_dict['subjects'] = layout.get_subjects(**config['search_terms'])
    inputs_config_dict['sessions'] = layout.get_sessions(**config['search_terms'])
    if len(inputs_config_dict['sessions'])  == 0:
        inputs_config_dict['subj_wildcards'] = { 'subject': '{subject}'}
    else:
        inputs_config_dict['subj_wildcards'] = { 'subject': '{subject}', 'session': '{session}' }

    #set snakemake_dir to '.' if not defined
    if not 'snakemake_dir' in config.keys():
        config['snakemake_dir'] = '.'


    #write updated config file
    inputs_config = os.path.join('config','inputs_config.yml')
    os.makedirs(os.path.dirname(inputs_config),exist_ok=True)

    with open(inputs_config, 'w') as outfile:
        yaml.dump(inputs_config_dict, outfile, default_flow_style=False)

    #add to config dict before returning
    config.update(inputs_config_dict)


    #copy pipeline_description.json to results/dataset_description.json
    pipeline_description = os.path.join(config['snakemake_dir'],'pipeline_description.json')
    dataset_description = os.path.join('results','dataset_description.json')
    if os.path.exists(pipeline_description):
        try:
            os.mkdir(os.path.dirname(dataset_description))
        except FileExistsError:
            pass
        shutil.copyfile(pipeline_description,dataset_description)

    return config

def get_wildcard_constraints(config):

    """ returns: dict containing the wildcard constraints for all the wildcards defined, using [a-zA-Z0-9]+ for all wildcards """
    bids_constraints = '[a-zA-Z0-9]+'
    return  { entity: bids_constraints for imgtype in config['input_lists'].keys() for entity in config['input_lists'][imgtype].keys() }


def write_derivative_json(snakemake, **kwargs):
    with open(snakemake.input.json,'r') as f:
      sidecar = json.load(f)

    sidecar.update({'Sources': [snakemake.input], 'Parameters': {key:val for (key,val) in snakemake.params.items()}, **kwargs })

    with open(snakemake.output.json, 'w') as outfile:
        json.dump(sidecar, outfile,indent=4)


