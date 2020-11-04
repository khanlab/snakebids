# snakebids
Snakemake + BIDS




If we want to do some analysis that takes data from all subjects with a single task, then
our output file only has the {task} wildcard, so for the input, you need to expand over all the other wildcards, except task

You can use the `filter_list()` function to take an `input_zip_list` and automatically remove the wildcards already defined by the output wildcards


```
rule group_analysis:

    input: lambda wildcards: expand(bids(root='results',desc='denoised',suffix='bold.nii.gz',**config['input_wildcards']['bold']),\
                zip, **filter_list(config['input_zip_lists']['bold'], wildcards)) 
    output: bids(root='results',desc='group',suffix='combined.nii.gz',task='{task}')
```

another example, averaging scans within a session:
```
rule aggregate_over_session:
    input: lambda wildcards: expand(bids(root='results',desc='denoised',suffix='bold.nii.gz',**config['input_wildcards']['bold']),\
                zip, **filter_list(config['input_zip_lists']['bold'], wildcards)) 
    output: bids(root='results',desc='sesavg',suffix='bold.nii.gz',**config['subj_wildcards'])
```



# if you want a rule to run for each subj/(session), then you can use the config['subj_wildcards'] 
rule aggregate_over_session:
    input: lambda wildcards: expand(bids(root='results',desc='denoised',suffix='bold.nii.gz',**config['input_wildcards']['bold']),\
                zip, **filter_list(config['input_zip_lists']['bold'], wildcards)) 
    output: bids(root='results',desc='sesavg',suffix='bold.nii.gz',**config['subj_wildcards'])


