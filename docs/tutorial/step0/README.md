In this example we will make a workflow to smooth `bold` scans from a bids dataset.

We will start by creating a simple rule, then make this more generalizable in each step. This is the base command we are running:

In this rule, we start by hard-coding the paths for input and output, the conversion from FWHM to sigma
