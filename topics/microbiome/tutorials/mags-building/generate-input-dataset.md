# Generate input datasets for the training

Data from [PRJNA977416](https://www.ebi.ac.uk/ena/browser/view/PRJNA977416) that has been run through the workflow: https://usegalaxy.eu/u/paulzierep/h/prjna977416-bee-use-case 

## Selection of samples to keep

- Identify the number of bins per sample

    Sample | Number of bins
    --- | ---
    SRR24759596 | 0
    SRR24759597 | 2
    SRR24759598 | 4
    SRR24759599 | 4
    SRR24759600 | 8
    SRR24759601 | 3
    SRR24759602 | 3
    SRR24759603 | 4
    SRR24759604 | 4
    SRR24759605 | 2
    SRR24759606 | 4
    SRR24759607 | 3
    SRR24759608 | 1
    SRR24759609 | 4
    SRR24759610 | 2
    SRR24759611 | 1
    SRR24759612 | 3
    SRR24759613 | 4
    SRR24759614 | 2
    SRR24759615 | 2
    SRR24759616 | 4

- Keep samples with at least 4 bins

## Prepare data

History:

- Get samples from PRJNA977416
- Filter to keep the 9 samples with at least 4 bins
- Import and run the [preprocessing-for-mags workflow](https://usegalaxy.eu/u/paulzierep/w/preprocessing-for-mags)
- Import and run the [mags-host-removal](https://usegalaxy.eu/u/paulzierep/w/host-contamination-removal)

    
