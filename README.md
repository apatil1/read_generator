# Read Simulator

Creates synthetic paired end DNA from given fraction of Human, Bacteria, Phix174 and Virus/Phage.
Note: List of Human, Bacteria, Phix174 and Virus/Phage fractions should be space separated.
Example: -hu 0.5 0.2 -b 0.3 0.01 -x 0.1 0.01.

The Virus/Phage fraction is taken 1 - Human - Bacteria - Phix174.

Usage:

sythetic.py [-h] -s S -p P [-hu HU [HU ...]] [-b B [B ...]] [-x X [X ...]]

-h, --help       show this help message and exit
-s S             Number of total reads
-p P             Path to directory for output files
-hu HU [HU ...]  Human DNA percentage. Default: [0.5, 0.1, 0.01, 0.001]
-b B [B ...]     Bacterial DNA percentage. Default: [0.4, 0.25, 0.1, 0.5]
-x X [X ...]     Phix174 DNA percentage. Default: [0.01, 0.001]

Example:

```bash
python sythetic.py -s 10000 -p /ouputdir
```
