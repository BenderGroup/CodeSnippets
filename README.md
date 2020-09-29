# CodeSnippets
Small snippets of code which don't require their own repository

## SCINS

#### Definition of SCINS [1]
SCINS describes a reduced graph of the scaffold of a molecule. It is characterized by a string of numbers in the format ABCDE-FGHI-JKLM.

| Character | Description | Details |
| :---: | :--- | :--- |
| A | Number of Chain Assemblies | Chain assemblies are contiguous linkers between ring assemblies. They are uncovered by removing all ring bonds in the molecule |
| B | Number of Chains| Chains are all unbranched linkers needed to cover all nonring bonds in the molecule|
| C | Number of Rings | |
| D | Number of Ring Assemblies| Ring assemblies are fragments remaining when all acyclic bonds have been removed|
| E | Number of Bridge Bonds| |
||||
| F | Number of Ring Assemblies Consisting of Exactly One Ring | A contiguous path of more than one bond shared between more than one rings counts as bridge bond|
| G | Number of Ring Assemblies Consisting of Exactly Two Rings| |
| H | Number of Ring Assemblies Consisting of three or more than three Rings| |
| I | Number of Macrocycles | |
||||
| J | Binned Length of Shortest Chain | If the binned length of the shortest chain exists, it is used; otherwise, it is zero|
| K | Binned Length of Second Shortest Chain | If the binned length of the second shortest chain exists, it is used; otherwise, it is zero|
| L | Binned Length of Third Shortest Chain | If the binned length of the third shortest chain exists, it is used; otherwise, it is zero|
| M | Binned Length of Fourth Shortest Chain | If thebinned length of the fourth shortest chain exists, it is used; otherwise, it is zero|


#### Implementation
Generation of SCINS based on SMILES. The code can be found in [SCINS.py](SCINS.py)


```python
from SCINS import SCINS_generator

smi = "COc1ncc(-c2c(N)ncnc2N[C@@H](C)c2nn3ccc(C)c3c(=O)n2-c2ccccc2)cc1NS(=O)(=O)c1ccc(O)cc1"

generator = SCINS_generator(smi)

SCINS_vector = generator.Calculate_SCINS("vec")
SCINS_string = generator.Calculate_SCINS("str")
SCINS_string_formatted = generator.Calculate_SCINS("code")

print(SCINS_string_formatted)
```

#### Relevant papers

[1] Schuffenhauer, A., Brown, N., Ertl, P., Jenkins, J. L., & Selzer, P. (2007). Biological Activity Space. Journal of Chemical Information and Modeling, 325–336. https://doi.org/10.1021/ci6004004

[2] Bender, A., Jenkins, J. L., Scheiber, J., Sukuru, S. C. K., Glick, M., & Davies, J. W. (2009). How similar are similarity searching methods? A principal component analysis of molecular descriptor space. Journal of Chemical Information and Modeling, 49(1), 108–119. https://doi.org/10.1021/ci800249s
