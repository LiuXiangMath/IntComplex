# IntComplex for high-order interactions

---
## Installation
1. Create the environment

```
conda create -n intcomplex python=3.11 -y
conda activate intcomplex
```
2. Install required packages

```
pip install -r requirements.txt
```

---
## Usage
Run the following command to try an example
```
python ./code/example.py
```



---
## Protein--Ligand Binding Affinity Prediction
To apply the IntComplex framework to protein--ligand complexes and perform the binding affinity prediction task, run the following command. Before execution, please download the PDBbind dataset from https://www.pdbbind-plus.org.cn
```
python ./code/pdbbind.py
```

---

## Citation
```shell
# BibTex
@article{liu2024intcomplex,
  title={Intcomplex for high-order interactions},
  author={Liu, Xiang and Liu, Ran and Li, Jingyan and Wu, Rongling and Wu, Jie},
  journal={arXiv preprint arXiv:2412.02806},
  year={2024}
}
```
