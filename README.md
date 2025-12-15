# Branch-and-Price for the k-Vertex Cut Problem
[![DOI](https://zenodo.org/badge/1069071398.svg)](https://doi.org/10.5281/zenodo.17897865)



<div align="center">
   <img src="./logo.png" alt="Logo">
</div>


This repository contains the implementation and experimental data for the branch-and-price algorithm described in:

> **Branch-and-price strikes back for the k-vertex cut problem**  
> F. Ciccarelli, F. Furini, C. Hojny, and M. Lübbecke


## 📁 Repository structure

### [`branch-and-price/`](branch-and-price/)
SCIP-based branch-and-price solver implementation. Contains source code, computational results and comprehensive documentation.

See [`branch-and-price/README.md`](branch-and-price/README.md) for detailed usage instructions, parameter descriptions, and output formats.

### [`data/`](data/)
Benchmark instances for the k-vertex cut problem:
- `raw_instances/`: original graph instances in DIMACS format.
- `preprocessed_instances/`: instances after preprocessing with fixed vertices removed.

### [`results/`](results/)
Computational results from experimental runs:
- `BPresults.xlsx`: complete branch-and-price computational results.
- `bestKnownValues.xlsx`: best known solutions for benchmark instances.

## 🛠️ Requirements

### Software dependencies
- **[SCIP](https://www.scipopt.org/)** (≥ 8.0): optimization framework for branch-and-price.
- **[Gurobi](https://www.gurobi.com/)**: used for preprocessing of the instances and ILP warm start.
- **[LEMON](https://lemon.cs.elte.hu/)**: graph library used for solving max flow problems.
- **[Graphviz](https://graphviz.org/)** (optional): for solution visualization.

### Build dependencies
- C++17 compatible compiler
- Make

## 🚀 Getting started

1. **Install dependencies**: ensure SCIP, Gurobi, LEMON, and optionally Graphviz are installed and accessible.

2. **Build the solver**:
   ```bash
   cd branch-and-price
   make
   ```

3. **Run on an instance**:
   ```bash
   ./bin/kvertexcut -f ../data/raw_instances/instance.dimacs -s params/params_k10.set
   ```

4. **View results**: check the generated `.res` and `.sol` files in the configured output directories.

For detailed parameter configuration and output format, see [`branch-and-price/README.md`](branch-and-price/README.md).

## 📄 License

This software is released for **academic and research purposes only**. See [LICENSE](LICENSE) for full terms.

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{ciccarelli2024branch,
  title={Branch-and-price strikes back for the k-vertex cut problem},
  author={Ciccarelli, Fabio and Furini, Fabio and Hojny, Christopher and L{\"u}bbecke, Marco},
  year={2024}
}
```

## 📧 Contact

For questions or issues, please contact: [f.ciccarelli@uniroma1.it](mailto:f.ciccarelli@uniroma1.it)
