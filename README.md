[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)

<p align="center">
  <a href="https://mqt.readthedocs.io">
   <picture>
     <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/munich-quantum-toolkit/.github/refs/heads/main/docs/_static/mqt-banner-dark.svg" width="90%">
     <img src="https://raw.githubusercontent.com/munich-quantum-toolkit/.github/refs/heads/main/docs/_static/mqt-banner-light.svg" width="90%" alt="MQT Banner">
   </picture>
  </a>
</p>

# MQT IonShuttler

<i>MQT IonShuttler</i> is a tool for generating shuttling schedules for grid-like trapped-ion quantum computers based on the Quantum CCD (QCCD) architecture. It supports

* **exact shuttling schedules** for small architectures with **a single processing zone (PZ)**, and
* **heuristic shuttling schedules** for larger devices with **one _or_ multiple processing zones**.

<p align="center">
  <figure>
    <object data="docs/figures/QCCD_device.pdf" type="application/pdf" width="48%"></object>
    <object data="docs/figures/graph.pdf" type="application/pdf" width="48%"></object>
    <figcaption><b>Figure&nbsp;1:</b> (<b>a</b>) Potential QCCD device with four processing zones; (<b>b</b>) corresponding graph used by <i>MQT IonShuttler</i>.</figcaption>
  </figure>
</p>

The exact solution guarantees optimality but is limited to a single PZ, while the heuristic method scales to many qubits and PZs. In the heuristic workflow, an optional **compilation** feature allows for dynamic rescheduling of gates based on the current ion positions and dependencies (leveraging Qiskit’s `DAGDependency`), potentially reducing shuttling overhead compared to executing a fixed sequence.

<i>MQT IonShuttler</i> is part of the [_Munich Quantum Toolkit_](https://mqt.readthedocs.io) (MQT) developed by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/).

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by creating an [issue](https://github.com/cda-tum/mqt-ion-shuttler/issues) on GitHub.

---

## Getting Started

We strongly recommend using [virtual environments](https://docs.python.org/3/library/venv.html) to set up the tool and install the dependencies.

### Exact Solution (single PZ)
```bash
python3 -m venv .venv
. .venv/bin/activate
pip install .
python3 run.py --help
# Example
python3 run.py algorithms/qft_05.json
```

The script supports an additional `--plot` argument to visualise the result. All parameters are documented via `--help`. Architectures and algorithms are specified in JSON files—see the examples in [`algorithms/`](algorithms/).

### Heuristic Solution (single & multiple PZs)
```bash
python3 -m venv .venv
. .venv/bin/activate
pip install .
python3 run_heuristic.py --help
# Example with 12 logical qubits executed on four PZs
python3 run_heuristic.py algorithms_heuristic/qft_12.json
```

In addition to the parameters shown by `--help`, the heuristic flow offers the optional `--compile` switch that triggers the dynamic compilation step described above.

Architectures and algorithms are specified in JSON files—see the examples in [`algorithms_heuristic/`](algorithms_heuristic/).

---

## References

This implementation is based on the following publications:

1. D. Schoenberger, S. Hillmich, M. Brandl, and R. Wille, “Using Boolean Satisfiability for Exact Shuttling in Trapped-Ion Quantum Computers,” *ASP-DAC 2024*.
2. D. Schoenberger, S. Hillmich, M. Brandl, and R. Wille, “Shuttling for Scalable Trapped-Ion Quantum Computers,” *IEEE TCAD*, 2024.

---

## Acknowledgements

The Munich Quantum Toolkit has been supported by the European Union’s Horizon 2020 research and innovation programme (DA QC, grant agreement No. 101001318 and MILLENION, grant agreement No. 101114305), the Bavarian State Ministry for Science and Arts through the Distinguished Professorship Program, as well as the Munich Quantum Valley, which is supported by the Bavarian state government with funds from the Hightech Agenda Bayern Plus.

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/munich-quantum-toolkit/.github/refs/heads/main/docs/_static/mqt-funding-footer-dark.svg" width="90%">
    <img src="https://raw.githubusercontent.com/munich-quantum-toolkit/.github/refs/heads/main/docs/_static/mqt-funding-footer-light.svg" width="90%" alt="MQT Funding Footer">
  </picture>
</p>
