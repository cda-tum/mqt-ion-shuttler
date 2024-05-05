[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)

<p align="center">
  <a href="https://mqt.readthedocs.io">
   <picture>
     <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cda-tum/mqt/main/docs/_static/mqt_light.png" width="60%">
     <img src="https://raw.githubusercontent.com/cda-tum/mqt/main/docs/_static/mqt_dark.png" width="60%">
   </picture>
  </a>
</p>

# MQT IonShuttler

This a tool for generating shuttling schedules for grid-like memory zones inside trapped-ion quantum computers. The tool supports optimal shuttling schedules for small architectures and also includes a preliminary option to produce heuristic shuttling schedules for large devices. The optimal solution is based on the paper *Using Boolean Satisfiability for Exact Shuttling in Trapped-Ion Quantum Computers* by D. Schoenberger, S. Hillmich, M. Brandl, and R. Wille in ASP-DAC 2024. The heuristic solution is based on the paper *Shuttling for Scalable Trapped-Ion Quantum Computers* by D. Schoenberger, S. Hillmich, M. Brandl, and R. Wille.

MQT IonShuttler is part of the [_Munich Quantum Toolkit_](https://mqt.readthedocs.io) (_MQT_) developed by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/).

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by creating an [issue](#) on GitHub.

## Getting started

We strongly recommend using [virtual environments](https://docs.python.org/3/library/venv.html) to set up the tool and install the dependencies

# Exact Solution
```commandline
$ python3 -m venv .venv
$ . .venv/bin/activate
(.venv) $ pip install .
(.venv) $ python3 run.py --help
[...]
(.venv) $ python3 run.py algorithms/qft_05.json
...
```

The script supports an additional `--plot` argument to plot the result.
An overview over all parameters is printed with the `--help` parameter.
The architecture and the algorithm are specified in json files.
You can find examples in the [`algorithms/`](algorithms/) folder.

# Heuristic Solution
```commandline
$ python3 -m venv .venv
$ . .venv/bin/activate
(.venv) $ pip install .
(.venv) $ python3 run_heuristic.py --help
[...]
(.venv) $ python3 run_heuristic.py algorithms_heuristic/qft_24.json
...
```

An overview over all parameters is printed with the `--help` parameter.
The architecture and the algorithm are specified in json files.
You can find examples in the [`algorithms_heuristic/`](algorithms_heuristic/) folder.

## Acknowledgements

The Munich Quantum Toolkit has been supported by the European
Research Council (ERC) under the European Union's Horizon 2020 research and innovation program (grant agreement
No. 101001318), the Bavarian State Ministry for Science and Arts through the Distinguished Professorship Program, as well as the
Munich Quantum Valley, which is supported by the Bavarian state government with funds from the Hightech Agenda Bayern Plus.

<p align="center">
<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cda-tum/mqt/main/docs/_static/tum_dark.svg" width="28%">
<img src="https://raw.githubusercontent.com/cda-tum/mqt/main/docs/_static/tum_light.svg" width="28%" alt="TUM Logo">
</picture>
<picture>
<img src="https://raw.githubusercontent.com/cda-tum/mqt/main/docs/_static/logo-bavaria.svg" width="16%" alt="Coat of Arms of Bavaria">
</picture>
<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cda-tum/mqt/main/docs/_static/erc_dark.svg" width="24%">
<img src="https://raw.githubusercontent.com/cda-tum/mqt/main/docs/_static/erc_light.svg" width="24%" alt="ERC Logo">
</picture>
<picture>
<img src="https://raw.githubusercontent.com/cda-tum/mqt/main/docs/_static/logo-mqv.svg" width="28%" alt="MQV Logo">
</picture>
</p>
