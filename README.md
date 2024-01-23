[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cda-tum/mqt-ion-shuttler/main/docs/source/_static/mqt_light.png" width="60%">
    <img src="https://raw.githubusercontent.com/cda-tum/mqt-ion-shuttler/main/docs/source/_static/mqt_dark.png" width="60%">
  </picture>
</p>

# MQT IonShuttler

This a tool for generating optimal shuttling schedules for grid-like memory zones inside trapped-ion quantum computers, based on the paper *Using Boolean Satisfiability for Exact Shuttling in Trapped-Ion Quantum Computers* by D. Schoenberger, S. Hillmich, M. Brandl, and R. Wille in ASP-DAC 2024.

MQT IonShuttler is part of the Munich Quantum Toolkit (MQT) developed by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/).

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by creating an [issue](#) on GitHub.

## Getting started

We strongly recommend using [virtual environments](https://docs.python.org/3/library/venv.html) to set up the tool and install the dependencies

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
