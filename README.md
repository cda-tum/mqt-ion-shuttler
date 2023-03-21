[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)

# MQT IonShuttler

This a tool for generating optimal shuttling schedules for grid-like memory zones inside trapped-ion quantum computers.
It is part of the Munich Quantum Toolkit (MQT; formerly known as JKQ and developed by the [Institute for Integrated Circuits](https://iic.jku.at/eda/) at the [Johannes Kepler University Linz](https://jku.at)).

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
