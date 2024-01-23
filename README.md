[![PyPI](https://img.shields.io/pypi/v/mqt.ddsim?logo=pypi&style=flat-square)](https://pypi.org/project/mqt.ddsim/)
![OS](https://img.shields.io/badge/os-linux%20%7C%20macos%20%7C%20windows-blue?style=flat-square)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://img.shields.io/github/actions/workflow/status/cda-tum/mqt-ddsim/ci.yml?branch=main&style=flat-square&logo=github&label=ci)](https://github.com/cda-tum/mqt-ddsim/actions/workflows/ci.yml)
[![CD](https://img.shields.io/github/actions/workflow/status/cda-tum/mqt-ddsim/cd.yml?style=flat-square&logo=github&label=cd)](https://github.com/cda-tum/mqt-ddsim/actions/workflows/cd.yml)
[![Documentation](https://img.shields.io/readthedocs/ddsim?logo=readthedocs&style=flat-square)](https://mqt.readthedocs.io/projects/ddsim)
[![codecov](https://img.shields.io/codecov/c/github/cda-tum/mqt-ddsim?style=flat-square&logo=codecov)](https://codecov.io/gh/cda-tum/mqt-ddsim)

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cda-tum/mqt-ddsim/main/docs/source/_static/mqt_light.png" width="60%">
    <img src="https://raw.githubusercontent.com/cda-tum/mqt-ddsim/main/docs/source/_static/mqt_dark.png" width="60%">
  </picture>
</p>

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
