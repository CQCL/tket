# Project Description

pytket is a python module for interfacing with tket - an optimising compiler for quantum circuits. In addition to pytket there are several extension modules for accessing a range of quantum hardware and classical simulators. The extension modules also provide integration with several widely used quntum software tools.

The source code for the tket compiler can be found in the relevant [github repository](https://github.com/CQCL/tket).

## Installation

Installation is supported for Linux, MacOS and Windows. Installation requires python 3.8 or newer.

To install run the pip command: 

`` pip install pytket``

See [Installation troubleshooting](https://cqcl.github.io/tket/pytket/api/install.html) for help with installation.

To install the pytket extension modules add a hyphen and the extension name to the command:

`` pip install pytket-quantinuum ``

For a complete list of pytket extensions see this page: https://cqcl.github.io/pytket-extensions/api/index.html

## Documentaion and Examples

API reference: https://cqcl.github.io/tket/pytket/api/getting_started.html

To get started using pytket see our [user manual](https://cqcl.github.io/pytket/manual/index.html).

Finally for worked examples using tket see our github [examples repository](https://github.com/CQCL/pytket/tree/main/examples).




## Support and Discussion

For bugs and feature requests we recommend creating an issue on the [github repository](https://github.com/CQCL/pytket).

User support: tket-support@cambridgequantum.com

Discussion: Join our public slack channel [here](https://join.slack.com/t/tketusers/shared_invite/zt-18qmsamj9-UqQFVdkRzxnXCcKtcarLRA).

Mailing list: join [here](https://list.cambridgequantum.com/cgi-bin/mailman/listinfo/tket-users).

## Citation

If you wish to cite tket in any academic publications, we generally recommend citing our [software overview](https://arxiv.org/abs/2003.10611) paper for most cases.

If your work is on the topic of specific compilation tasks, it may be more appropriate to cite one of our other papers:

- "On the qubit routing problem" for qubit placement (a.k.a. allocation) and routing (a.k.a. swap network insertion, connectivity solving).
- "Phase Gadget Synthesis for Shallow Circuits" for representing exponentiated Pauli operators in the ZX calculus and their circuit decompositions.
- "A Generic Compilation Strategy for the Unitary Coupled Cluster Ansatz" for sequencing of terms in Trotterisation and Pauli diagonalisation.

We are also keen for others to benchmark their compilation techniques against us. We recommend checking our benchmark repository for examples on how to run basic benchmarks with the latest version of pytket. Please list the release version of pytket with any benchmarks you give, and feel free to get in touch for any assistance needed in setting up fair and representative tests.

## Licence

pytket is distributed under the apache licence 2.0. To read the licence click [here](https://github.com/CQCL/pytket/blob/main/LICENCE).
