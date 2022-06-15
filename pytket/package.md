Pytket is a python module for interfacing with TKET, an optimising compiler for quantum circuits developed by Cambridge Quantum. In addition to pytket there are several extension modules for accessing a range of quantum hardware and classical simulators. The extension modules also provide integration with several widely used quantum software tools.

The source code for the TKET compiler can be found in [this github repository](https://github.com/CQCL/tket).

## Installation

Installation is supported for Linux, MacOS and Windows. Installation requires python 3.8, 3.9 or 3.10.

To install run the pip command: 

`` pip install pytket``

See [Installation troubleshooting](https://cqcl.github.io/tket/pytket/api/install.html) for help with installation.

To install the pytket extension modules add a hyphen and the extension name to the command:

`` pip install pytket-quantinuum ``

For a list of pytket extensions see this page: https://cqcl.github.io/pytket-extensions/api/index.html.

## Documentation and Examples

API reference: https://cqcl.github.io/tket/pytket/api/

To get started using pytket see the [user manual](https://cqcl.github.io/pytket/manual/index.html).

For worked examples using TKET see our [examples repository](https://github.com/CQCL/pytket/tree/main/examples).

## Support and Discussion

For bugs and feature requests we recommend creating an issue on the [github repository](https://github.com/CQCL/tket).

User support: tket-support@cambridgequantum.com

For discussion, join the public slack channel [here](https://join.slack.com/t/tketusers/shared_invite/zt-18qmsamj9-UqQFVdkRzxnXCcKtcarLRA).

Mailing list: join [here](https://list.cambridgequantum.com/cgi-bin/mailman/listinfo/tket-users).

## Citation

If you wish to cite TKET in any academic publications, we generally recommend citing our [software overview](https://arxiv.org/abs/2003.10611) paper for most cases.

If your work is on the topic of specific compilation tasks, it may be more appropriate to cite one of our other papers:

- "On the qubit routing problem" for qubit placement (a.k.a. allocation) and routing (a.k.a. swap network insertion, connectivity solving). https://arxiv.org/abs/1902.08091 .
- "Phase Gadget Synthesis for Shallow Circuits" for representing exponentiated Pauli operators in the ZX calculus and their circuit decompositions. https://arxiv.org/abs/1906.01734 .
- "A Generic Compilation Strategy for the Unitary Coupled Cluster Ansatz" for sequencing of terms in Trotterisation and Pauli diagonalisation. https://arxiv.org/abs/2007.10515 .