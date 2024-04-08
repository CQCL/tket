import typer
from enum import Enum
from pytket_benchmarking.compiler_passes.qiskit_compilers import QiskitDefault
from pytket_benchmarking.compiler_passes.pytket_compilers import PytketDefault
from pytket_benchmarking.compiler_passes.architectures import HeavyHexagon
from pytket.passes import BasePass
from pytket.architecture import Architecture
from pytket import Circuit

app = typer.Typer()

class Compilers(Enum):

    QiskitDefault = "Qiskit"
    PytketDefault = "Pytket"

compilers_dict = {
    Compilers.QiskitDefault: QiskitDefault,
    Compilers.PytketDefault: PytketDefault,
}

class Architectures(Enum):

    heavy_hexagon = "HeavyHexagon"

architectures_dict = {
    Architectures.heavy_hexagon: HeavyHexagon(n_rows=10, n_columns=10)
}

@app.command()
def compile(compiler: Compilers, architecture: Architectures):
    circuit = Circuit(2)
    compilers_dict[compiler](
        architecture=architectures_dict[architecture],
        optimisation_level=0,
    ).apply(circuit)


@app.command()
def hello(name: str):
    print(f"Hello {name}")


@app.command()
def goodbye(name: str, formal: bool = False):
    if formal:
        print(f"Goodbye Ms. {name}. Have a good day.")
    else:
        print(f"Bye {name}!")


if __name__ == "__main__":
    app()