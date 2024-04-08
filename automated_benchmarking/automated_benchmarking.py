import typer
from enum import Enum
from pytket_benchmarking.compiler_passes.qiskit_compilers import QiskitDefault
from pytket_benchmarking.compiler_passes.pytket_compilers import PytketDefault
from pytket_benchmarking.compiler_passes.architectures import HeavyHexagon
from pytket.passes import BasePass
from pytket.architecture import Architecture
from pytket import Circuit
from typing_extensions import Annotated

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
def compile(
    compiler: Annotated[Compilers, typer.Argument(help="Compiler to use")],
    architecture: Annotated[Architectures, typer.Argument(help="Architecture to compile to")]
):
    """Compile circuit
    """
    circuit = Circuit(2)
    compilers_dict[compiler](
        architecture=architectures_dict[architecture],
        optimisation_level=0,
    ).apply(circuit)

if __name__ == "__main__":
    app()