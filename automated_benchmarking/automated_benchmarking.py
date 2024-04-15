import typer
from enum import Enum
from pytket_benchmarking.compiler_passes.qiskit_compilers import QiskitDefault
from pytket_benchmarking.compiler_passes.pytket_compilers import PytketDefault
from pytket_benchmarking.compiler_passes.architectures import HeavyHexagon
from pytket.passes import BasePass
from pytket.architecture import Architecture
from pytket import Circuit
from typing_extensions import Annotated
from pathlib import Path
from pytket_benchmarking.compiler_benchmarking.circuit_suite_manager import CircuitSuiteManager
from pytket_benchmarking.utils.storage_manager import LocalStorage
from pytket_benchmarking.compiler_benchmarking.pass_runner import TimedPassRunner
from pytket_benchmarking.compiler_benchmarking import CompiledCircuitsManager

app = typer.Typer(
    # pretty_exceptions_enable=False,
)

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
    compiler: Annotated[Compilers, typer.Argument(help="Compiler to use.")],
    architecture: Annotated[Architectures, typer.Argument(help="Architecture to compile to.")],
    suite_path: Annotated[Path, typer.Argument(help="Path to circuit suite.")],
    compiled_path: Annotated[Path, typer.Argument(help="Path to compiled circuits.")]
):
    """Compile circuit
    """

    storage_manager = LocalStorage(directory_path=suite_path)
    circuit_suite = CircuitSuiteManager.from_storage_manager(
        storage_manager=storage_manager
    )

    compiler_pass=compilers_dict[compiler](
        architecture=architectures_dict[architecture],
        optimisation_level=0,
    )
    pass_runner = TimedPassRunner(
        compiler_pass=compiler_pass,
        label=compiler.value,
    )

    storage_manager = LocalStorage(
        directory_path=compiled_path,
        write_permission=True,
    )
    compiled_circuit_mgr = CompiledCircuitsManager(
        storage_manager=storage_manager,
    )

    for original_circuit in circuit_suite:
        compiled_circuit_mgr.run_circuit(
            pass_runner=pass_runner,
            original_circuit=original_circuit,
        )

@app.command()
def plot():
    pass

if __name__ == "__main__":
    app()