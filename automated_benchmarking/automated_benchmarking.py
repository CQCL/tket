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
from pytket_benchmarking.compiler_benchmarking.experiment_analyser import CompiledCircuitExperiment
import matplotlib.pyplot as plt
from typing import Optional


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
    compiler: Annotated[Compilers, typer.Argument(help="Compiler to use.")],
    architecture: Annotated[Architectures, typer.Argument(help="Architecture to compile to.")],
    suite_path: Annotated[Path, typer.Argument(help="Path to circuit suite.")],
    compiled_path: Annotated[Path, typer.Argument(help="Path to compiled circuits.")],
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
def plot(
    suite_path: Annotated[Path, typer.Argument(help="Path to circuit suite.")],
    compiled_path: Annotated[Path, typer.Argument(help="Path to compiled circuits.")],
    show: Annotated[bool, typer.Option("--show", "-s", help='Display plot.')] = False,
    plot_path: Annotated[
        Optional[Path],
        typer.Argument(
            help=(
                "Path where plot should be saved. "
                + "If none is given it will not be saved."
            )
        )
    ] = None,
):
    storage_manager = LocalStorage(directory_path=suite_path)
    circuit_suite = CircuitSuiteManager.from_storage_manager(
        storage_manager=storage_manager
    )

    storage_manager = LocalStorage(directory_path=compiled_path)
    compiled_circuit_mgr = CompiledCircuitsManager.from_storage_manager(
        storage_manager=storage_manager,
    )

    compiled_func = lambda circuit : circuit.n_2qb_gates()
    original_func = lambda circuit : circuit.n_2qb_gates()

    benchmarking_experiment = CompiledCircuitExperiment(
        compiled_circuit_mgr_list = [compiled_circuit_mgr],
        circuit_suite_mgr_list = [circuit_suite],
    )

    fig, _ = benchmarking_experiment.plot_circuit_comparison(
        compiled_func=compiled_func,
        original_func=original_func,
    )
    if show:
        plt.show()

    if plot_path is not None:
        fig.savefig(fname=plot_path)

@app.command()
def return_test():
    return print(2)

if __name__ == "__main__":
    app()