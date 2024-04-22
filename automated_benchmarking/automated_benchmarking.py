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
from typing import Optional, List
from pytket.passes import SequencePass, DecomposeBoxes
from rich.progress import track


app = typer.Typer()

class Compilers(Enum):

    QiskitIBMQ = "QiskitIBMQ"
    PytketIBMQ = "PytketIBMQ"

ibmq_architecture = HeavyHexagon(n_rows=3, n_columns=3)

compilers_dict = {
    Compilers.QiskitIBMQ: SequencePass(
        pass_list=[
            DecomposeBoxes(),
            QiskitDefault(
                architecture=ibmq_architecture,
                optimisation_level=3,
            )
        ]
    ),
    Compilers.PytketIBMQ: PytketDefault(architecture=ibmq_architecture, optimisation_level=2),
}

@app.command()
def compile(
    compiler: Annotated[Compilers, typer.Argument(help="Compiler to use.")],
    suite_path: Annotated[Path, typer.Argument(help="Path to circuit suite.")],
    compiled_path: Annotated[Path, typer.Argument(help="Path to compiled circuits.")],
):
    """Compile circuit
    """

    storage_manager = LocalStorage(directory_path=suite_path)
    circuit_suite = CircuitSuiteManager.from_storage_manager(
        storage_manager=storage_manager
    )

    pass_runner = TimedPassRunner(
        compiler_pass=compilers_dict[compiler],
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
    plot_path: Annotated[
        Path,
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

    fig.savefig(fname=plot_path)

@app.command()
def percentage_better(
    circuit_suite_path: Annotated[Path, typer.Argument(help="Path circuit suite.")],
    compiled_path: Annotated[Path, typer.Argument(help="Path to compiled circuits.")],
    compilers: Annotated[List[Compilers], typer.Argument(help="Compilers to compare.")],
):    

    storage_manager = LocalStorage(
        directory_path=compiled_path,
        write_permission=True,
    )
    compiled_circuit_mgr = CompiledCircuitsManager(
        storage_manager=storage_manager,
    )

    storage_manager = LocalStorage(directory_path=circuit_suite_path)
    circuit_suite_mgr = CircuitSuiteManager(
        storage_manager=storage_manager
    )

    for compiler in compilers:

        pass_runner = TimedPassRunner(
            compiler_pass=compilers_dict[compiler],
            label=compiler.value,
        )
        
        for original_circuit in track(
            circuit_suite_mgr,
            description="Processing..."
        ):
            compiled_circuit_mgr.run_circuit(
                pass_runner=pass_runner,
                original_circuit=original_circuit,
            )
    
    benchmarking_experiment = CompiledCircuitExperiment(
        compiled_circuit_mgr_list = [compiled_circuit_mgr],
        circuit_suite_mgr_list = [circuit_suite_mgr],
    )

    data = benchmarking_experiment.compiler_comparison(
        comparison_func=lambda circuit : circuit.n_2qb_gates(),
    )

    data['pytket better'] = data['PytketIBMQ'] <= data['QiskitIBMQ']
    print(f"pytket is okay at {100 * len(data.loc[data['pytket better']]) / len(data)}")


@app.command()
def return_test():
    return print(2)

if __name__ == "__main__":
    app()