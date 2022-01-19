# Copyright 2019-2022 Cambridge Quantum Computing
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Display a circuit as html."""

import os
from typing import Dict, Union, Optional
import jinja2
from pytket.circuit import Circuit  # type: ignore
from .utils import (
    format_register,
    format_mapping,
    format_raw_matrix,
    format_bool_matrix,
    format_logic_exp,
    get_gate_colour,
    has_gate_info,
    get_op_params,
    has_sub_circuit,
    get_sub_circuit,
    is_classical_gate,
    get_box_matrix,
    get_expbox_details,
    is_control_gate,
    get_target_args,
    get_op_name,
    get_op_display_name,
    format_op_params,
    parse_circuit,
)


# Set up jinja to access our templates
dirname = os.path.dirname(__file__)

loader = jinja2.FileSystemLoader(searchpath=dirname)
env = jinja2.Environment(loader=loader)

# Register the filters we need to use
env.filters["format_register"] = format_register
env.filters["format_mapping"] = format_mapping
env.filters["format_raw_matrix"] = format_raw_matrix
env.filters["format_bool_matrix"] = format_bool_matrix
env.filters["format_logic_exp"] = format_logic_exp
env.filters["get_gate_colour"] = get_gate_colour
env.filters["has_gate_info"] = has_gate_info
env.filters["get_op_params"] = get_op_params
env.filters["has_sub_circuit"] = has_sub_circuit
env.filters["get_sub_circuit"] = get_sub_circuit
env.filters["is_classical_gate"] = is_classical_gate
env.filters["get_box_matrix"] = get_box_matrix
env.filters["get_expbox_details"] = get_expbox_details
env.filters["is_control_gate"] = is_control_gate
env.filters["get_target_args"] = get_target_args
env.filters["get_op_name"] = get_op_name
env.filters["get_op_display_name"] = get_op_display_name
env.filters["format_op_params"] = format_op_params
env.filters["parse_circuit"] = parse_circuit


def render_circuit_as_html(
    circuit: Union[Dict[str, Union[str, float, dict]], Circuit],
    recursive: bool = False,
    condensed: bool = True,
    jupyter: bool = False,
) -> Optional[str]:
    """
    Render a circuit as HTML for inline display.

    :param circuit: the circuit to render.
    :param recursive: whether to display nested circuits as circuits themselves,
        or as generic boxes
    :param condensed: whether to render the circuit on one line only
        (may require scrolling),
        or allow the circuit to wrap around onto multiple lines
    :param jupyter: set to true to render generated HTML in cell output
    """
    if not isinstance(circuit, Circuit):
        circuit = Circuit.from_dict(circuit)

    options = {
        "recursive": recursive,
        "condensed": condensed,
    }

    template = env.get_template("circuit.html")
    html = template.render(
        {
            "circuit": circuit,
            "display_options": options,
        }
    )

    if jupyter:
        # If we are in a notebook, we can tell jupyter to display the html.
        # We don't import at the top in case we are not in a notebook environment.
        from IPython.core.display import (  # type: ignore
            HTML,
            display,
        )  # pylint: disable=C0415

        display(HTML(html))
        return None

    return html


def render_circuit_jupyter(
    circuit: Union[Dict[str, Union[str, float, dict]], Circuit],
    recursive: bool = False,
    condensed: bool = True,
) -> None:
    """Render a circuit as jupyter cell output.

    :param circuit: the circuit to render.
    :param recursive: whether to display nested circuits as circuits themselves,
        or as generic boxes
    :param condensed: whether to render the circuit on one line only
        (may require scrolling), or allow the circuit to wrap around
        onto multiple lines
    """
    render_circuit_as_html(circuit, recursive, condensed, True)
