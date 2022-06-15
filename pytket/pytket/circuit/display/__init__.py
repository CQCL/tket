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
import uuid
import json
import jinja2
from pytket.circuit import Circuit  # type: ignore


# Set up jinja to access our templates
dirname = os.path.dirname(__file__)

loader = jinja2.FileSystemLoader(searchpath=dirname)
env = jinja2.Environment(loader=loader)


def render_circuit_as_html(
    circuit: Union[Dict[str, Union[str, float, dict]], Circuit],
    jupyter: bool = False,
) -> Optional[str]:
    """
    Render a circuit as HTML for inline display.

    :param circuit: the circuit to render.
    :param jupyter: set to true to render generated HTML in cell output.
    """
    if not isinstance(circuit, Circuit):
        circuit = Circuit.from_dict(circuit)

    uid = uuid.uuid4()
    html_template = env.get_template("circuit.html")
    html = html_template.render(
        {
            "circuit_json": json.dumps(circuit.to_dict()),
            "uid": uid,
            "jupyter": jupyter,
        }
    )
    if jupyter:
        # If we are in a notebook, we can tell jupyter to display the html.
        # We don't import at the top in case we are not in a notebook environment.
        from IPython.display import (  # type: ignore
            HTML,
            display,
        )  # pylint: disable=C0415

        display(HTML(html))
        return None

    return html


def render_circuit_jupyter(
    circuit: Union[Dict[str, Union[str, float, dict]], Circuit],
) -> None:
    """Render a circuit as jupyter cell output.

    :param circuit: the circuit to render.
    """
    render_circuit_as_html(circuit, True)
