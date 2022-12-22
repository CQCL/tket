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

import json
import os
import tempfile
import time
import uuid
import webbrowser
from typing import Dict, Optional, Union, cast

import jinja2

from pytket.circuit import Circuit  # type: ignore

# Set up jinja to access our templates
dirname = os.path.dirname(__file__)

loader = jinja2.FileSystemLoader(searchpath=dirname)
env = jinja2.Environment(loader=loader)

RenderCircuit = Union[Dict[str, Union[str, float, dict]], Circuit]


def render_circuit_as_html(
    circuit: RenderCircuit,
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
    circuit: RenderCircuit,
) -> None:
    """Render a circuit as jupyter cell output.

    :param circuit: the circuit to render.
    """
    render_circuit_as_html(circuit, True)


def view_browser(circuit: RenderCircuit, browser_new: int = 2, sleep: int = 5) -> None:
    """Write circuit render html to a tempfile and open in browser.

    Waits for some time for browser to load then deletes tempfile.

    :param circuit: the Circuit or serialized Circuit to render.
    :param browser_new: ``new`` parameter to ``webbrowser.open``, default 2.
    :param sleep: Number of seconds to sleep before deleting file, default 5.

    """

    fp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".html", delete=False, dir=os.getcwd()
    )
    try:
        fp.write(cast(str, render_circuit_as_html(circuit)))
        fp.close()

        webbrowser.open("file://" + os.path.realpath(fp.name), new=browser_new)

        # give browser enough time to open before deleting file
        time.sleep(sleep)
    finally:
        os.remove(fp.name)
