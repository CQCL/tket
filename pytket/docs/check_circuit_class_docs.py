import inspect
import re

from pytket._tket.circuit import Circuit


with open("circuit_class.rst", "r", encoding="utf-8") as file:
    rst = file.read()

in_methods = {p[0] for p in inspect.getmembers(Circuit)}
pattern = re.compile(r"(?:\bautomethod\b|\bautoproperty\b):: (\w+)\n")
in_docs = set(pattern.findall(rst))
missing = in_methods - in_docs
print([s for s in missing if not s.startswith("_")])
