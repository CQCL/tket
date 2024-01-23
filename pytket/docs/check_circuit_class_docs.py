import inspect
import re

from pytket._tket.circuit import Circuit

# paste contents of https://github.com/CQCL/tket/blob/7b0df5cf68503944cd34e593fb094ba5aa4357f4/pytket/docs/circuit_class.rst


with open("circuit_class.rst", "r", encoding="") as file:
    rst = file.read()

in_methods = {p[0] for p in inspect.getmembers(Circuit)}
pattern = re.compile(r"(?:\bautomethod\b|\bautoproperty\b):: (\w+)\n")
in_docs = set(pattern.findall(rst))
missing = in_methods - in_docs
print([s for s in missing if not s.startswith("_")])
