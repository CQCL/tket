import inspect
import re

from pytket._tket.circuit import Circuit

# This script is used to check that no methods or properties of
# the Circuit class are left out in the circuit_class.rst file.


class MissingCircuitDocsError(Exception):
    "Raised when the Circuit class has missing methods or properties."

    def __init__(self, missing_list: list[str]) -> None:
        super().__init__(
            f"Missing Circuit class docs for {missing_list}! Please"
            " update circuit_class.rst.",
        )


def main() -> None:
    with open("circuit_class.rst", encoding="utf-8") as file:
        rst = file.read()

    in_methods = {p[0] for p in inspect.getmembers(Circuit)}
    pattern = re.compile(r"(?:\bautomethod\b|\bautoproperty\b):: (\w+)\n")
    in_docs = set(pattern.findall(rst))
    missing = in_methods - in_docs
    missing_public_docs = [s for s in missing if not s.startswith("_")]
    if len(missing_public_docs) != 0:
        raise MissingCircuitDocsError(missing_list=missing_public_docs)


if __name__ == "__main__":
    main()
