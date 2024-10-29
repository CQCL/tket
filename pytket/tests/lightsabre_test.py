from pytket.passes import LightSABRE
from pytket import Circuit
from pytket.passes import BasePass
from pytket.architecture import Architecture

if __name__ == "__main__":
    a = Architecture([(0, 1), (1, 2)])
    c = Circuit(3).CX(0, 2).CX(0, 1)
    d = LightSABRE(a, 45, 2).to_dict()
    print(d)
    bp = BasePass.from_dict(d)

    bp.apply(c)
    for x in c:
        print(x)
