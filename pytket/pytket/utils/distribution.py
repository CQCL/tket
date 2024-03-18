# Copyright 2019-2024 Cambridge Quantum Computing
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

from collections import defaultdict
from typing import (
    Any,
    Callable,
    DefaultDict,
    Dict,
    List,
    Set,
    Tuple,
    Union,
    Generic,
    TypeVar,
    Counter,
)
import warnings
import numpy as np
from scipy.stats import rv_discrete

Number = Union[float, complex]
T0 = TypeVar("T0")
T1 = TypeVar("T1")


class EmpiricalDistribution(Generic[T0]):
    """Represents an empirical distribution of values.

    Supports methods for combination, marginalization, expectation value, etc.

    >>> dist1 = EmpiricalDistribution(Counter({(0, 0): 3, (0, 1): 2, (1, 0): 4, (1, 1):
    ... 0}))
    >>> dist2 = EmpiricalDistribution(Counter({(0, 0): 1, (0, 1): 0, (1, 0): 2, (1, 1):
    ... 1}))
    >>> dist1.sample_mean(lambda x : x[0] + 2*x[1])
    0.8888888888888888
    >>> dist3 = dist2.condition(lambda x: x[0] == 1)
    >>> dist3
    EmpiricalDistribution(Counter({(1, 0): 2, (1, 1): 1}))
    >>> dist4 = dist1 + dist3
    >>> dist4
    EmpiricalDistribution(Counter({(1, 0): 6, (0, 0): 3, (0, 1): 2, (1, 1): 1}))
    """

    def __init__(self, C: Counter[T0]):
        self._C: Counter[T0] = Counter({x: c for x, c in C.items() if c > 0})

    def as_counter(self) -> Counter[T0]:
        """Return the distribution as a :py:class:`collections.Counter` object."""
        return self._C

    @property
    def total(self) -> int:
        """Return the total number of observations."""
        return sum(self._C.values())  # Counter.total() new in 3.10

    @property
    def support(self) -> Set[T0]:
        """Return the support of the distribution (set of all observations)."""
        return set(self._C.keys())

    def __eq__(self, other: object) -> bool:
        """Compare distributions for equality."""
        if not isinstance(other, EmpiricalDistribution):
            return NotImplemented
        return self._C == other._C

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self._C)})"

    def __getitem__(self, x: T0) -> int:
        """Get the count associated with an observation."""
        return self._C[x]

    def __add__(
        self, other: "EmpiricalDistribution[T0]"
    ) -> "EmpiricalDistribution[T0]":
        """Combine two distributions."""
        return EmpiricalDistribution(self._C + other._C)

    def condition(self, criterion: Callable[[T0], bool]) -> "EmpiricalDistribution[T0]":
        """Return a new distribution conditioned on the given criterion.

        :param criterion: A boolean function defined on all possible observations.
        """
        return EmpiricalDistribution(
            Counter({x: c for x, c in self._C.items() if criterion(x)})
        )

    def map(self, mapping: Callable[[T0], T1]) -> "EmpiricalDistribution[T1]":
        """Return a distribution over a transformed domain.

        The provided function maps elements in the original domain to new elements. If
        it is not injective, counts are combined.

        :param mapping: A function defined on all possible observations, mapping them
            to another domain.
        """
        C: Counter[T1] = Counter()
        for x, c in self._C.items():
            C[mapping(x)] += c
        return EmpiricalDistribution(C)

    def sample_mean(self, f: Callable[[T0], Number]) -> Number:
        """Compute the sample mean of a functional.

        The provided function maps observations to numerical values.

        :return: Estimate of the mean of the functional based on the observations."""
        return sum(c * f(x) for x, c in self._C.items()) / self.total

    def sample_variance(self, f: Callable[[T0], Number]) -> Number:
        """Compute the sample variance of a functional.

        The provided function maps observations to numerical values.

        The sample variance is an unbiased estimate of the variance of the underlying
        distribution.

        :return: Estimate of the variance of the functional based on the
            observations."""
        if self.total < 2:
            raise RuntimeError(
                "At least two samples are required in order to compute the sample "
                "variance."
            )
        fs = [(f(x), c) for x, c in self._C.items()]
        M0 = self.total
        M1 = sum(c * v for v, c in fs)
        M2 = sum(c * v**2 for v, c in fs)
        return (M2 - M1**2 / M0) / (M0 - 1)


class ProbabilityDistribution(Generic[T0]):
    """Represents an exact probability distribution.

    Supports methods for combination, marginalization, expectation value, etc. May be
    derived from an :py:class:`EmpriricalDistribution`.
    """

    def __init__(self, P: Dict[T0, float], min_p: float = 0.0):
        """Initialize with a dictionary of probabilities.

        :param P: Dictionary of probabilities.
        :param min_p: Optional probability below which to ignore values. Default
            0. Distribution is renormalized after removing these values.

        The values must be non-negative. If they do not sum to 1, a warning is
        emitted; the distribution will contain normalized values.
        """
        if any(x < 0 for x in P.values()):
            raise ValueError("Distribution contains negative probabilities")
        s0 = sum(P.values())
        if np.isclose(s0, 0):
            raise ValueError("Distribution has zero weight")
        if not np.isclose(s0, 1):
            warnings.warn(
                "Probabilities used to initialize ProbabilityDistribution do "
                "not sum to 1: renormalizing."
            )
        newP = {x: p for x, p in P.items() if p > min_p}
        s = sum(newP.values())
        self._P = {x: p / s for x, p in newP.items()}

    def as_dict(self) -> Dict[T0, float]:
        """Return the distribution as a :py:class:`dict` object."""
        return self._P

    def as_rv_discrete(self) -> Tuple[rv_discrete, List[T0]]:
        """Return the distribution as a :py:class:`scipy.stats.rv_discrete` object.

        This method returns an RV over integers {0, 1, ..., k-1} where k is the size of
        the support, and a list whose i'th member is the item corresponding to the value
        i of the RV.
        """
        X = list(self._P.keys())
        return (rv_discrete(values=(range(len(X)), [self._P[x] for x in X])), X)

    @property
    def support(self) -> Set[T0]:
        """Return the support of the distribution (set of all possible outcomes)."""
        return set(self._P.keys())

    def __eq__(self, other: object) -> bool:
        """Compare distributions for equality."""
        if not isinstance(other, ProbabilityDistribution):
            return NotImplemented
        keys0 = frozenset(self._P.keys())
        keys1 = frozenset(other._P.keys())
        if keys0 != keys1:
            return False
        return all(np.isclose(self._P[x], other._P[x]) for x in keys0)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self._P)})"

    def __getitem__(self, x: T0) -> float:
        """Get the probability associated with a possible outcome."""
        return self._P.get(x, 0.0)

    @classmethod
    def from_empirical_distribution(
        cls, ed: EmpiricalDistribution[T0]
    ) -> "ProbabilityDistribution[T0]":
        """Estimate a probability distribution from an empirical distribution."""
        S = ed.total
        if S == 0:
            raise ValueError("Empirical distribution has no values")
        f = 1 / S
        return cls({x: f * c for x, c in ed.as_counter().items()})

    def condition(
        self, criterion: Callable[[T0], bool]
    ) -> "ProbabilityDistribution[T0]":
        """Return a new distribution conditioned on the given criterion.

        :param criterion: A boolean function defined on all possible outcomes.
        """
        S = sum(c for x, c in self._P.items() if criterion(x))
        if np.isclose(S, 0):
            raise ValueError("Condition has probability zero")
        f = 1 / S
        return ProbabilityDistribution(
            {x: f * c for x, c in self._P.items() if criterion(x)}
        )

    def map(self, mapping: Callable[[T0], T1]) -> "ProbabilityDistribution[T1]":
        """Return a distribution over a transformed domain.

        The provided function maps elements in the original domain to new elements. If
        it is not injective, probabilities are combined.

        :param mapping: A function defined on all possible outcomes, mapping them to
            another domain.
        """
        P: DefaultDict[Any, float] = defaultdict(float)
        for x, p in self._P.items():
            P[mapping(x)] += p
        return ProbabilityDistribution(P)

    def expectation(self, f: Callable[[T0], Number]) -> Number:
        """Compute the expectation value of a functional.

        The provided function maps possible outcomes to numerical values.

        :return: Expectation of the functional.
        """
        return sum(p * f(x) for x, p in self._P.items())

    def variance(self, f: Callable[[T0], Number]) -> Number:
        """Compute the variance of a functional.

        The provided function maps possible outcomes to numerical values.

        :return: Variance of the functional.
        """
        fs = [(f(x), p) for x, p in self._P.items()]
        return sum(p * v**2 for v, p in fs) - (sum(p * v for v, p in fs)) ** 2


def convex_combination(
    dists: List[Tuple[ProbabilityDistribution[T0], float]]
) -> ProbabilityDistribution[T0]:
    """Return a convex combination of probability distributions.

    Each pair in the list comprises a distribution and a weight. The weights must be
    non-negative and sum to 1.

    >>> dist1 = ProbabilityDistribution({0: 0.25, 1: 0.5, 2: 0.25})
    >>> dist2 = ProbabilityDistribution({0: 0.5, 1: 0.5})
    >>> dist3 = convex_combination([(dist1, 0.25), (dist2, 0.75)])
    >>> dist3
    ProbabilityDistribution({0: 0.4375, 1: 0.5, 2: 0.0625})
    >>> dist3.expectation(lambda x : x**2)
    0.75
    """
    P: DefaultDict[T0, float] = defaultdict(float)
    S = 0.0
    for pd, a in dists:
        if a < 0:
            raise ValueError("Weights must be non-negative.")
        for x, p in pd._P.items():
            P[x] += a * p
        S += a
    if not np.isclose(S, 1):
        raise ValueError("Weights must sum to 1.")
    return ProbabilityDistribution(P)
