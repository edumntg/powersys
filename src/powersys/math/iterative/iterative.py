from dataclasses import dataclass
from simple_parsing.helpers import Serializable
from typing import Optional

@dataclass
class IterativeArgs(Serializable):
    max_iters: Optional[int] = 500
    tol: Optional[float] = 1E-8

class Iterative():
    def __init__(self, args: IterativeArgs):
        self.args = args
        self.max_iters = args.max_iters
        self.tol = args.tol

    def step(self):
        pass
