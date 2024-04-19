from dataclasses import dataclass
from simple_parsing.helpers import Serializable

@dataclass
class IterativeArgs(Serializable):
    max_iters: int
    tol: float

class Iterative():
    def __init__(self, args: IterativeArgs):
        self.args = args
        self.max_iters = args.max_iters
        self.tol = args.tol

    def step(self):
        pass
