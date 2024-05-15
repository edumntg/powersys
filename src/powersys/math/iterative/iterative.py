from dataclasses import dataclass
from simple_parsing.helpers import Serializable
from typing import Optional
from ...models.powersystem import *

@dataclass
class IterativeArgs(Serializable):
    max_iters: Optional[int] = 500
    tol: Optional[float] = 1E-6
    verbose: Optional[bool] = False

class Iterative():
    def __init__(self, model: PowerSystem, args: IterativeArgs):
        self.model = model
        self.args = args
        self.max_iters = args.max_iters
        self.tol = args.tol

    def step(self):
        pass
