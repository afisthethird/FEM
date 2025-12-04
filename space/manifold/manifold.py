# Libraries
from abc import ABC, abstractmethod
from dataclasses import dataclass
# Scripts
from types import *


# Goal: order-based node insertion into Topologies
# Goal: abstract/templated assembly of CoordinateSpaces into ManifoldSpaces


class ManifoldSpace(Space, ABC)