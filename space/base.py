# Libraries
from abc import ABC
# Scripts
from ..types import *


class Space(ABC):
    DIMALTY : NumericIntegerValueType

    def __contains__(self, item):
        # For class-level checking
        if issubclass(item, self):
            return True
        # Class attributes
        if hasattr(item, "HOST_SPCE"):
            if item.HOST_SPCE == self:
                return True
            return item.HOST_SPCE in self
        # Instance attributes
        if hasattr(item, "host_spce"):
            if item.host_spce == self:
                return True
            return item.host_spce in self
        # Not found
        return False

class Space0D(Space, ABC):
    DIMALTY = 0
class Space1D(Space, ABC):
    DIMALTY = 1
class Space2D(Space, ABC):
    DIMALTY = 2
class Space3D(Space, ABC):
    DIMALTY = 3

