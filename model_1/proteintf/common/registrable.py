"""
Code taken from https://github.com/allenai/allennlp/blob/master/allennlp/common/registrable.py
"""

from collections import defaultdict
from typing import TypeVar, Type, Dict, List

T = TypeVar('T')


class Registrable:
    _registry: Dict[Type, Dict[str, Type]] = defaultdict(dict)
    default_implementation: str = None

    @classmethod
    def register(cls: Type[T], name: str):
        registry = Registrable._registry[cls]

        def add_subclass_to_registry(subclass: Type[T]):
            # Add to registry, raise an error if key has already been used.
            if name in registry:
                message = "Cannot register %s as %s; name already in use for %s" % (
                        name, cls.__name__, registry[name].__name__)
                raise ValueError(message)
            registry[name] = subclass
            return subclass
        return add_subclass_to_registry

    @classmethod
    def by_name(cls: Type[T], name: str) -> Type[T]:
        # print(f"instantiating registered subclass {name} of {cls}")
        if name not in Registrable._registry[cls]:
            raise KeyError("%s is not a registered name for %s" % (name, cls.__name__))
        return Registrable._registry[cls].get(name)

    @classmethod
    def list_available(cls) -> List[str]:
        """List default first if it exists"""
        keys = list(Registrable._registry[cls].keys())
        default = cls.default_implementation

        if default is None:
            return keys
        elif default not in keys:
            message = "Default implementation %s is not registered" % default
            raise KeyError(message)
        else:
            return [default] + [k for k in keys if k != default]

    @classmethod
    def from_params(cls, params):
        pass
