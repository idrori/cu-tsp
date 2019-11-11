import torch
from torch import nn
from typing import Sequence, Union

from allennlp.common import Params
from allennlp.nn import Activation
from allennlp.modules import FeedForward


class FeedForwardWithResNLN(FeedForward):
    def __init__(self,
                 input_dim: int,
                 num_layers: int,
                 hidden_dims: Union[int, Sequence[int]],
                 activations: Union[Activation, Sequence[Activation]],
                 dropout: Union[float, Sequence[float]] = 0.0,
                 residual: bool = False,
                 layer_norm: bool = False) -> None:
        if not isinstance(hidden_dims, list):
            hidden_dims = [hidden_dims] * num_layers  # type: ignore
        super(FeedForwardWithResNLN, self).__init__(input_dim, num_layers, hidden_dims, activations, dropout)
        self._residual = residual
        if residual:
            assert len(set(hidden_dims)) == 1 and hidden_dims[0] == input_dim
        if layer_norm:
            self._layer_norms = nn.ModuleList([nn.LayerNorm(dim) for dim in [input_dim] + hidden_dims[:-1]])
        else:
            self._layer_norms = [lambda x: x for _ in range(len(hidden_dims))]

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        # pylint: disable=arguments-differ
        for layer_norm, layer, activation, dropout in zip(self._layer_norms,
                                                          self._linear_layers,
                                                          self._activations,
                                                          self._dropout):
            outputs = dropout(activation(layer(layer_norm(inputs))))
            if self._residual:
                outputs += inputs
            inputs = outputs
        return inputs

    # Requires custom logic around the activations (the automatic `from_params`
    # method can't currently instantiate types like `Union[Activation, List[Activation]]`)
    @classmethod
    def from_params(cls, params: Params):
        input_dim = params.pop_int('input_dim')
        num_layers = params.pop_int('num_layers')
        hidden_dims = params.pop('hidden_dims')
        activations = params.pop('activations')
        dropout = params.pop('dropout', 0.0)
        residual = params.pop('residual', False)
        layer_norm = params.pop('layer_norm', False)
        if isinstance(activations, list):
            activations = [Activation.by_name(name)() for name in activations]
        else:
            activations = Activation.by_name(activations)()
        params.assert_empty(cls.__name__)
        return cls(input_dim=input_dim,
                   num_layers=num_layers,
                   hidden_dims=hidden_dims,
                   activations=activations,
                   dropout=dropout,
                   residual=residual,
                   layer_norm=layer_norm)
