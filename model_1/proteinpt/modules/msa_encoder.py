import torch.nn as nn
from allennlp.common import Params
from allennlp.modules import FeedForward


class MSAEncoder(nn.Module):
    def __init__(self,
                 feedforward: FeedForward,
                 output_dim: int,
                 layer_norm: bool = True):
        super(MSAEncoder, self).__init__()
        self.feedforward = feedforward
        self.output_dim = output_dim
        self.output_layer = nn.Linear(feedforward.get_output_dim(), output_dim)
        self.layer_norm = lambda x: x
        if layer_norm:
            self.layer_norm = nn.LayerNorm(output_dim)

    def get_output_dim(self):
        return self.output_dim

    def forward(self, msa):
        msa = self.feedforward(msa).mean(dim=1)
        return self.layer_norm(self.output_layer(msa))

    @classmethod
    def from_params(cls, params: Params):
        feedforward = FeedForward.from_params(params.pop('feedforward'))
        output_dim = params.pop_int('output_dim')
        layer_norm = params.pop('layer_norm', True)
        params.assert_empty(cls.__name__)
        return cls(feedforward=feedforward, output_dim=output_dim, layer_norm=layer_norm)
