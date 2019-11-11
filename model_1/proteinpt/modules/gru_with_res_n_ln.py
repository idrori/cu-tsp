import torch
from torch import nn

from allennlp.modules import Seq2SeqEncoder
from allennlp.modules.seq2seq_encoders import PytorchSeq2SeqWrapper


@Seq2SeqEncoder.register('gru_with_res_n_ln')
class GRUWithResNLN(Seq2SeqEncoder):
    def __init__(self,
                 input_size: int,
                 hidden_size: int,
                 num_layers: int,
                 bias: bool = True,
                 bidirectional: bool = False,
                 residual: bool = False,
                 layer_norm: bool = False):
        super(GRUWithResNLN, self).__init__()

        self.input_size = input_size
        self.output_size = hidden_size
        self.bidirectional = bidirectional
        if bidirectional:
            self.output_size *= 2

        self.residual = residual
        if residual:
            assert self.input_size == self.output_size

        if layer_norm:
            self.layer_norms = nn.ModuleList([
                nn.LayerNorm(dim) for dim in [self.input_size] + [self.output_size] * (num_layers - 1)
            ])
        else:
            self.layer_norms = [lambda x: x for _ in range(len(num_layers))]

        self.gru_layers = nn.ModuleList([PytorchSeq2SeqWrapper(
            nn.modules.GRU(input_size=input_size,
                           hidden_size=hidden_size,
                           num_layers=1,
                           bias=bias,
                           batch_first=True,
                           dropout=0.,
                           bidirectional=bidirectional)
        ) for _ in range(num_layers)])

    def get_input_dim(self) -> int:
        return self.input_size

    def get_output_dim(self) -> int:
        return self.output_size

    def is_bidirectional(self) -> bool:
        raise self.bidirectional

    def forward(self,
                inputs: torch.Tensor,
                mask: torch.Tensor,
                hidden_state: torch.Tensor = None) -> torch.Tensor:
        for layer_norm, gru_layer in zip(self.layer_norms, self.gru_layers):
            outputs = gru_layer(layer_norm(inputs), mask, hidden_state)
            if self.residual:
                outputs += inputs
            inputs = outputs
        return inputs
