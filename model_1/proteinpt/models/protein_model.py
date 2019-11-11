from typing import Dict, Optional

import torch
from torch import nn

from allennlp.data import Vocabulary
from allennlp.models.model import Model
from allennlp.training.metrics import Average
from allennlp.nn import InitializerApplicator, RegularizerApplicator
from allennlp.nn.util import get_text_field_mask, add_positional_features
from allennlp.modules import Seq2SeqEncoder, TextFieldEmbedder

from proteinpt.modules import MSAEncoder, FeedForwardWithResNLN, PairwiseBilinearDecoder, ElementwiseLinearDecoder


@Model.register('protein_model')
class ProteinModel(Model):
    def __init__(self,
                 vocab: Vocabulary,
                 target: str,
                 aa_embedder: TextFieldEmbedder,
                 ss_embedder: TextFieldEmbedder,
                 encoder: Seq2SeqEncoder,
                 feedforward: FeedForwardWithResNLN,
                 msa_encoder: MSAEncoder = None,
                 use_ss: bool = True,
                 use_positional_encoding: bool = False,
                 input_dropout: float = 0.,
                 dim_pssm: int = 21,
                 initializer: InitializerApplicator = InitializerApplicator(),
                 regularizer: Optional[RegularizerApplicator] = None) -> None:
        super(ProteinModel, self).__init__(vocab, regularizer)

        self.aa_embedder = aa_embedder
        self.ss_embedder = ss_embedder
        self.pssm_embedder = nn.Linear(dim_pssm, aa_embedder.get_output_dim())

        self.use_ss = use_ss
        self.use_positional_encoding = use_positional_encoding
        self.input_dropout = nn.Dropout(input_dropout)

        self.encoder = encoder
        self.msa_encoder = msa_encoder
        self.feedforward = feedforward

        self.target = target
        if target == 'dcalpha':
            self.decoder = PairwiseBilinearDecoder(feedforward.get_output_dim())
        elif target == 'angles':
            self.decoder = ElementwiseLinearDecoder(feedforward.get_output_dim(), 2)
        elif target == 'coords':
            self.decoder = ElementwiseLinearDecoder(feedforward.get_output_dim(), 3)
        else:
            raise ValueError('Target must be one of dcalpha, angles, or coords, but got {}'.format(target))

        self.metrics = {'rmse': Average()}

        initializer(self)

    def forward(self, **inputs):
        x = self.aa_embedder(inputs['aa']) + self.pssm_embedder(inputs['pssm'])
        if 'ss' in inputs and self.use_ss:
            x += self.ss_embedder(inputs['ss'])
        if self.use_positional_encoding:
            add_positional_features(x)
        x = self.input_dropout(x)

        mask = get_text_field_mask(inputs['aa'])
        x = self.encoder(x, mask)

        if 'msa' in inputs and self.msa_encoder is not None:
            x += self.msa_encoder(inputs['msa'])

        x = self.decoder(self.feedforward(x))

        outputs = {'predictions': x, 'protein_id': inputs['protein_id'], 'length': inputs['length']}
        if 'dcalpha' in inputs:
            mask = mask.unsqueeze(-1).float()
            if self.target == 'dcalpha':
                target = inputs['dcalpha']
                mask = mask.matmul(mask.transpose(dim0=-2, dim1=-1))
                mask_triu = torch.triu(torch.ones_like(mask[0]), diagonal=1).unsqueeze(0).to(mask.device)
                mask *= mask_triu
            elif self.target == 'angles':
                target = torch.stack([inputs['psi'], inputs['phi']], dim=2)
            else:
                target = inputs['coords']

            mse = ((mask * (x - target)) ** 2).sum() / mask.sum()

            if torch.isnan(mse):
                while True:
                    expr = input('\nInput = ')
                    if expr == 'q':
                        exit(0)
                    try:
                        print(eval(expr))
                    except Exception as e:
                        print(e)

            self.metrics['rmse']((mse ** 0.5).detach().item())
            outputs['loss'] = mse

        return outputs

    def get_metrics(self, reset: bool = False) -> Dict[str, float]:
        return {metric_name: metric.get_metric(reset) for metric_name, metric in self.metrics.items()}
