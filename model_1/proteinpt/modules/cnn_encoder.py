import torch.nn as nn
from allennlp.common import Params


def ensure_is_list(value, list_len):
    if not isinstance(value, list):
        value = [value] * list_len
    return value


class CNNEncoder(nn.Module):
    def __init__(self, output_channels, conv_kernel_sizes,
                 conv_strides=1, conv_paddings=0, pool_kernel_sizes=0, pool_strides=None, pool_paddings=0):
        super(CNNEncoder, self).__init__()

        list_len = len(output_channels)
        conv_strides = ensure_is_list(conv_strides, list_len)
        conv_paddings = ensure_is_list(conv_paddings, list_len)
        pool_kernel_sizes = ensure_is_list(pool_kernel_sizes, list_len)
        pool_strides = ensure_is_list(pool_strides, list_len)
        pool_paddings = ensure_is_list(pool_paddings, list_len)
        assert len(output_channels) == \
               len(conv_kernel_sizes) == len(conv_strides) == len(conv_paddings) == \
               len(pool_kernel_sizes) == len(pool_strides) == len(pool_paddings)

        ic = 2
        conv_layers = []
        for oc, cks, cs, cp, pks, ps, pp in zip(output_channels,
                                                conv_kernel_sizes, conv_strides, conv_paddings,
                                                pool_kernel_sizes, pool_strides, pool_paddings):
            conv_layers += [nn.Conv2d(in_channels=ic,
                                      out_channels=oc,
                                      kernel_size=cks,
                                      stride=cs,
                                      padding=cp),
                            nn.BatchNorm2d(oc),
                            nn.ReLU()]
            ic = oc
            if pks > 0:
                conv_layers += [nn.MaxPool2d(kernel_size=pks,
                                             stride=ps,
                                             padding=pp)]
        conv_layers += [nn.AdaptiveAvgPool2d((None, 1))]
        self.conv = nn.Sequential(*conv_layers)
        self.output_dim = ic

        self.output_channels = output_channels
        self.conv_kernel_sizes = conv_kernel_sizes
        self.conv_strides = conv_strides
        self.conv_paddings = conv_paddings
        self.pool_kernel_sizes = pool_kernel_sizes
        self.pool_strides = pool_strides
        self.pool_paddings = pool_paddings

    def get_output_dim(self):
        return self.output_dim

    def forward(self, inputs, labels=None):
        return self.conv(inputs).squeeze(dim=3).transpose(dim0=1, dim1=2)

    @classmethod
    def from_params(cls, params: Params):
        output_channels = params.pop('output_channels')
        conv_kernel_sizes = params.pop('conv_kernel_sizes')
        conv_strides = params.pop('conv_strides', 1)
        conv_paddings = params.pop('conv_paddings', 0)
        pool_kernel_sizes = params.pop('pool_kernel_sizes', 0)
        pool_strides = params.pop('pool_strides', None)
        pool_paddings = params.pop('pool_paddings', 0)
        params.assert_empty(cls.__name__)
        return cls(output_channels=output_channels,
                   conv_kernel_sizes=conv_kernel_sizes,
                   conv_strides=conv_strides,
                   conv_paddings=conv_paddings,
                   pool_kernel_sizes=pool_kernel_sizes,
                   pool_strides=pool_strides,
                   pool_paddings=pool_paddings)
