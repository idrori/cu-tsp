from torch import nn


class PairwiseBilinearDecoder(nn.Module):
    def __init__(self, dim_in):
        super(PairwiseBilinearDecoder, self).__init__()
        self.weight = nn.Linear(dim_in, dim_in)

    def forward(self, x):
        return x.matmul(self.weight(x).transpose(dim0=-2, dim1=-1))


class ElementwiseLinearDecoder(nn.Module):
    def __init__(self, dim_in, dim_out):
        super(ElementwiseLinearDecoder, self).__init__()
        self.weight = nn.Linear(dim_in, dim_out)

    def forward(self, x):
        return self.weight(x)
