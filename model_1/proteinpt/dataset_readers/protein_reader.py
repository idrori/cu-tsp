import pickle
import logging
import numpy as np
from typing import Any, Dict
from overrides import overrides

from allennlp.data.tokenizers import Token
from allennlp.data.instance import Instance
from allennlp.data.dataset_readers import DatasetReader
from allennlp.data.fields import TextField, ArrayField, MetadataField
from allennlp.data.token_indexers import TokenIndexer, SingleIdTokenIndexer

logger = logging.getLogger(__name__)


@DatasetReader.register('protein_reader')
class ProteinDatasetReader(DatasetReader):
    def __init__(self,
                 lazy: bool = False,
                 aa_indexers: Dict[str, TokenIndexer] = None,
                 ss_indexers: Dict[str, TokenIndexer] = None) -> None:
        super(ProteinDatasetReader, self).__init__(lazy)
        self._aa_indexers = aa_indexers or {'aa': SingleIdTokenIndexer(namespace='aa')}
        self._ss_indexers = ss_indexers or {'ss': SingleIdTokenIndexer(namespace='ss')}

    @overrides
    def _read(self, file_path):
        with open(file_path, 'r') as fin:
            for f_name in fin:
                with open(f_name.rstrip('\n'), 'rb') as fin:
                    data = pickle.load(fin)
                    for protein_id, values in data.items():
                        yield self.text_to_instance(protein_id, values)

    @overrides
    def text_to_instance(self,
                         protein_id: str,
                         values: Dict[str, Any]) -> Instance:
        instance = {'protein_id': MetadataField(protein_id), 'length': MetadataField(values['length'])}

        aa = [Token(l) for l in values['aa']]
        instance['aa'] = TextField(aa, self._aa_indexers)

        if len(values.get('ss', [])) > 0:
            ss = [Token(l) for l in values['ss']]
            instance['ss'] = TextField(ss, self._ss_indexers)

        instance['pssm'] = ArrayField(np.array(values['pssm'], dtype=np.float32).transpose())
        if 'msacovonehot' in values:
            msacovonehot = values['msacovonehot']
            msacovembedding = values['msacovembedding']
            if not isinstance(msacovonehot, np.ndarray):
                msacovonehot = msacovembedding = np.array([[0.]])
            msa = np.stack([msacovonehot.astype(np.float32), msacovembedding.astype(np.float32)], axis=2)
            instance['msa'] = ArrayField(msa)

        if 'dcalpha' in values:
            for field in ['dcalpha', 'phi', 'psi']:
                target = np.array(values[field], dtype=np.float32)
                target[np.isnan(target)] = 0.
                instance[field] = ArrayField(target)
            instance['coords'] = ArrayField(np.stack(values['coords']).astype(np.float32))

        return Instance(instance)
