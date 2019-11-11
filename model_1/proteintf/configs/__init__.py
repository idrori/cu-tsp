from proteintf.configs.baseline import baseline as _baseline_config
from proteintf.configs.birnn import birnn as _birnn_config
from proteintf.configs.transformer import transformer as _transformer_config

CONFIGS = {
    'baseline': _baseline_config,
    'birnn': _birnn_config,
    'transformer': _transformer_config
}
