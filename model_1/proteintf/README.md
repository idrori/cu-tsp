Model 1: BiGRUs with skip connections and embeddings

| Author     | Affiliation         | Email               |
| ---------- | ------------------- | ------------------- |
| Jinhao Lei | Columbia University | jl5207@columbia.edu |
| Weiyi Lu | Columbia University | wl2681@columbia.edu |

All the commands below should be executed in the root directory `model_1/`. Assume all the train and test pickle files are under `../data/35k/` or `../data/100k/`

To train the model that predicts dcalphas,
```bash
python -m proteintf.trainer birnn \
	-s <model_path> \
	-o '{"targets": ["dcalphas"], "dataset_reader": {"root": "../data/35k"}}'
```
  
where `<model_path>` is the path to the folder where the model will be saved.

To make predictions using the trained model,
```bash
python -m proteintf.predictor <model_path>
```

where `<model_path>` should be the same as the one used above. The predictions will be saved at `<model_path>/predictions/dcalphas.pkl`.

Similarly, to train the model that predicts psis and phis,
```bash
python -m proteintf.trainer birnn \
	-s ../outputs/psis_phis \
	-o '{"targets": ["psis", "phis"], "dataset_reader": {"root": ""../data/35k""}}'
```
  
And to make predictions using the trained model,
```bash
python -m proteintf.predictor <model_path>
```

The predictions will be saved at `<model_path>/predictions/{psis,phis}.pkl`.

Finally to combine all the outputs into one pickle file, run
```bash
python -m proteintf.utils.make_output \
	--testpkl-path <path to test.pkl> \
	--dcalphas-path <path to dcalphas.pkl> \
	--psis-path <path to psis.pkl> \
	--phis-path <path to phis.pkl> \
	--output-dir <output directory path>
```
The output will be put under `<output directory path>` as `output.pkl`.

	

