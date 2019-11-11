| Model \ Dataset                                    | 100k without q8 | 100k | 100k with MSA |
| -------------------------------------------------- | ----- | ----- | ----- |
| Model 1: BiGRUs with skip connections and embeddings |  |   |   |

| Author     | Affiliation         | Email               |
| ---------- | ------------------- | ------------------- |
| Jinhao Lei | Columbia University | jl5207@columbia.edu |
| Weiyi Lu | Columbia University | wl2681@columbia.edu |

The model is written using [allennlp](https://allennlp.org/). All the commands below should be executed in the root directory `model_1/`.

Assuming the train/validation/test data are in `../data/100k/{train,validation,test}`, run the following command to create the input files
```bash
for dat in {train,validation,test}
do
find ../data/100k/$dat -name "*pkl" > 100k_$dat.txt
done
```

We train separate models to predict the distance matrix, the torsion angles, and the 3D coordinates. Use the following command to train a model that predicts one of the targets
```bash
allennlp train proteinpt/configs/100k/bigru.json \
	-o '{"train_data_path": "100k_train.txt", "validation_data_path": "100k_validation.txt", "model": {"target": "$TARGET"}}' \
	-f -s outputs/100k_bigru_$TARGET \
	--include proteinpt
```
where `$TARGET` should be replaced by `dcalpha`, `angles`, or `coords`, corresponding to the distance matrix, the torsion angles, and the 3D coordinates, respectively. The trained model will be saved in `outputs/100k_bigru_$TARGET`.
 
Use the following command to make predictions using a trained model,
```bash
python -m proteinpt.utils.predictor \
	--model-path outputs/100k_bigru_$TARGET/model.tar.gz \
	--input-path 100k_test.txt \
	--output-path outputs/100k_bigru_$TARGET/predictions.pkl
```
where again `$TARGET` is one of the three targets. The predictions are saved in `outputs/100k_bigru_$TARGET/predictions.pkl`.

Finally to combine all the predictions into one pickle file, run
```bash
python -m proteinpt.utils.make_output \
	--dcalphas-path outputs/100k_bigru_dcalpha/predictions.pkl \
	--angles-path outputs/100k_bigru_angles/predictions.pkl \
	--coords-path outputs/100k_bigru_coords/predictions.pkl \
	--output-path outputs/100k_bigru_predictions.pkl
```
The outputs are saved in `outputs/100k_bigru_predictions.pkl`.
